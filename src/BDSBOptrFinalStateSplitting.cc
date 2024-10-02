/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2023.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
// this class needs headers from Geant4 10.1 onwards
#include "G4Version.hh"

#if G4VERSION_NUMBER > 1009

#include "BDSBOptrFinalStateSplitting.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSPhysicsUtilities.hh"

#include "globals.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4BiasingProcessSharedData.hh"
#include "BDSBOptnFinalStateSplitting.hh"
#include "G4InteractionLawPhysical.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Proton.hh"
#include "G4VProcess.hh"

#include <limits>
#include <regex>

BDSBOptrFinalStateSplitting::BDSBOptrFinalStateSplitting(const G4String& particleNameIn,
						       const G4String& name):
  G4VBiasingOperator(name),
  fSetup(true),
  particleName(particleNameIn),
  particleIsIon(false)
{
  fParticleToSplit = G4ParticleTable::GetParticleTable()->FindParticle(particleName);

  if (!fParticleToSplit)
    {throw BDSException(__METHOD_NAME__, "Particle \"" + particleName + "\" not found");}
  
  particleIsIon = BDS::IsIon(fParticleToSplit);
}

BDSBOptrFinalStateSplitting::~BDSBOptrFinalStateSplitting()
{
  for (auto change : fChangeFinalStateOperations)
    {delete change.second;}
}

void BDSBOptrFinalStateSplitting::StartRun()
{
  // Setup stage:
  // Start by collecting processes under biasing, create needed biasing
  // operations and associate these operations to the processes:
  if (fSetup)
    {
      const G4ProcessManager*           processManager = fParticleToSplit->GetProcessManager();
      const G4BiasingProcessSharedData* sharedData     = G4BiasingProcessInterface::GetSharedData(processManager);
      
      if (sharedData)
	{
	  // sharedData tested, as is can happen a user attaches an operator to a
	  // volume but without defined BiasingProcessInterface processes.
    for (const auto& wrapperProcess : sharedData->GetPhysicsBiasingProcessInterfaces())
      {
	      G4String operationName = "FSchange-"+wrapperProcess->GetWrappedProcess()->GetProcessName();
	      fChangeFinalStateOperations[wrapperProcess] = new BDSBOptnFinalStateSplitting(operationName);
	      fSplittingFactor[wrapperProcess]    = 1;
	      fSplittingThreshold[wrapperProcess] = 0.0;
	    }
	}
      fSetup = false;
    }
}

void BDSBOptrFinalStateSplitting::SetBias(const G4String& processName,
                                          G4int    factor,
                                          G4double threshold)
{
  const G4ProcessManager*           processManager = fParticleToSplit->GetProcessManager();
  const G4BiasingProcessSharedData* sharedData     = G4BiasingProcessInterface::GetSharedData(processManager);

  G4bool allProcesses = false;
  if (processName == "all")
    {allProcesses = true;}
  
  G4bool processFound = false;
  for (const auto& wrapperProcess : sharedData->GetPhysicsBiasingProcessInterfaces())
    {
      G4String currentProcess = wrapperProcess->GetWrappedProcess()->GetProcessName();
      
      // check if the name is already wrapped for biasing of some kind or splitting
      std::regex braces("[\\w\\-\\+_$]*\\((\\w+)\\)");
      //std::regex braces("[\\w\\-\\_\\+]*\\((\\w+)\\)");
      std::smatch match;
      if (std::regex_search(currentProcess, match, braces))
        {currentProcess = match[1];} // overwrite the variable to match (in this scope)
      
      if (allProcesses || processName == currentProcess)
	{
      fSplittingFactor[wrapperProcess]    = factor;
      fSplittingThreshold[wrapperProcess] = threshold;
	  processFound                  = true; // the process was found at some point

      fChangeFinalStateOperations[wrapperProcess]->SetSplittingFactor(factor);
      fChangeFinalStateOperations[wrapperProcess]->SetSplittingThreshold(threshold);
	}
    }
  if (!processFound)
    {
      G4cout << "\nCouldn't find process by name. Available processes are:" << G4endl;
      for (const auto wrapperProcess : sharedData->GetPhysicsBiasingProcessInterfaces())
	{
	  G4String currentProcessName = wrapperProcess->GetWrappedProcess()->GetProcessName();
	  G4cout << "\"" << currentProcessName << "\"" << G4endl;
	}
      throw BDSException(__METHOD_NAME__, "Process \"" + processName +
			 "\" not found registered to particle \"" + particleName + "\"");
    }
}

G4VBiasingOperation* BDSBOptrFinalStateSplitting::ProposeFinalStateBiasingOperation(const G4Track*                   track,
										  const G4BiasingProcessInterface* callingProcess)
{
  // Check if current particle type is the one to bias:
  const G4ParticleDefinition* definition = track->GetDefinition();
  if (particleIsIon)
    {// we're looking for an ion and this generally isn't an ion
      if (!G4IonTable::IsIon(definition) || definition == G4Proton::Definition())
        {return nullptr;}
      // else it's generally an ion so continue - ie apply GenericIon to any ion
    }
  else if (definition != fParticleToSplit)
    {return nullptr;}
  
  // fetch the operation associated to this callingProcess:
  BDSBOptnFinalStateSplitting* operation = fChangeFinalStateOperations[callingProcess];
  
  // get the operation that was proposed to the process in the previous step:
  G4VBiasingOperation* previousOperation = callingProcess->GetPreviousOccurenceBiasingOperation();

  if (previousOperation != operation)
	{// should not happen !
	  //G4cout << __METHOD_NAME__ << "Logic Problem" << G4endl;
      return nullptr;
    }

  return operation;  
}

void BDSBOptrFinalStateSplitting::OperationApplied(const G4BiasingProcessInterface* callingProcess,
						  G4BiasingAppliedCase,
						  G4VBiasingOperation* occurenceOperationApplied,
						  G4double,
						  G4VBiasingOperation*,    
						  const G4VParticleChange*)
{
  BDSBOptnFinalStateSplitting* operation = fChangeFinalStateOperations[callingProcess];
  if (operation == occurenceOperationApplied)
    {operation->SetInteractionOccured();}
}

#endif

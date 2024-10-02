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
#include "G4Version.hh"
#if G4VERSION_NUMBER > 1009 // consistent with BDSBOptChangeCrossSection

#include "BDSBOptrFinalStateSplitting.hh"
#include "BDSBOptrMultiParticleFinalStateSplitting.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "G4BiasingProcessInterface.hh"

#include "globals.hh"
#include "G4GenericIon.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Proton.hh"


BDSBOptrMultiParticleFinalStateSplitting::BDSBOptrMultiParticleFinalStateSplitting():
  G4VBiasingOperator("BDSIM Final State Splitting"),
  fCurrentOperator(nullptr),
  fnInteractions(0)
{
#ifdef BDSDEBUG
  debug = true;
#else
  debug = false;
#endif
}

BDSBOptrMultiParticleFinalStateSplitting::~BDSBOptrMultiParticleFinalStateSplitting()
{}

void BDSBOptrMultiParticleFinalStateSplitting::AddParticle(const G4String& particleName)
{
  const G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  
  if (!particle) 
    {throw BDSException(__METHOD_NAME__, "Particle \"" + particleName + "\" not found");}

  BDSBOptrFinalStateSplitting* optr = new BDSBOptrFinalStateSplitting(particleName,particleName);
  fParticlesToBias.push_back(particle);
  fBOptrForParticle[particle] = optr;
  optr->StartRun();
}

void BDSBOptrMultiParticleFinalStateSplitting::SetBias(const G4String& biasObjectName,
						      const G4String& particleName,
						      const G4String& process,
						      G4double        dBias,
                              G4double        threshold,
                              const G4String& product,
                              std::vector<std::string> productList)
{
  // important feedback for the user
  G4cout << "Bias> Biasing process \"" << process << "\" for particle \"" << particleName << "\" by factor " << dBias;
  G4cout << ", producing " << product << " particles" << G4endl;
  
  const G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if (!particle)
    {throw BDSException(__METHOD_NAME__, "Particle \"" + particleName + "\" not found");}
  
  try
    {fBOptrForParticle[particle]->SetBias(process, dBias, threshold);}
  catch (BDSException& e)
    {
      e.AppendToMessage("in bias definition \"" + biasObjectName + "\"");
      throw e;
    }
}

G4VBiasingOperation* BDSBOptrMultiParticleFinalStateSplitting::ProposeFinalStateBiasingOperation(const G4Track* track,
											       const G4BiasingProcessInterface* callingProcess)
{
  // -- examples of limitations imposed to apply the biasing:
  // -- limit application of biasing to primary particles only:
  //  if ( track->GetParentID() != 0 ) return 0;
  // -- limit to at most 5 biased interactions:
  //  if ( fnInteractions > 4 )        return 0;
  // -- and limit to a weight of at least 0.05:
  //  if ( track->GetWeight() < 0.05 ) return 0;
  
  if (fCurrentOperator)
    {return fCurrentOperator->GetProposedOccurenceBiasingOperation(track, callingProcess);}
  else
    {return nullptr;}
}

void BDSBOptrMultiParticleFinalStateSplitting::StartTracking(const G4Track* track)
{
  // -- fetch the underneath biasing operator, if any, for the current particle type:
  const G4ParticleDefinition* definition = track->GetParticleDefinition();
  std::map <const G4ParticleDefinition*,BDSBOptrFinalStateSplitting*>::iterator it = fBOptrForParticle.find(definition);
  fCurrentOperator = 0;
  if (it != fBOptrForParticle.end())
    {fCurrentOperator = (*it).second;}

  // try again for ions as they have a generic and specific definition
  // processes are attached to the generic one
  if (G4IonTable::IsIon(definition) && definition != G4Proton::Definition())
    {
      auto search = fBOptrForParticle.find(G4GenericIon::Definition());
      if (search != fBOptrForParticle.end())
	{fCurrentOperator = search->second;}
    }
  // -- reset count for number of biased interactions:
  fnInteractions = 0;
}

void  BDSBOptrMultiParticleFinalStateSplitting::
OperationApplied(const G4BiasingProcessInterface*               callingProcess, 
		 G4BiasingAppliedCase                              biasingCase,
		 G4VBiasingOperation*                occurenceOperationApplied, 
		 G4double                        weightForOccurenceInteraction,
		 G4VBiasingOperation*               finalStateOperationApplied, 
		 const G4VParticleChange*               particleChangeProduced)
{
  // -- count number of biased interactions:
  fnInteractions++;
  
  // -- inform the underneath biasing operator that a biased interaction occurred:
  if (fCurrentOperator)
    {
      fCurrentOperator->ReportOperationApplied(callingProcess,
					       biasingCase,
					       occurenceOperationApplied,
					       weightForOccurenceInteraction,
					       finalStateOperationApplied,
					       particleChangeProduced);
    }
}

#endif

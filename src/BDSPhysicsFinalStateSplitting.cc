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
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSPhysicsFinalStateSplitting.hh"
#include "BDSWrapperFinalStateSplitting.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4String.hh"
#include "G4TrackFastVector.hh"
#include "G4Types.hh"
#include "G4Version.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include <map>
#include <set>

// constructor : main object/class
BDSPhysicsFinalStateSplitting::BDSPhysicsFinalStateSplitting(G4String splittingParentParticleIn,
                                                             G4String splittingProcessIn,
                                                             G4int splittingFactorIn,
                                                             G4double splittingThresholdEKIn,
                                                             std::vector<std::string> splittingProductParticlesIn):
  G4VPhysicsConstructor("BDSPhysicsFinalStateSplitting"),
  splittingParentParticle(splittingParentParticleIn),
  splittingProcess(splittingProcessIn),
  splittingFactor(splittingFactorIn),
  splittingThresholdEK(splittingThresholdEKIn),
  splittingProductParticles(splittingProductParticlesIn)
{
  // Must have a splitting factor above 1
  if (splittingFactorIn < 1)
    {throw BDSException(__METHOD_NAME__, "the splitting factor must be an integer 1 or greater.");}
  G4int maxSize = G4TrackFastVectorSize/2;
  if (splittingFactorIn > maxSize)
  {
    G4String msg = "the maximum safe splitting factor is " + std::to_string(maxSize);
    msg += " based on the G4TrackFastVectorSize in Geant4";
    throw BDSException(__METHOD_NAME__, msg);
  }
}

// destructor
BDSPhysicsFinalStateSplitting::~BDSPhysicsFinalStateSplitting()
{;}

// Not using it. Stays at default
void BDSPhysicsFinalStateSplitting::ConstructParticle()
{;}

// find given particles and processes happening on said particles
void BDSPhysicsFinalStateSplitting::ConstructProcess()
{
    // Single use check
    if (Activated())
    {return;}

    // define input particles and associated input processes


    // Get particles from G4 table
#if G4VERSION_NUMBER > 1029
    auto aParticleIterator =  G4ParticleTable::GetParticleTable()->GetIterator();
#endif
    aParticleIterator->reset();

    // Get Helper to add and remove processes in PhysicsList
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

    // Loop on G4 particles
    while( (*aParticleIterator)() ) {
      // Get particle information
      G4ParticleDefinition *particle = aParticleIterator->value();
      G4String particleName = particle->GetParticleName();

      // Select the input particle to bias
      if (particleName == splittingParentParticle) {

        // Get process list associated with this particle
        G4ProcessManager *pManager = particle->GetProcessManager();
        G4ProcessVector *processVector = pManager->GetProcessList();

        // Loop on number of G4 processes associated
        for (G4int i = 0; i < (G4int) processVector->entries(); ++i) {
          // Get process information
          G4VProcess *process = (*processVector)[i];
          G4String processName = process->GetProcessName();

          if (processName == splittingProcess) {
            // Wrapper. Do the splitting.
            auto wrappedProcess = new BDSWrapperFinalStateSplitting(process, splittingFactor,
                                                                    splittingThresholdEK, splittingProductParticles);

            // Remove process from list to avoid duplicate data from next line
            pManager->RemoveProcess(process);

            // Add wrapped process to the other processes
            ph->RegisterProcess(wrappedProcess, particle);

            // Prints stuff
            G4cout << "Bias> final state splitting> wrapping \"" << process->GetProcessName()
                   << "\" for particle \"" << particle->GetParticleName() << "\": x" << splittingFactor
                   << " for parent Ek > " << splittingThresholdEK / CLHEP::GeV << " GeV" << G4endl;
          }
        }
      }
    }


    // Single use activation
    SetActivated();
}
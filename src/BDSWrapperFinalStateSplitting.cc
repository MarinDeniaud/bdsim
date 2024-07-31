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
#include "BDSPhysicsVectorLinear.hh"
#include "BDSWrapperFinalStateSplitting.hh"

#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Types.hh"
#include "G4VParticleChange.hh"
#include "G4VProcess.hh"

#include <cmath>
#include <limits>
#include <vector>

// Counter for understanding occurence.
G4int BDSWrapperFinalStateSplitting::nCallsThisEvent = 0;

// constructor : main object/class
BDSWrapperFinalStateSplitting::BDSWrapperFinalStateSplitting(G4VProcess* originalProcess,
                                                             G4int splittingFactorIn,
                                                             G4double splittingThresholdEKIn,
                                                             std::vector<std::string> splittingProductParticlesIn):
  BDSWrapperProcess("FinalStateSplittingWrapper"),
  splittingFactor(splittingFactorIn),
  splittingThresholdEK(splittingThresholdEKIn),
  splittingProductParticles(splittingProductParticlesIn)

{
    // Get infos on the original processes
  RegisterProcess(originalProcess);
  theProcessSubType = originalProcess->GetProcessSubType();
  theProcessName = "FinalStateSplittingWrapper("+originalProcess->GetProcessName()+")";
}

// destructor
BDSWrapperFinalStateSplitting::~BDSWrapperFinalStateSplitting()
{;}

// Perform the splitting
G4VParticleChange* BDSWrapperFinalStateSplitting::PostStepDoIt(const G4Track& track,
							 const G4Step& step)
{
    // G4cout << "BDSWrapperFinalStateSplitting::PostStepDoIt" << G4endl;
    // Get the original particle changes for a given track and step
  G4VParticleChange* particleChange = pRegProcess->PostStepDoIt(track, step);

  // Cases were we return the initial particleChange

  if (splittingFactor == 1)
  {return particleChange;}

  G4double parentEk = track.GetKineticEnergy();
  if (parentEk < 0.8*splittingThresholdEK) // smaller of the two thresholds by design
    {return particleChange;}

  G4int nSecondaries = particleChange->GetNumberOfSecondaries();
  if (nSecondaries == 0)
    {return particleChange;}

  // Test if wanted particle in secondaries. If not, return initial particleChange
  G4bool present = false;
  std::vector<G4int> secondaryPDGIDs;
  for (G4int i = 0; i < nSecondaries; i++)
    {
      G4String secondaryName = particleChange->GetSecondary(i)->GetDefinition()->GetParticleName();
      auto IsParticleToSplit = std::count(splittingProductParticles.begin(), splittingProductParticles.end(), secondaryName) != 0;
      present = IsParticleToSplit || present;
    }
  if (!present)
    {return particleChange;}

  // we keep hold of the tracks and manage their memory
  std::vector<G4Track*> originalSecondaries;
  std::vector<G4Track*> originalWantedParticles;
  for (G4int i = 0; i < nSecondaries; i++)
    {
      G4Track* secondary = particleChange->GetSecondary(i);
      G4String secondaryName = secondary->GetDefinition()->GetParticleName();
      auto IsParticleToSplit = std::count(splittingProductParticles.begin(), splittingProductParticles.end(), secondaryName) != 0;
      if (IsParticleToSplit)
      {originalWantedParticles.push_back(secondary);}
      else
      {originalSecondaries.push_back(secondary);}
    }

  G4int nOriginalSecondaries = nSecondaries;

  particleChange->Clear(); // doesn't delete the secondaries

  // Attempt to generate more of the desired particle.
  // We set a maximum number of trials to get to the given splitting factor
  G4int maxTrials = 10 * splittingFactor;
  G4int nSuccessfulSplits = 0;
  G4int iTry = 0;
  std::vector<G4Track*> newParticles;
  while (iTry < maxTrials && nSuccessfulSplits < splittingFactor-1)
    {
      iTry++;
      particleChange->Clear(); // wipes the vector of tracks, but doesn't delete them
      particleChange = pRegProcess->PostStepDoIt(track, step);
      G4bool success = false;
      for (G4int i = 0; i < particleChange->GetNumberOfSecondaries(); i++)
        {
          G4Track* secondary = particleChange->GetSecondary(i);
          G4String secondaryName = secondary->GetDefinition()->GetParticleName();
          auto IsParticleToSplit = std::count(splittingProductParticles.begin(), splittingProductParticles.end(), secondaryName) != 0;
          if (IsParticleToSplit)
            {
              newParticles.push_back(secondary);
              success = true;
            }
          else
            {delete secondary;}
        }
      particleChange->Clear();
      if (success)
        {nSuccessfulSplits++;}
    }

  // Now we put back the secondaries and the splitted particles in particleChange
  particleChange->Clear();
  particleChange->SetNumberOfSecondaries(nOriginalSecondaries + static_cast<G4int>(newParticles.size()));
  particleChange->SetSecondaryWeightByProcess(true);

  // When no splitted particles: add everything back
  if (nSuccessfulSplits == 0)
    {// we've cleared the original ones, so we have to put them back
      for (auto secondary : originalSecondaries)
        {particleChange->AddSecondary(secondary);}
      for (auto muon : originalWantedParticles)
        {particleChange->AddSecondary(muon);}
      return particleChange;
    }

  // When we do have splitted particles :
  // the original particle(s) count as 1 even if there are 2 of them as it's 1x call to the process
  G4double weightFactor = 1.0 / (static_cast<G4double>(nSuccessfulSplits) + 1.0);
  // put in the original secondaries with an unmodified weight
  for (auto aSecondary : originalSecondaries)
    {particleChange->AddSecondary(aSecondary);}
  for (auto originalWantedParticle : originalWantedParticles)
    {
      G4double existingWeight = originalWantedParticle->GetWeight();
      G4double newWeight = existingWeight * weightFactor;
      originalWantedParticle->SetWeight(newWeight);
      particleChange->AddSecondary(originalWantedParticle);
    }
  for (auto newParticle : newParticles)
    {
      G4double existingWeight = newParticle->GetWeight();
      G4double newWeight = existingWeight * weightFactor;
      newParticle->SetWeight(newWeight);
      particleChange->AddSecondary(newParticle);
    }

  nCallsThisEvent++;
  return particleChange;
}

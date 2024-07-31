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
#ifndef BDSWrapperFinalStateSplitting_H
#define BDSWrapperFinalStateSplitting_H
#include "BDSWrapperProcess.hh"

#include "G4Types.hh"

class BDSPhysicsVectorLinear;
class G4Step;
class G4Track;
class G4VParticleChange;
class G4VProcess;

/**
 * @brief Wrapper process to produce more specific particles by resampling the process.
 *
 * Wrap a process. If that process post-step do-it particle change produces a
 * specified particle, then resample it until we get the desired number of those.
 * Keep only the new particles. Put together the other original secondaries with the
 * now N specified particles and weight each of them by w_i * 1/N.
 *
 * @author Laurie Nevay and Marin Deniaud
 */

class BDSWrapperFinalStateSplitting: public BDSWrapperProcess
{
public:
  BDSWrapperFinalStateSplitting() = delete;
  BDSWrapperFinalStateSplitting(G4VProcess* originalProcess,
                                G4int splittingFactorIn,
                                G4double splittingThresholdEKIn,
                                std::vector<std::string> splittingProductParticlesIn);
  virtual ~BDSWrapperFinalStateSplitting();
  
  /// Do the splitting operation.
  virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
                                          const G4Step& step);
  
  /// Counter for understanding occurence.
  static G4int nCallsThisEvent;
  
private:
    G4int splittingFactor;
    G4double splittingThresholdEK;
    std::vector<std::string> splittingProductParticles;
};

#endif

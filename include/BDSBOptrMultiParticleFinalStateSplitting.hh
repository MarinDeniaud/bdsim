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
#ifndef BDSBOPTRMULTIPARTICLEFINALSTATESPLITTING_H
#define BDSBOPTRMULTIPARTICLEFINALSTATESPLITTING_H

#include "G4Version.hh"
#if G4VERSION_NUMBER > 1009

#include "G4VBiasingOperator.hh"
class BDSBOptrFinalStateSplitting;
class G4ParticleDefinition; 

#include <map>

/**
 * @brief Multi-particle final state splitting.
 *
 * This class uses classes from Geant4 10.0 and is dependent on 
 * BDSBOptFinalStateSplitting which needs 10.1.
 *
 * Largely based on the Geant4 example of this feature.
 *
 * @author Stewart Boogert and Marin Deniaud
 */

class BDSBOptrMultiParticleFinalStateSplitting: public G4VBiasingOperator
{
public:
    BDSBOptrMultiParticleFinalStateSplitting();
  virtual ~BDSBOptrMultiParticleFinalStateSplitting();
  
  void AddParticle(const G4String& particleName);
  void SetBias(const G4String& biasObjectName,
	       const G4String& particleName,
	       const G4String& process,
	       G4double        dBias,
           G4double        threshold,
           const G4String& product,
           std::vector<std::string> productList);
  void StartTracking(const G4Track* track) override;

private: 
  virtual G4VBiasingOperation* ProposeOccurenceBiasingOperation(const G4Track*, const G4BiasingProcessInterface*) override
  {return 0;}

  // -- Method used:
  virtual G4VBiasingOperation* ProposeFinalStateBiasingOperation(const G4Track*                   track,
                                                                 const G4BiasingProcessInterface* callingProcess) override;

  virtual G4VBiasingOperation* ProposeNonPhysicsBiasingOperation(const G4Track*, const G4BiasingProcessInterface*) override
  {return 0;}
  virtual void OperationApplied(const G4BiasingProcessInterface* callingProcess,
				G4BiasingAppliedCase             biasingCase,
 				G4VBiasingOperation*             occurenceOperationApplied,
				G4double                         weightForOccurenceInteraction,
				G4VBiasingOperation*             finalStateOperationApplied, 
				const G4VParticleChange*         particleChangeProduced) override;
  // prevent compiler warning (since second G4VBiasingOperator::OperationApplied is hidden)
  using G4VBiasingOperator::OperationApplied;
  
  std::map<const G4ParticleDefinition*, BDSBOptrFinalStateSplitting*> fBOptrForParticle;
  std::vector<const G4ParticleDefinition*>                            fParticlesToBias;
  BDSBOptrFinalStateSplitting*                                        fCurrentOperator;

  G4int  fnInteractions; ///< Count number of biased interactions for current track.
  G4bool debug;
};

#endif

#endif

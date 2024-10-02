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
#ifndef BDSBOPTNFINALSTATESPLITTING_H
#define BDSBOPTNFINALSTATESPLITTING_H

#include "G4VBiasingOperation.hh"
#include "G4ParticleChange.hh"

class BDSBOptnFinalStateSplitting : public G4VBiasingOperation {
public:
    // -- Constructor :
    BDSBOptnFinalStateSplitting(G4String name);
    // -- destructor:
    virtual ~BDSBOptnFinalStateSplitting();

public:
    // ----------------------------------------------
    // -- Methods from G4VBiasingOperation interface:
    // ----------------------------------------------
    // -- Unused:
    virtual const G4VBiasingInteractionLaw*
    ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*,
                                           G4ForceCondition& )
    { return 0; }

    // --Used:
    virtual G4VParticleChange*   ApplyFinalStateBiasing( const G4BiasingProcessInterface*,
                                                         const G4Track*,
                                                         const G4Step*,
                                                         G4bool&                          );

    // -- Unsued:
    virtual G4double           DistanceToApplyOperation( const G4Track*,
                                                         G4double,
                                                         G4ForceCondition*)
    {return DBL_MAX;}
    virtual G4VParticleChange* GenerateBiasingFinalState( const G4Track*,
                                                          const G4Step*   )
    {return 0;}


public:
    // ----------------------------------------------
    // -- Additional methods, specific to this class:
    // ----------------------------------------------
    // -- Splitting factor:
    void     SetSplittingFactor(G4int splittingFactor)
    { fSplittingFactor = splittingFactor; }
    G4int  GetSplittingFactor() const { return fSplittingFactor; }

    // -- Splitting threshold:
    void     SetSplittingThreshold(G4int splittingThreshold)
    { fSplittingThreshold = splittingThreshold; }
    G4int  GetSplittingThreshold() const { return fSplittingThreshold; }

    G4bool GetInteractionOccured() const { return fInteractionOccured; }
    void   SetInteractionOccured() { fInteractionOccured = true; }

private:
    G4int            fSplittingFactor;
    G4double         fSplittingThreshold;
    G4ParticleChange fParticleChange;
    G4bool           fInteractionOccured;
};

#endif

/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/

#ifndef BDSBeamPipe_h
#define BDSBeamPipe_h 1

#include"globals.hh"
#include "BDSMaterials.hh"
#include "G4LogicalVolume.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4UniformMagField.hh"
#include "G4IntersectionSolid.hh"
#include "G4VSolid.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4FieldManager.hh"

#include "BDSEnergyCounterSD.hh"

class BDSBeamPipe
{
public:
  BDSBeamPipe(const G4String& aName, G4double aLength, G4double aRadius,
             G4double angle=0);
  ~BDSBeamPipe();

  G4LogicalVolume* GetLogicalVolume();
  G4LogicalVolume* GetInnerLogicalVolume();
  G4ThreeVector GetPos();
  G4RotationMatrix* GetRot();

  void SetBPFieldManager(G4FieldManager* aFieldManager);
    void SetCoarseFieldManager(G4FieldManager* aFieldManager);

protected:

private:
  G4LogicalVolume* itsLogicalVolume;
  G4LogicalVolume* itsInnerLogicalVolume;
    G4LogicalVolume* itsCoarseInnerLogicalVolume;
  G4VisAttributes* SetVisAttributes();


  G4UserLimits* itsUserLimits;
  G4VisAttributes* itsVisAttributes;
  
  G4ThreeVector itsPos;
  G4RotationMatrix* itsRot;

  G4Trd* itsTrd1;
  G4Trd* itsTrd2;
  G4IntersectionSolid* itsTubeInTrd;
  G4IntersectionSolid* itsInnerTubeInTrd;

  G4Tubs* itsTube;
  G4Tubs* itsInnerTube;

  BDSEnergyCounterSD* itsECounter;
};

#endif

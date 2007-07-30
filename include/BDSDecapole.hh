/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/

#ifndef BDSDecapole_h
#define BDSDecapole_h 1

#include "globals.hh"
#include "BDSMaterials.hh"
#include "G4LogicalVolume.hh"
#include "BDSDecStepper.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"               

#include "BDSMultipole.hh"
#include "BDSDecMagField.hh"

class BDSDecapole :public BDSMultipole
{
  public:
  BDSDecapole(G4String& aName, G4double aLength,
	      G4double bpRad,G4double FeRad,
	      G4double BQuadPrime, G4String aMaterial = "");
    ~BDSDecapole();

  protected:

  private:
  G4double itsBQuadPrime;

  //  void BuildOuterLogicalVolume();
  void BuildBPFieldAndStepper();
  //void BuildMarkerLogicalVolume();

  G4VisAttributes* SetVisAttributes();

  // field related objects:
  BDSDecStepper* itsStepper;
  BDSDecMagField* itsMagField;
  G4Mag_UsualEqRhs* itsEqRhs;

};

#endif

/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/

#ifndef BDSSampler_h
#define BDSSampler_h 1

#include "globals.hh"
#include "BDSAcceleratorComponent.hh"
#include "BDSMaterials.hh"
#include "G4LogicalVolume.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

class BDSSampler :public BDSAcceleratorComponent
{
public:
  BDSSampler(G4String& aName,G4double aLength);
  ~BDSSampler();

protected:

private:
  void SamplerLogicalVolume();
  G4VisAttributes* SetVisAttributes();

  // field related objects:
  G4UserLimits* itsUserLimits;
  G4VisAttributes* itsVisAttributes;
};

#endif

/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/
//      ------------ BDSAcceleration physics process --------
//                     by Grahame Blair, 18 October 2001
#include "BDSGlobalConstants.hh" // must be first in include list

#include "BDSAcceleration.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

extern G4int event_number;

BDSAcceleration::BDSAcceleration(const G4String& processName)
  : G4VDiscreteProcess(processName)
     // initialization
{} 
 

G4VParticleChange* 
BDSAcceleration::PostStepDoIt(const G4Track& trackData,
				     const G4Step& stepData)
{
  aParticleChange.Initialize(trackData);

G4double z=(trackData.GetPosition()).z()+BDSGlobals->GetWorldSizeZ();

G4double NewKinEnergy;
   NewKinEnergy = 10.*GeV + 490.*GeV*z/(16500.*m);
  //
  // Update the incident particle 
  //
  aParticleChange.SetMomentumChange(trackData.GetMomentumDirection());
  aParticleChange.SetEnergyChange(NewKinEnergy);
  aParticleChange.SetLocalEnergyDeposit (0.); 

  return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
}


BDSAcceleration::~BDSAcceleration(){
}




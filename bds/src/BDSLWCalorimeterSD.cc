/* BDSIM code for LW Calorimeter.    Version 1.0
   Author: John Carter, Royal Holloway, Univ. of London.
   Last modified 26.7.2004
   Copyright (c) 2004 by J.C.Carter.  ALL RIGHTS RESERVED. 
*/
#include "BDSGlobalConstants.hh" // must be first in include list

#include "BDSLWCalorimeterSD.hh"
#include "BDSLWCalorimeterHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4Navigator.hh"
#include "G4AffineTransform.hh"

#include "G4RunManager.hh"
#include <vector.h>

typedef std::vector<G4int> MuonTrackVector;
extern MuonTrackVector* theMuonTrackVector;

extern G4double initial_x,initial_xp,initial_y,initial_yp,initial_z,initial_E;


BDSLWCalorimeterSD::BDSLWCalorimeterSD(G4String name)
:G4VSensitiveDetector(name),StoreHit(true)
{
  collectionName.insert("LWCalorimeterCollection");
}

BDSLWCalorimeterSD::~BDSLWCalorimeterSD()
{;}

void BDSLWCalorimeterSD::Initialize(G4HCofThisEvent*HCE)
{
  LWCalorimeterCollection = 
    new BDSLWCalorimeterHitsCollection(SensitiveDetectorName,collectionName[0]);
  itsTotalEnergy = 0.;
}

G4bool BDSLWCalorimeterSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4Track* theTrack=aStep->GetTrack();
  itsCopyNumber=theTrack->GetVolume()->GetCopyNo()+1;
  AddEnergy(aStep->GetTotalEnergyDeposit());
  //G4cout<<"edep="<<aStep->GetTotalEnergyDeposit()/GeV<<"Total so far="<<itsTotalEnergy/GeV<<G4endl;
  return true;
  
}

void BDSLWCalorimeterSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  G4int nEvent= 
	G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  if(itsCopyNumber > 0){
    BDSLWCalorimeterHit* smpHit=
      new BDSLWCalorimeterHit(itsCopyNumber,itsTotalEnergy,nEvent);
    
    LWCalorimeterCollection->insert(smpHit);
  }
  //if(theTrack->GetVolume()!=theTrack->GetNextVolume())StoreHit=true;
  //else StoreHit=false;
  
  
  static G4int HCID = -1;
  if(HCID<0)
    { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection( HCID, LWCalorimeterCollection );
}

void BDSLWCalorimeterSD::clear(){} 

void BDSLWCalorimeterSD::DrawAll(){} 

void BDSLWCalorimeterSD::PrintAll(){} 


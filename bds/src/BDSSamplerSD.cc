/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/
#include "BDSGlobalConstants.hh" // must be first in include list

#include "BDSSamplerSD.hh"
#include "BDSSamplerHit.hh"
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


BDSSamplerSD::BDSSamplerSD(G4String name)
:G4VSensitiveDetector(name),StoreHit(true)
{
  collectionName.insert("SamplerCollection");
}

BDSSamplerSD::~BDSSamplerSD()
{;}

void BDSSamplerSD::Initialize(G4HCofThisEvent*HCE)
{
  SamplerCollection = 
    new BDSSamplerHitsCollection(SensitiveDetectorName,collectionName[0]); 
}

G4bool BDSSamplerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4Track* theTrack=aStep->GetTrack();

  //  // tmp - only store muons
  //     G4String pName=theTrack->GetDefinition()->GetParticleName();
  //    if(pName=="mu+"||pName=="mu-")
  // 	{ // tmp

  if(StoreHit)
    {
      G4double time=theTrack->GetGlobalTime();

      G4double energy=theTrack->GetKineticEnergy()+
	theTrack->GetDefinition()->GetPDGMass();

      G4ThreeVector pos=theTrack->GetPosition();
      G4ThreeVector momDir=theTrack->GetMomentumDirection();

      const G4RotationMatrix* Rot=theTrack->GetVolume()->GetFrameRotation();
      const G4ThreeVector Trans=theTrack->GetVolume()->GetFrameTranslation();

      G4ThreeVector LocalPosition=pos+Trans; 
      G4ThreeVector LocalDirection=(*Rot)*momDir; 

      G4double x=LocalPosition.x();
      G4double y=LocalPosition.y();
      G4double xPrime=LocalDirection.x();
      G4double yPrime=LocalDirection.y();

      G4double z=-(time*c_light-(pos.z()+BDSGlobals->GetWorldSizeZ()));
      // apply a correction that takes ac... gab to do later!

      G4int nEvent= 
	G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
      G4int nSampler=theTrack->GetVolume()->GetCopyNo()+1;
      G4int PDGtype=theTrack->GetDefinition()->GetPDGEncoding();


      // G4cout<<"particle="<<theTrack->GetDefinition()->GetParticleName()<<
      //	" encoding="<<PDGtype<<G4endl;

      G4double start_z,start_E,start_x,start_xp,start_y,start_yp;
      G4String pName=theTrack->GetDefinition()->GetParticleName();

      G4ThreeVector vtx=theTrack->GetVertexPosition();
      G4ThreeVector dir=theTrack->GetVertexMomentumDirection();


      if(pName!="e+" && pName!="e-") 
	{
	  // store production point
	  start_z =  vtx.z()+BDSGlobals->GetWorldSizeZ();
	  start_E   =  theTrack->GetVertexKineticEnergy();
	  start_x   =  vtx.x();
	  start_xp  =  dir.x();
	  start_y   =  vtx.y();
	  start_yp  =  dir.y();	  
	}
      else
	{// for electrons (and positrons) store the point of last scatter
	  start_z  = initial_z + BDSGlobals->GetWorldSizeZ();
	  start_E  = initial_E;
	  start_x  = initial_x;
	  start_xp = initial_xp;
	  start_y  = initial_y;
	  start_yp = initial_yp;
	}

      G4double weight=theTrack->GetWeight();

      if(BDSGlobals->GetStoreMuonTrajectories())
	if(pName=="mu+"||pName=="mu-") 
	  theMuonTrackVector->push_back(theTrack->GetTrackID());

      /*
      else
	{
	  start_z=initial_z;
	  start_E=initial_E;
	  start_x=initial_x;
	  start_xp=initial_xp;
	  start_y=initial_y;
	  start_yp=initial_yp;
	}
      */

      BDSSamplerHit* smpHit
	= new BDSSamplerHit(
			 nSampler,
			 start_z,start_E,start_x,start_xp,
			 start_y,start_yp,
			 z,energy,x,xPrime,y,yPrime,weight,
			 PDGtype,nEvent);

      SamplerCollection->insert(smpHit);
      if(theTrack->GetVolume()!=theTrack->GetNextVolume())StoreHit=true;
      else StoreHit=false;
      return true;

    }
  else
    {
      if(theTrack->GetVolume()!=theTrack->GetNextVolume())StoreHit=true;
      else StoreHit=false;
      return false;
    }

}

void BDSSamplerSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection( HCID, SamplerCollection );
}

void BDSSamplerSD::clear(){} 

void BDSSamplerSD::DrawAll(){} 

void BDSSamplerSD::PrintAll(){} 


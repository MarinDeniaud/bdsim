/* BDSIM code.    
AWAKE spectrometer.
*/

#include "BDSGlobalConstants.hh" 
#include "BDSAwakeSpectrometer.hh"
#include "BDSMaterials.hh"
#include "BDSSampler.hh"
#include "BDSSamplerSD.hh"
#include "BDSCCDCamera.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"               
#include "G4UserLimits.hh"
#include "G4SubtractionSolid.hh"
#include "BDSDebug.hh"
#include "BDSAwakeMultilayerScreen.hh"
//#include "UltraFresnelLens.hh"
//#include "UltraFresnelLensParameterisation.hh"
#include "G4NystromRK4.hh"

#include "G4Trap.hh"
//#include "BDSOutputBase.hh"
#include "BDSDipoleStepper.hh"
#include "BDSUniformMagField.hh"

typedef std::map<G4String,G4LogicalVolume*> LogVolMap;
extern LogVolMap* LogVol;
#define DEBUG 1
//============================================================
BDSAwakeSpectrometer::BDSAwakeSpectrometer (G4String aName, G4double length=2.7*CLHEP::m, G4String bmapFile = "", G4double BField=0, G4double poleStartZ=62.733*CLHEP::cm, G4String material="lanex", G4double thickness = 0.3 * CLHEP::mm, G4double windowScreenGap=0, G4double angle = -45*CLHEP::pi/180.0, G4double windowThickness=0, G4String windowMaterial="G4_Al", G4double screenEndZ = (258-62.733)*CLHEP::cm, G4String spec=""):
  BDSAcceleratorComponent(aName, length, 0, 0, 0,"","",0,0,0,0,0,0,"",bmapFile), _mlScreen(NULL), _camera(NULL), _BField(BField), _poleStartZ(poleStartZ), _screenEndZ(screenEndZ), _material(material), _thickness(thickness), _screenAngle(angle), _windowScreenGap(windowScreenGap), _windowThickness(windowThickness), _windowMaterial(windowMaterial)
{
  try{
    _vacuumChamberType=getParameterValueInt(spec,"vacuumChamberType");
  } catch(boost::bad_lexical_cast&){
    //If the cast fails, set vacuum chamber type to 1
    _vacuumChamberType=1;
  }

  try{
    _magnetGeometryType=getParameterValueInt(spec,"magnetGeometryType");
  } catch(boost::bad_lexical_cast&){
    //If the cast fails, set magnet geometry type to 1
    _magnetGeometryType=1;
  }


  //Set as part of precision region (for energy loss monitoring)
  itsPrecisionRegion=1;

  //Set the rotation of the screen
  _screenRotationMatrix = new G4RotationMatrix();
    _screenRotationMatrix->rotateY(_screenAngle);

  _vacRotationMatrix = new G4RotationMatrix();

  _magRotationMatrix = new G4RotationMatrix();
}

void BDSAwakeSpectrometer::MagnetDefaults(){
  //Initialise geometry pointers.
  itsCoilLog=NULL;
  //Gap between the pole faces.
  itsPoleAperture=65*CLHEP::mm;
  //Gap between the coil faces.
  itsCoilAperture=180*CLHEP::mm;
  //Part of pole extending below coil
  G4double outerPoleSize=0.5*50*CLHEP::mm;
  //Initialise dimensions.
  itsYokeSize.set(
		  902*CLHEP::mm,	
		  1168*CLHEP::mm,   
		  1*CLHEP::m
		  );
  
  itsCoilSize.set(
		  //  (320+320+262)*CLHEP::mm,
		  262*CLHEP::mm,
		  180*CLHEP::mm,
		  (1000+320+285)*CLHEP::mm
		  );
  
  itsPoleSize.set(
		  320*CLHEP::mm,
		  //itsCoilSize.y()+
		  outerPoleSize,
		  itsYokeSize.z()
		  );
  
  itsYokeMiddleSize=itsYokeSize;
  itsYokeMiddleSize.setY(2*itsCoilSize.y()+itsCoilAperture);
  G4double yokeCoilGap=30*CLHEP::mm;
  itsYokeMiddleSize.setX(itsYokeSize.x()-itsPoleSize.x()-itsCoilSize.x()-yokeCoilGap);
  itsYokeUpperSize=itsYokeSize;
  itsYokeUpperSize.setY((itsYokeSize.y()-itsYokeMiddleSize.y())/2.0);
  itsYokeLowerSize=itsYokeUpperSize;
  
  itsMiddleCoilSize.set(
		  itsPoleSize.x(),
		  itsCoilSize.y(),
		  itsCoilSize.x()
		  );
  



  //Yoke apertures.
  itsAperture1Size.set(
		       262*CLHEP::mm+yokeCoilGap,
		       itsCoilAperture+2*itsCoilSize.y(),
		       itsYokeSize.z()
		       );
  
  itsAperture2Size.set(
		       itsPoleSize.x(),
		       itsPoleAperture,
		       itsYokeSize.z()
		       );

  //the position of the magnet centre relative to the marker volume.
  itsPolePos.set(
		 13*CLHEP::cm,
		 0,
		 _poleStartZ+itsYokeSize.z()/2.0
		 );

  //The field map centre corresponds with the pole centre.
  //  itsBmapZOffset=(-itsLength/2.0 + 62.733*CLHEP::cm)*0.5;
  itsBmapZOffset=0.5*itsPolePos.z();
  itsBmapXOffset=0.5*itsPolePos.x();


  //The position of the yoke relative to the marker volume
  itsYokePos.set(
		 itsPolePos.x()+(itsPoleSize.x()-itsYokeSize.x())/2.0,
		 itsPolePos.y(),
		 itsPolePos.z()
		 );

  itsYokeUpperPos=itsYokePos;
  itsYokeUpperPos.setX(itsPolePos.x()+(itsPoleSize.x()-itsYokeUpperSize.x())/2.0);
  itsYokeUpperPos.setY(itsYokePos.y()+itsYokeMiddleSize.y()/2.0+itsYokeUpperSize.y()/2.0);

  itsYokeMiddlePos=itsYokePos;
  //  itsYokeMiddlePos.setX(-itsYokeSize.x()/2.0+itsYokeMiddleSize.x()/2.0);
  itsYokeMiddlePos.setX(itsYokePos.x()-itsYokeSize.x()/2.0+itsYokeMiddleSize.x()/2.0);

  itsYokeLowerPos=itsYokePos;
  itsYokeLowerPos.setX(itsPolePos.x()+(itsPoleSize.x()-itsYokeUpperSize.x())/2.0);
  itsYokeLowerPos.setY(itsYokePos.y()-itsYokeMiddleSize.y()/2.0-itsYokeUpperSize.y()/2.0);
		 
  itsUpstreamCoilLength=285*CLHEP::mm; //The part of the could that sticks out of the upstream/downstream end;
  itsDownstreamCoilLength=320*CLHEP::mm;
  //The coil position relative to the marker volume.
  itsCoilPos.set(
		 // (itsYokeSize.x() + itsCoilSize.x())/2.0 - (320+262)*CLHEP::mm,
		 itsYokeSize.x()/2.0 -itsAperture1Size.x() - itsAperture2Size.x() + yokeCoilGap + itsCoilSize.x()/2.0,
		 0,
		 (itsDownstreamCoilLength-itsUpstreamCoilLength)/2.0
		 );
  itsCoilPos += itsYokePos;
  
  
  //Aperture position relative to the centre of the yoke.
  itsAperture2Pos.set(
		     (itsYokeSize.x()-itsAperture2Size.x())/2.0, //C-type magnet, aperture is to the side.
		     0, //Aperture is centred vertically around the magnet centre.
		     0
		     );

  itsAperture1Pos.set(
		      itsYokeSize.x()/2.0-itsAperture2Size.x()-itsAperture1Size.x()/2.0, //C-type magnet, aperture is to the side.
		      0, //Aperture is centred vertically around the magnet centre.
		      0
		      );

  //Upper coil position relative to marker volume centre.
  itsUpperCoilPos.set(
		      itsCoilPos.x(),
		      (itsCoilAperture+itsCoilSize.y())/2.0,
		      itsCoilPos.z()
		      );
  itsLowerCoilPos.set(
		      itsCoilPos.x(),
		      -(itsCoilAperture+itsCoilSize.y())/2.0,
		      itsCoilPos.z()
		      );

  itsUpperFrontCoilPos.set(
			   itsUpperCoilPos.x()+itsCoilSize.x()/2.0+itsMiddleCoilSize.x()/2.0,
			   itsUpperCoilPos.y(),
			   itsCoilPos.z()+itsCoilSize.z()/2.0-itsMiddleCoilSize.z()/2.0
			   );

  itsLowerFrontCoilPos.set(
			   itsUpperCoilPos.x()+itsCoilSize.x()/2.0+itsMiddleCoilSize.x()/2.0,
			   itsLowerCoilPos.y(),
			   itsCoilPos.z()+itsCoilSize.z()/2.0-itsMiddleCoilSize.z()/2.0
			   );

  itsUpperRearCoilPos.set(
			  itsUpperCoilPos.x()+itsCoilSize.x()/2.0+itsMiddleCoilSize.x()/2.0,
			  itsUpperCoilPos.y(),
			  itsCoilPos.z()-itsCoilSize.z()/2.0+itsMiddleCoilSize.z()/2.0
			  );
  
  itsLowerRearCoilPos.set(
			  itsUpperCoilPos.x()+itsCoilSize.x()/2.0+itsMiddleCoilSize.x()/2.0,
			  itsLowerCoilPos.y(),
			  itsCoilPos.z()-itsCoilSize.z()/2.0+itsMiddleCoilSize.z()/2.0
			  );
  

  itsUpperLeftCoilPos.set(
			  itsUpperCoilPos.x()+itsAperture2Size.x()+itsCoilSize.x(),
			  itsUpperCoilPos.y(),
			  itsUpperCoilPos.z()
			  );
  
  itsLowerLeftCoilPos.set(
			  itsLowerCoilPos.x()+itsAperture2Size.x()+itsCoilSize.x(),
			  itsLowerCoilPos.y(),
			  itsLowerCoilPos.z()
			  );
  
  itsUpperPolePos.set(
		      itsPolePos.x(),
		      itsPolePos.y()+(itsPoleSize.y()+itsPoleAperture)/2.0,
		      itsPolePos.z()
		      );

  itsLowerPolePos.set(
		      itsPolePos.x(),
		      itsPolePos.y()-(itsPoleSize.y()+itsPoleAperture)/2.0,
		      itsPolePos.z()
		      );
  
}

void BDSAwakeSpectrometer::SetVisAttributes()
{
  itsVisAttributes=new G4VisAttributes(G4Colour(0.3,0.4,0.2));
  itsVisAttributes->SetForceWireframe(true);

  _visAttFront=new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.5));
  _visAttScint=new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.5));
  _visAttBase =new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.5));
  _visAttSampler=new G4VisAttributes(G4Colour(0.2,0.2,0.0,0.5));
  

  _visAttFront->SetForceSolid(true);
  _visAttScint->SetForceSolid(true);
  _visAttBase->SetForceSolid(true);
  _visAttSampler->SetForceSolid(true);
  _visAttSampler->SetVisibility(true);
}

void BDSAwakeSpectrometer::BuildMagnet(){
    BuildYoke();
    BuildCoils();
}

void BDSAwakeSpectrometer::BuildCoils(){
  G4ThreeVector size = itsCoilSize;
  G4VSolid*  coilSolid = new G4Box("coilSolid",size.x()/2.0,size.y()/2.0,size.z()/2.0);
  G4VSolid*  middleCoilSolid = new G4Box("middleCoilSolid",itsMiddleCoilSize.x()/2.0,itsMiddleCoilSize.y()/2.0,itsMiddleCoilSize.z()/2.0);
  /*
  G4ThreeVector subtractPos;
  G4VSolid* coilSubtract = new G4Box("coilSubtract", itsAperture2Size.x()/2.0, size.y()/2.0,size.z()/2.0); 
  subtractPos.set(
		  itsPolePos.x()-itsCoilPos.x(),
		  0,
		  size.z()/2.0
		  );
  G4VSolid* coilSolid = new G4SubtractionSolid("coilSolid", coilSolid1, coilSubtract, nullRotationMatrix, 
					       subtractPos);
  */
  
  itsCoilLog = new G4LogicalVolume(coilSolid, BDSMaterials::Instance()->GetMaterial("G4_Cu"),"itsCoilLog",0,0,0);

  itsMiddleCoilLog = new G4LogicalVolume(middleCoilSolid, BDSMaterials::Instance()->GetMaterial("G4_Cu"),"itsMiddleCoilLog",0,0,0);
  
  G4VisAttributes* CoilVisAtt = new G4VisAttributes(G4Color(0.0,0.5,0.5,0.5));
  CoilVisAtt->SetForceSolid(true);
  CoilVisAtt->SetVisibility(true);
  itsCoilLog->SetVisAttributes(CoilVisAtt);
  itsMiddleCoilLog->SetVisAttributes(CoilVisAtt);
}

void BDSAwakeSpectrometer::PlaceCoils(){
  if(itsCoilLog == NULL){
    BuildCoils();
  }
  new G4PVPlacement(_magRotationMatrix,itsUpperCoilPos,itsCoilLog,"CoilUpper",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  new G4PVPlacement(_magRotationMatrix,itsLowerCoilPos,itsCoilLog,"CoilLower",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  new G4PVPlacement(_magRotationMatrix,itsUpperLeftCoilPos,itsCoilLog,"CoilUpperLeft",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  new G4PVPlacement(_magRotationMatrix,itsLowerLeftCoilPos,itsCoilLog,"CoilLowerLeft",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());

  new G4PVPlacement(_magRotationMatrix,itsUpperFrontCoilPos,itsMiddleCoilLog,"CoilUpperFront",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  new G4PVPlacement(_magRotationMatrix,itsLowerFrontCoilPos,itsMiddleCoilLog,"CoilLowerFront",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());

  new G4PVPlacement(_magRotationMatrix,itsUpperRearCoilPos,itsMiddleCoilLog,"CoilUpperRear",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  new G4PVPlacement(_magRotationMatrix,itsLowerRearCoilPos,itsMiddleCoilLog,"CoilLowerRear",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
}



void BDSAwakeSpectrometer::BuildYoke(){
  G4VSolid* YokeUpperSolid = new G4Box("YokeSolid1",itsYokeUpperSize.x()/2.0,itsYokeUpperSize.y()/2.0,itsYokeUpperSize.z()/2.0);

  itsYokeUpperLog = new G4LogicalVolume(YokeUpperSolid, BDSMaterials::Instance()->GetMaterial("G4_Fe"),"itsYokeUpperLog",0,0,0);

  G4VSolid* YokeMiddleSolid = new G4Box("YokeSolid1",itsYokeMiddleSize.x()/2.0,itsYokeMiddleSize.y()/2.0,itsYokeMiddleSize.z()/2.0);

  itsYokeMiddleLog = new G4LogicalVolume(YokeMiddleSolid, BDSMaterials::Instance()->GetMaterial("G4_Fe"),"itsYokeMiddleLog",0,0,0);

  G4VSolid* YokeLowerSolid = new G4Box("YokeSolid1",itsYokeLowerSize.x()/2.0,itsYokeLowerSize.y()/2.0,itsYokeLowerSize.z()/2.0);

  itsYokeLowerLog = new G4LogicalVolume(YokeLowerSolid, BDSMaterials::Instance()->GetMaterial("G4_Fe"),"itsYokeLowerLog",0,0,0);

  G4VisAttributes* YokeVisAtt = new G4VisAttributes(G4Color(0,1,0,0.5));
  YokeVisAtt->SetForceSolid(true);
  YokeVisAtt->SetVisibility(true);
  itsYokeUpperLog->SetVisAttributes(YokeVisAtt);
  itsYokeMiddleLog->SetVisAttributes(YokeVisAtt);
  itsYokeLowerLog->SetVisAttributes(YokeVisAtt);
}

void BDSAwakeSpectrometer::PlaceYoke(){
  if(itsYokeMiddleLog == NULL){
    BuildYoke();
  }
  new G4PVPlacement(_magRotationMatrix,itsYokeUpperPos,itsYokeUpperLog,"YokeUpper",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  new G4PVPlacement(_magRotationMatrix,itsYokeMiddlePos,itsYokeMiddleLog,"YokeMiddle",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  new G4PVPlacement(_magRotationMatrix,itsYokeLowerPos,itsYokeLowerLog,"YokeLower",
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
}

void BDSAwakeSpectrometer::BuildBPFieldAndStepper(){
  // set up the magnetic field and stepper
  G4cout << __METHOD_NAME__ << " - Bfield = " << _BField << G4endl;
  itsMagField=new BDSUniformMagField(G4ThreeVector(0,-_BField,0));
  G4ThreeVector pmin;
  //Set the extent of the magnetic field to within the magnet pole region.
  itsMagField->SetFieldExtent(itsPolePos-itsAperture2Size/2.0, itsPolePos+itsAperture2Size/2.0);
  itsEqRhs=new G4Mag_UsualEqRhs(itsMagField);  
  itsStepper = new G4NystromRK4(itsEqRhs);
  /*
  itsStepper = new BDSDipoleStepper(itsEqRhs);
  */
  //  BDSDipoleStepper* dipoleStepper = dynamic_cast<BDSDipoleStepper*>(itsStepper);

  /*
  BDSDipoleStepper* dipoleStepper = (BDSDipoleStepper*)itsStepper;
  dipoleStepper->SetBField(-_BField); // note the - sign...
  dipoleStepper->SetBGrad(0);
  */
}

void BDSAwakeSpectrometer::SetBPFieldMgr(){
  BuildFieldMgr(itsStepper, itsMagField);
  itsMarkerLogicalVolume->SetFieldManager(itsFieldMgr,true);
}


void BDSAwakeSpectrometer::BuildVacuumChamber(){
  switch(_vacuumChamberType){
  case 0:
    _vacChamb=NULL;
  case 1:
    _vacChamb = new BDSSpectrVacChamb(itsName + "_vacChamb",
				      itsLength,
				      _poleStartZ,
				      _screenEndZ-std::abs(cos(_screenAngle)*_totalThickness),
				      _screenWidth,
				      _screenAngle,
				      _vacInnerWidth,
				      _vacInnerHeight,
				      _vacThickness);
    break;
  default:
    G4String exceptionString = (G4String)"vacuumChamberType: " + _vacuumChamberType + (G4String)" unknown.";
    G4Exception(exceptionString.c_str(), "-1", FatalErrorInArgument, "");
  }
}


void BDSAwakeSpectrometer::PlaceVacuumChamber(){
  if(_vacuumChamberType!=0){
    if(!_vacChamb){
      BuildVacuumChamber();
    }
    _vacChamb->Place(itsMarkerLogicalVolume);
  }
}

void BDSAwakeSpectrometer::BuildCameraScoringPlane(){
  G4String tmp = "_cameraScoringPlane";
  _scoringPlaneName=itsName+tmp;
  int nThisSampler= BDSSampler::GetNSamplers() + 1;
  G4String ident="_camera";
  _samplerName = ("Sampler_"+BDSGlobalConstants::Instance()->StringFromInt(nThisSampler)+"_"+_scoringPlaneName);
  _samplerName2 = ("Sampler_"+BDSGlobalConstants::Instance()->StringFromInt(nThisSampler)+"_"+_scoringPlaneName+"_2");

  
  //Build and place the volume...
  itsCameraScoringPlaneSolid = new G4Box("CameraScoringPlaneSolid",100*CLHEP::mm/2.0,500*CLHEP::mm/2.0,_scoringPlaneThickness/2.0);

  itsCameraScoringPlaneLog = new G4LogicalVolume(itsCameraScoringPlaneSolid,BDSMaterials::Instance()->GetMaterial("vacuum"),"CameraScoringPlaneLog",0,0,0);
  itsCameraScoringPlaneLog->SetVisAttributes(_visAttSampler);

  G4double dispX=_cameraScreenDist-_scoringPlaneThickness/2.0;
  G4double dispY=0;
  G4double dispZ=-_cameraScreenDist/2.0;;


  new G4PVPlacement(BDSGlobalConstants::Instance()->RotY90(),G4ThreeVector(dispX,dispY,dispZ),itsCameraScoringPlaneLog,_samplerName,
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  itsCameraScoringPlaneLog2 = new G4LogicalVolume(itsCameraScoringPlaneSolid,BDSMaterials::Instance()->GetMaterial("vacuum"),"CameraScoringPlaneLog2",0,0,0);
  itsCameraScoringPlaneLog2->SetVisAttributes(_visAttSampler);

  G4double dispX2=-sin(_screenAngle)*_cameraScreenDist;
  G4double dispY2=0;
  G4double dispZ2=cos(_screenAngle)*_cameraScreenDist-_cameraScreenDist/2.0;


  new G4PVPlacement(_screenRotationMatrix,G4ThreeVector(dispX2,dispY2,dispZ2),itsCameraScoringPlaneLog2,_samplerName2,
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  
  (*LogVol)[_samplerName]=itsCameraScoringPlaneLog;
  (*LogVol)[_samplerName2]=itsCameraScoringPlaneLog2;

  itsCameraScoringPlaneLog->SetSensitiveDetector(BDSSampler::GetSensitiveDetector());
  itsCameraScoringPlaneLog2->SetSensitiveDetector(BDSSampler::GetSensitiveDetector());
  //SPM bdsOutput->nSamplers++;
  BDSSampler::AddExternalSampler(_samplerName+"_1");
  BDSSampler::AddExternalSampler(_samplerName2+"_1");

  _samplerName3 = ("Sampler_"+BDSGlobalConstants::Instance()->StringFromInt(nThisSampler)+"_"+_scoringPlaneName+"_3");
  _samplerName4 = ("Sampler_"+BDSGlobalConstants::Instance()->StringFromInt(nThisSampler)+"_"+_scoringPlaneName+"_4");

  
  //Build and place the volume...
  itsCameraScoringPlaneLog3 = new G4LogicalVolume(itsCameraScoringPlaneSolid,BDSMaterials::Instance()->GetMaterial("vacuum"),"CameraScoringPlaneLog3",0,0,0);
  itsCameraScoringPlaneLog3->SetVisAttributes(_visAttSampler);

  G4double dispX3=_cameraScreenDist/2.0-_scoringPlaneThickness/2.0;
  G4double dispY3=0;
  G4double dispZ3=-_cameraScreenDist/2.0;;

  new G4PVPlacement(BDSGlobalConstants::Instance()->RotY90(),G4ThreeVector(dispX3,dispY3,dispZ3),itsCameraScoringPlaneLog3,_samplerName3,
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  itsCameraScoringPlaneLog4 = new G4LogicalVolume(itsCameraScoringPlaneSolid,BDSMaterials::Instance()->GetMaterial("vacuum"),"CameraScoringPlaneLog4",0,0,0);
  itsCameraScoringPlaneLog4->SetVisAttributes(_visAttSampler);

  G4double dispX4=-sin(_screenAngle)*_cameraScreenDist/2.0;
  G4double dispY4=0;
  G4double dispZ4=cos(_screenAngle)*_cameraScreenDist/2.0-_cameraScreenDist/2.0;


  new G4PVPlacement(_screenRotationMatrix,G4ThreeVector(dispX4,dispY4,dispZ4),itsCameraScoringPlaneLog4,_samplerName4,
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  
  (*LogVol)[_samplerName3]=itsCameraScoringPlaneLog3;
  (*LogVol)[_samplerName4]=itsCameraScoringPlaneLog4;
  itsCameraScoringPlaneLog3->SetSensitiveDetector(BDSSampler::GetSensitiveDetector());
  itsCameraScoringPlaneLog4->SetSensitiveDetector(BDSSampler::GetSensitiveDetector());
  BDSSampler::AddExternalSampler(_samplerName3+"_1");
  BDSSampler::AddExternalSampler(_samplerName4+"_1");

  _samplerName5 = ("Sampler_"+BDSGlobalConstants::Instance()->StringFromInt(nThisSampler)+"_"+_scoringPlaneName+"_5");
  _samplerName6 = ("Sampler_"+BDSGlobalConstants::Instance()->StringFromInt(nThisSampler)+"_"+_scoringPlaneName+"_6");

  
  //Build and place the volume...
  itsCameraScoringPlaneLog5 = new G4LogicalVolume(itsCameraScoringPlaneSolid,BDSMaterials::Instance()->GetMaterial("vacuum"),"CameraScoringPlaneLog5",0,0,0);
  itsCameraScoringPlaneLog5->SetVisAttributes(_visAttSampler);

  G4double dispX5=_cameraScreenDist/4.0-_scoringPlaneThickness/2.0;
  G4double dispY5=0;
  G4double dispZ5=-_cameraScreenDist/2.0;;

  new G4PVPlacement(BDSGlobalConstants::Instance()->RotY90(),G4ThreeVector(dispX5,dispY5,dispZ5),itsCameraScoringPlaneLog5,_samplerName5,
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  itsCameraScoringPlaneLog6 = new G4LogicalVolume(itsCameraScoringPlaneSolid,BDSMaterials::Instance()->GetMaterial("vacuum"),"CameraScoringPlaneLog6",0,0,0);
  itsCameraScoringPlaneLog6->SetVisAttributes(_visAttSampler);

  G4double dispX6=-sin(_screenAngle)*_cameraScreenDist/4.0;
  G4double dispY6=0;
  G4double dispZ6=cos(_screenAngle)*_cameraScreenDist/4.0-_cameraScreenDist/2.0;


  new G4PVPlacement(_screenRotationMatrix,G4ThreeVector(dispX6,dispY6,dispZ6),itsCameraScoringPlaneLog6,_samplerName6,
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  
  (*LogVol)[_samplerName5]=itsCameraScoringPlaneLog5;
  (*LogVol)[_samplerName6]=itsCameraScoringPlaneLog6;
  itsCameraScoringPlaneLog5->SetSensitiveDetector(BDSSampler::GetSensitiveDetector());
  itsCameraScoringPlaneLog6->SetSensitiveDetector(BDSSampler::GetSensitiveDetector());
  BDSSampler::AddExternalSampler(_samplerName5+"_1");
  BDSSampler::AddExternalSampler(_samplerName6+"_1");

#ifndef NOUSERLIMITS
  G4double maxStepFactor=0.5;
  G4UserLimits* itsScoringPlaneUserLimits =  new G4UserLimits();
  itsScoringPlaneUserLimits->SetMaxAllowedStep(_scoringPlaneThickness*maxStepFactor);
  itsCameraScoringPlaneLog->SetUserLimits(itsScoringPlaneUserLimits);
#endif
}

//void BDSAwakeSpectrometer::BuildFresnelLens(){
  ////////////////////////////////////////////////////////////////////////////////////////////////////////                                                    
  /*
  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#           Building the Fresnel lens ...            #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;

  G4double      LensDiameter        = 457*CLHEP::mm ; // Size of the optical active area of the lens.                                                                
  G4int      LensNumOfGrooves    = 13 ;
  //G4int      LensNumOfGrooves    = 129 ;                                                                                                                    
  //G4int      LensNumOfGrooves    = 1287 ;                                                                                                                   

  G4double      LensBorderThickness = 2.8*CLHEP::mm ;     // Thickness of the border area.                                                                           
  G4double      LensFocalLength     = 441.973*CLHEP::mm ; // This parameter depends on the lens geometry, etc !!                                                     
  G4Material   *LensMaterial        = G4Material::GetMaterial(name = "Acrylic") ;
  G4ThreeVector LensPosition        = UVscopePosition+G4ThreeVector(0.0*CLHEP::mm,0.0*CLHEP::mm,UVscopeHeight/2.0-UVscopeBaffle) ;

  UltraFresnelLens *FresnelLens = new UltraFresnelLens(LensDiameter,LensNumOfGrooves,LensMaterial,_log) ;
  */
//}


void BDSAwakeSpectrometer::BuildScreenScoringPlane(){
  G4String tmp = "_screenScoringPlane";
  _screenScoringPlaneName=itsName+tmp;
  int nThisSampler= BDSSampler::GetNSamplers() + 1;
  G4String ident="_screen";
  _screenSamplerName = ("Sampler_"+BDSGlobalConstants::Instance()->StringFromInt(nThisSampler)+"_"+_screenScoringPlaneName);
  _screenSamplerName2 = ("Sampler_"+BDSGlobalConstants::Instance()->StringFromInt(nThisSampler)+"_"+_screenScoringPlaneName+"_2");
  
  //Build and place the volume...
  itsScreenScoringPlaneSolid = new G4Box("ScreenScoringPlaneSolid",_screenWidth/2.0,_screenHeight/2.0,_scoringPlaneThickness/2.0);
  itsScreenScoringPlaneLog = new G4LogicalVolume(itsScreenScoringPlaneSolid,BDSMaterials::Instance()->GetMaterial("vacuum"),"ScreenScoringPlaneLog",0,0,0);
  itsScreenScoringPlaneLog->SetVisAttributes(_visAttSampler);
  itsScreenScoringPlaneLog2 = new G4LogicalVolume(itsScreenScoringPlaneSolid,BDSMaterials::Instance()->GetMaterial("vacuum"),"ScreenScoringPlaneLog2",0,0,0);
  itsScreenScoringPlaneLog2->SetVisAttributes(_visAttSampler);
  G4double dispX=0;
  G4double dispY=0;
  G4double dispZ=2*std::cos(std::abs(_screenAngle))*(_screenThickness/2.0+_scoringPlaneThickness/2.0)-_cameraScreenDist/2.0;
  G4double dispZ2=-2*std::cos(std::abs(_screenAngle))*(_screenThickness/2.0+_scoringPlaneThickness/2.0)-_cameraScreenDist/2.0;
  new G4PVPlacement(_screenRotationMatrix,G4ThreeVector(dispX,dispY,dispZ),itsScreenScoringPlaneLog,_screenSamplerName,
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());

  new G4PVPlacement(_screenRotationMatrix,G4ThreeVector(dispX,dispY,dispZ2),itsScreenScoringPlaneLog2,_screenSamplerName2,
		    itsMarkerLogicalVolume,false,0,BDSGlobalConstants::Instance()->GetCheckOverlaps());
  
  //--
  (*LogVol)[_screenSamplerName]=itsScreenScoringPlaneLog;
  
  (*LogVol)[_screenSamplerName2]=itsScreenScoringPlaneLog2;
  
  //--
  itsScreenScoringPlaneLog->SetSensitiveDetector(BDSSampler::GetSensitiveDetector());
  //-----------
  itsScreenScoringPlaneLog2->SetSensitiveDetector(BDSSampler::GetSensitiveDetector());
  //SPM bdsOutput->nSamplers++;
  //--
  BDSSampler::AddExternalSampler(_screenSamplerName+"_1");
  //----------
  BDSSampler::AddExternalSampler(_screenSamplerName2+"_1");
#ifndef NOUSERLIMITS
  G4double maxStepFactor=0.5;
  G4UserLimits* itsScoringPlaneUserLimits =  new G4UserLimits();
  itsScoringPlaneUserLimits->SetMaxAllowedStep(_scoringPlaneThickness*maxStepFactor);
  itsScreenScoringPlaneLog->SetUserLimits(itsScoringPlaneUserLimits);
#endif
}

void BDSAwakeSpectrometer::Build(){
      SetVisAttributes(); 
      BuildScreen();
      BuildCamera();	
      CalculateLengths();
      BuildMarkerLogicalVolume();
      //      BuildScreenScoringPlane();
      //      BuildCameraScoringPlane();
      PlaceScreen();
      //      PlaceCamera();
      //      if(BDSGlobalConstants::Instance()->GetBuildTunnel()){
      //	BuildTunnel();
      //      }
      AddSensitiveVolume(itsMarkerLogicalVolume);
      switch(_magnetGeometryType){
      case 0:
	break;
      case 1:
	BuildMagnet();
	PlaceMagnet();
	break;
      default:
	G4String exceptionString = (G4String)"magnetGeometryType: " + _magnetGeometryType + (G4String)" unknown.";
	G4Exception(exceptionString.c_str(), "-1", FatalErrorInArgument, "");
      }
      BuildVacuumChamber();
      PlaceVacuumChamber();
      BuildFieldAndStepper();
}

void BDSAwakeSpectrometer::PlaceMagnet(){
    PlaceYoke();
  PlaceCoils();
}

void BDSAwakeSpectrometer::BuildCamera(){
  _camera=new BDSCCDCamera();
}
void BDSAwakeSpectrometer::PlaceCamera(){
  _camera->phys(new G4PVPlacement(_screenRotationMatrix,
				  G4ThreeVector(-1*_cameraScreenDist*sin(_screenAngle),0,1*_cameraScreenDist*cos(_screenAngle)),
				  _camera->log(),
				  _camera->name()+"_phys",
				  itsMarkerLogicalVolume,
				  false,
				  0,
				  true)
		);
}

void BDSAwakeSpectrometer::BuildScreen()
{
  G4cout << "Building BDSAwakeMultilayerScreen...." << G4endl;
  G4double grainSize = 25*1e-6*CLHEP::m;
  _mlScreen = new BDSAwakeMultilayerScreen(_material,_thickness, _windowScreenGap ,grainSize, _windowThickness, _windowMaterial);
  
  G4cout << "finished." << G4endl;
  //  if(BDSGlobalConstants::Instance()->GetSensitiveComponents()){
  //    for(int i=0; i<_mlScreen->nLayers(); i++){
  //      AddSensitiveVolume(_mlScreen[i].log());
  //    }
  //  } 
  G4cout << "BDSAwakeSpectrometer: finished building screen" << G4endl;
}

void BDSAwakeSpectrometer::PlaceScreen(){
  _mlScreen->place(_screenRotationMatrix,
		   G4ThreeVector(_screenCentreX,0,_screenCentreZ),
		   itsMarkerLogicalVolume
		   );
}

void BDSAwakeSpectrometer::CalculateLengths(){
  std::cout << __METHOD_NAME__ << std::endl;
  BDSAcceleratorComponent::CalculateLengths();
  //-------
  //Screen dimensions.
  _screenWidth=_mlScreen->size().x();
  _screenHeight=_mlScreen->size().y();
  std::cout << "... got screen dimensions... " << std::endl;
  //The scoring plane...
  _scoringPlaneThickness=1*CLHEP::um;
  _screenThickness = _mlScreen->size().z();
  _totalThickness = _screenThickness + 2*_scoringPlaneThickness;
  G4double z_wid = _screenWidth * std::sin(std::abs(_screenAngle));//Length due to the screen width and angle
  G4double z_thi = _totalThickness * std::cos(std::abs(_screenAngle));//Length due to the screen thickness
  _screen_z_dim = z_wid+z_thi;
  G4double x_wid = _screenWidth * std::cos(std::abs(_screenAngle));//Length due to the screen width and angle
  G4double x_thi = _totalThickness * std::sin(std::abs(_screenAngle));//Length due to the screen thickness
  _screen_x_dim = x_wid+x_thi;
  //Vacuum chamber dimensions.
  _vacThickness=2*CLHEP::mm;
  _vacHeight=6.7*CLHEP::cm;
  _vacInnerWidth=_vacHeight-2*_vacThickness;
  _vacInnerHeight=_vacHeight-2*_vacThickness;

  _startZPos = -itsLength/2.0;
  //Pole position
  _poleStartZ += _startZPos;
  //Screen position
  _screenEndZ += _poleStartZ;
  _screenCentreZ = _screenEndZ - _screen_z_dim/2.0;
  _screenCentreX = _screen_x_dim/2.0 + _vacInnerWidth/2.0 + _vacThickness;
  
  /*
  itsXLength = itsYLength = BDSGlobalConstants::Instance()->GetComponentBoxSize()/2;
  itsXLength = std::max(itsXLength, this->GetTunnelRadius()+2*std::abs(this->GetTunnelOffsetX()) + BDSGlobalConstants::Instance()->GetTunnelThickness()+BDSGlobalConstants::Instance()->GetTunnelSoilThickness() + 4*BDSGlobalConstants::Instance()->GetLengthSafety() );   
  itsYLength = std::max(itsYLength, this->GetTunnelRadius()+2*std::abs(BDSGlobalConstants::Instance()->GetTunnelOffsetY()) + BDSGlobalConstants::Instance()->GetTunnelThickness()+BDSGlobalConstants::Instance()->GetTunnelSoilThickness()+4*BDSGlobalConstants::Instance()->GetLengthSafety() );
  */

  _cameraScreenDist=(4.0)*CLHEP::m;
  //  _cameraScreenDist=4*213*CLHEP::mm;


  
  //  G4double thi=_totalThickness+2*_cameraScreenDist+2*_camera->size().z()+2*_scoringPlaneThickness;

  MagnetDefaults();  


  //  itsXLength = (_screen_x_dim + 2*_vacWidth1)+2*_cameraScreenDist;
  //  itsXLength = std::max(2*(std::abs(_screenCentreX)+_screen_x_dim), itsYokeSize.x()+2*std::abs(itsPolePos.x()));
  //  itsYLength = std::max(std::max(_screenHeight,_camera->size().y()),itsYokeSize.y());
  //  itsYLength = std::max(itsYLength,50*CLHEP::cm);
  std::cout << __METHOD_NAME__ << " " << itsLength << " " << itsXLength << " " << itsYLength << std::endl;

  _vacDispZ2=(-itsLength)/2.0+(_vacWidth2)/2.0;


  _vacLength=itsLength;


  std::cout << __METHOD_END__ << std::endl;
}

void BDSAwakeSpectrometer::BuildMarkerLogicalVolume(){
  itsMarkerSolidVolume=new G4Box( itsName+"_marker_solid",
				  itsXLength/2.0,
				  itsYLength/2.0,
				  itsLength/2.0); //z half length 

  itsMarkerLogicalVolume=new G4LogicalVolume
    (itsMarkerSolidVolume, 
     BDSMaterials::Instance()->GetMaterial("vacuum"),
     itsName+"_marker_log");
  G4VisAttributes* visAtt = new G4VisAttributes(G4Color(0,1,0));
  visAtt->SetForceWireframe(true);
  visAtt->SetVisibility(true);
  itsMarkerLogicalVolume->SetVisAttributes(visAtt);
#ifndef NOUSERLIMITS
  G4double maxStepFactor=0.5;
  itsMarkerUserLimits =  new G4UserLimits();
  itsMarkerUserLimits->SetMaxAllowedStep(itsLength*maxStepFactor);
  itsMarkerUserLimits->SetUserMinEkine(BDSGlobalConstants::Instance()->GetThresholdCutCharged());
  itsMarkerLogicalVolume->SetUserLimits(itsMarkerUserLimits);
#endif
}

BDSAwakeSpectrometer::~BDSAwakeSpectrometer()
{
  delete _mlScreen;
  delete _camera;
  delete _vacRotationMatrix;
}

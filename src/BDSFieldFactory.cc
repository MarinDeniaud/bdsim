#include "BDSFieldType.hh"
#include "BDSDebug.hh"
#include "BDSExecOptions.hh"
#include "BDSFieldFactory.hh"
#include "BDSFieldMagSolenoid.hh"
#include "BDSGeometry.hh"
#include "BDSGeometrySQL.hh"
#include "BDSGeometryType.hh"
#include "BDSFieldMagDecapole.hh"
#include "BDSFieldMagMultipole.hh"
#include "BDSFieldMagMuonSpoiler.hh"
#include "BDSFieldMagOctupole.hh"
#include "BDSFieldMagQuadrupole.hh"
#include "BDSFieldMagSBend.hh"
#include "BDSFieldMagSextupole.hh"
#include "BDSFieldMagSolenoid.hh"
#include "BDSIntegratorDecapole.hh"
#include "BDSIntegratorDipole.hh"
#include "BDSIntegratorOctupole.hh"
#include "BDSIntegratorQuadrupole.hh"
#include "BDSIntegratorSextupole.hh"
#include "BDSIntegratorSolenoid.hh"
#include "BDSIntegratorTeleporter.hh"
#include "BDSMagFieldMesh.hh"
#include "BDSMagFieldSQL.hh"
#include "BDSMagFieldXY.hh"
#include "BDSMagField3D.hh"
#include "BDSMagnetType.hh"
#include "BDSUtilities.hh"

#include "globals.hh" // geant4 types / globals
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4SimpleRunge.hh"
#include "G4UniformMagField.hh"

#include <typeinfo>
#include <utility>

BDSFieldFactory* BDSFieldFactory::instance = nullptr;

BDSFieldFactory* BDSFieldFactory::Instance()
{
  if (instance == nullptr)
    {instance = new BDSFieldFactory();}
  return instance;
}

BDSFieldFactory::BDSFieldFactory()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  offset            = G4ThreeVector(0,0,0);
  format            = BDSFieldType::zero;
  fileName          = "";
  geometry          = nullptr;
  cacheLength       = 1*CLHEP::um;

  CleanUp();
}

void BDSFieldFactory::CleanUp()
{
  field         = nullptr;
  bField        = nullptr;
  eqOfMotion    = nullptr;
  bEqOfMotion   = nullptr;
  integrator    = nullptr;
  completeField = nullptr;
}

BDSMagFieldMesh* BDSFieldFactory::CreateMagneticField(G4String      formatAndFilePath,
						      G4ThreeVector offsetIn,
						      BDSGeometry*  geometryIn)
  
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  offset            = offsetIn;
  geometry          = geometryIn;

  std::pair<G4String, G4String> ff = BDS::SplitOnColon(formatAndFilePath);
  fileName = BDS::GetFullPath(ff.first);
  format   = BDS::DetermineFieldType(ff.second);
  
  if (format == BDSFieldType::threed)
    {return CreateMagField3D();}
  
  else if (format == BDSFieldType::xy)
    {return CreateMagFieldXY();}
  
  else if (format == BDSFieldType::zero)
    {
    if(geometry)
      {
	if (geometry->Format() == BDSGeometryType::lcdd)
	  {return CreateMagFieldLCDD();}
      }
    }

  // this logic seems unclear
  else if (format == BDSFieldType::mokka || format == BDSFieldType::zero)
    {
      if(geometry)
	{
	  if (geometry->Format() == BDSGeometryType::mokka)
	    {return CreateMagFieldSQL();}
	}
    }
  return nullptr;
}

BDSFieldObjects* BDSFieldFactory::CreateFieldMagEquation(const BDSMagnetType      type,
							 BDSMagnetStrength* const strength,
							 const G4double           brho)
{
  CleanUp();
  
  // switch on type and build correct field
  switch (type.underlying())
    {
    case BDSMagnetType::solenoid:
      CreateSolenoid(strength,   brho); break;
    case BDSMagnetType::sectorbend:
    case BDSMagnetType::rectangularbend: // TBC need to think about length for this
      CreateDipole(strength,     brho); break;
    case BDSMagnetType::quadrupole:
      CreateQuadrupole(strength, brho); break;
    case BDSMagnetType::sextupole:
      CreateSextupole(strength,  brho); break;
    case BDSMagnetType::octupole:
      CreateOctupole(strength,   brho); break;
    case BDSMagnetType::decapole:
      CreateDecapole(strength,   brho); break;
    case BDSMagnetType::multipole:
      CreateMultipole(strength,  brho); break;
    case BDSMagnetType::vkicker:
      CreateKicker(strength, brho, true); break;
    case BDSMagnetType::hkicker:
      CreateKicker(strength, brho, false); break;
    case BDSMagnetType::muonspoiler:
      CreateMuonSpoiler(strength, brho); break;
    default:
      G4cerr << __METHOD_NAME__ << "not an equation based field type - not building field" << G4endl;
      break; // this will return nullptr
    }
  return completeField;
}

BDSFieldObjects* BDSFieldFactory::CreateFieldMagOuter(const BDSMagnetType      type,
						      BDSMagnetStrength* const strength,
						      const G4double           brho)
{
  CleanUp();

  // switch on the type and build correct field
  switch (type.underlying())
    {
    case BDSMagnetType::sectorbend:
    case BDSMagnetType::rectangularbend:
    case BDSMagnetType::quadrupole:
    case BDSMagnetType::sextupole:
    case BDSMagnetType::octupole:
    case BDSMagnetType::decapole:
    case BDSMagnetType::multipole:
      CreateOuterMultipole(type, strength, brho); break;
    case BDSMagnetType::solenoid:
    case BDSMagnetType::vkicker:
    case BDSMagnetType::hkicker:
      break; // return a nullptr
    case BDSMagnetType::muonspoiler:
      CreateMuonSpoiler(strength, brho); break;
    default:
      G4cerr << __METHOD_NAME__ << "no outer field defined for this type of magnet" << G4endl;
      break; // return a nullptr
    }
  return completeField;
}

BDSFieldObjects* BDSFieldFactory::CreateTeleporter(G4ThreeVector teleporterDelta)
{
  bField      = new BDSMagFieldMesh(); //Zero magnetic field.
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new BDSIntegratorTeleporter(bEqOfMotion, teleporterDelta);
  CommonConstructor();
  return completeField;
}

void BDSFieldFactory::CreateOuterMultipole(const BDSMagnetType      type,
					   BDSMagnetStrength* const strength,
					   const G4double           brho)
{;}
  
void BDSFieldFactory::CommonConstructor()
{
  completeField = new BDSFieldObjects(bField, bEqOfMotion, integrator);
}

void BDSFieldFactory::CreateSolenoid(BDSMagnetStrength* strength,
				     G4double           brho)
{
  bField      = new BDSFieldMagSolenoid(strength, brho);
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new BDSIntegratorSolenoid(strength, brho, bEqOfMotion);
  CommonConstructor();
}

void BDSFieldFactory::CreateDipole(BDSMagnetStrength* strength,
				   G4double           brho)
{
  bField      = new BDSFieldMagSBend(strength, brho);
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new BDSIntegratorDipole(strength, brho, bEqOfMotion);
  CommonConstructor();
}

void BDSFieldFactory::CreateQuadrupole(BDSMagnetStrength* strength,
				       G4double           brho)
{
  bField      = new BDSFieldMagQuadrupole(strength, brho);
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new BDSIntegratorQuadrupole(strength, brho, bEqOfMotion);
  CommonConstructor();
}

void BDSFieldFactory::CreateSextupole(BDSMagnetStrength* strength,
				      G4double           brho)
{
  bField      = new BDSFieldMagSextupole(strength, brho);
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new BDSIntegratorSextupole(strength, brho, bEqOfMotion);
  CommonConstructor();
}

void BDSFieldFactory::CreateOctupole(BDSMagnetStrength* strength,
				     G4double           brho)
{
  bField      = new BDSFieldMagOctupole(strength, brho);
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new BDSIntegratorOctupole(strength, brho, bEqOfMotion);
  CommonConstructor();
}

void BDSFieldFactory::CreateDecapole(BDSMagnetStrength* strength,
				     G4double           brho)
{
  bField      = new BDSFieldMagDecapole(strength, brho);
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new BDSIntegratorDecapole(strength, brho, bEqOfMotion);
  CommonConstructor();
}

void BDSFieldFactory::CreateMultipole(BDSMagnetStrength* strength,
				      G4double           brho)
{
  bField      = new BDSFieldMagMultipole(strength, brho);
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new G4SimpleRunge(bEqOfMotion);
  CommonConstructor();
}

void BDSFieldFactory::CreateKicker(BDSMagnetStrength* strength,
				   G4double           brho,
				   G4bool             isVertical)
{
  if (isVertical)
    {bField = new BDSFieldMagSBend(strength, brho, G4ThreeVector(-1,0,0));}
  else
    {bField = new BDSFieldMagSBend(strength, brho);}
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new G4ClassicalRK4(bEqOfMotion);
  //integrator  = new BDSIntegratorDipole(strength, brho, bEqOfMotion);
  CommonConstructor();
}

void BDSFieldFactory::CreateMuonSpoiler(BDSMagnetStrength* strength,
					G4double           brho)
{
  bField      = new BDSFieldMagMuonSpoiler(strength, brho);
  bEqOfMotion = new G4Mag_UsualEqRhs(bField);
  integrator  = new G4ClassicalRK4(bEqOfMotion);
  CommonConstructor();
}

BDSMagFieldMesh* BDSFieldFactory::CreateMagFieldXY()
{return new BDSMagFieldXY(fileName);}

BDSMagFieldMesh* BDSFieldFactory::CreateMagField3D()
{return new BDSMagField3D(fileName, offset);}

BDSMagFieldMesh* BDSFieldFactory::CreateMagFieldSQL()
{
  if(geometry->HasFields() || !fileName.empty())
    {
      // Check for field file or volumes with fields
      // as there may be cases where there are no formats given
      // in gmad file but fields might be set to volumes in SQL files
      return new BDSMagFieldSQL(fileName,
				geometry->Length(),
				geometry->QuadVolBgrad(),
				geometry->SextVolBgrad(),
				geometry->OctVolBgrad(),
				geometry->UniformFieldVolField(),
				geometry->NPoleField(),
				geometry->HasUniformField());
    }
  else
    {return nullptr;}
}

BDSMagFieldMesh* BDSFieldFactory::CreateMagFieldLCDD()
{return (BDSMagFieldMesh*)(geometry->Field());}

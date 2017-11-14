#include "BDSMagnetOuterInfo.hh"

#include "BDSMagnetGeometryType.hh"
#include "globals.hh"             // geant4 globals / types
#include "G4Material.hh"

BDSMagnetOuterInfo::BDSMagnetOuterInfo():
  name("not_specified"),
  geometryType(BDSMagnetGeometryType::cylindrical),
  outerDiameter(0),
  outerMaterial(nullptr),
  vhRatio(1.0),
  angleIn(0),
  angleOut(0),
  yokeOnLeft(false),
  hStyle(false),
  buildEndPieces(true),
  coilWidthFraction(0.65),
  coilHeightFraction(0.8),
  geometryTypeAndPath("")
{;}

BDSMagnetOuterInfo::BDSMagnetOuterInfo(G4String              nameIn,
				       BDSMagnetGeometryType geometryTypeIn,
				       G4double              outerDiameterIn,
				       G4Material*           outerMaterialIn,
				       G4double              vhRatioIn,
				       G4double              angleInIn,
				       G4double              angleOutIn,
				       G4bool                yokeOnLeftIn,
				       G4bool                hStyleIn,
				       G4bool                buildEndPiecesIn,
				       G4double              coilWidthFractionIn,
				       G4double              coilHeightFractionIn,
				       G4String              geometryTypeAndPathIn):
  name(nameIn),
  geometryType(geometryTypeIn),
  outerDiameter(outerDiameterIn),
  outerMaterial(outerMaterialIn),
  vhRatio(vhRatioIn),
  angleIn(angleInIn),
  angleOut(angleOutIn),
  yokeOnLeft(yokeOnLeftIn),
  hStyle(hStyleIn),
  buildEndPieces(buildEndPiecesIn),
  coilWidthFraction(coilWidthFractionIn),
  coilHeightFraction(coilHeightFractionIn),
  geometryTypeAndPath(geometryTypeAndPathIn)
{;}
  
std::ostream& operator<< (std::ostream& out, BDSMagnetOuterInfo const& info)
{
  out << "Magnet Outer Info:  \"" << info.name << "\""             << G4endl;
  out << "Geometry Type:       "  << info.geometryType             << G4endl;
  out << "Diameter:            "  << info.outerDiameter            << G4endl;
  out << "Material:            "  << info.outerMaterial->GetName() << G4endl;
  out << "V / H Ratio:         "  << info.vhRatio                  << G4endl;
  out << "Angle In:            "  << info.angleIn                  << G4endl;
  out << "Angle Out:           "  << info.angleOut                 << G4endl;
  out << "Yoke on Left:        "  << info.yokeOnLeft               << G4endl;
  out << "H Style (C if not):  "  << info.hStyle                   << G4endl;
  out << "Build end pieces:    "  << info.buildEndPieces           << G4endl;
  out << "Coil Width Fraction  "  << info.coilWidthFraction        << G4endl;
  out << "Coil Height Fraction "  << info.coilHeightFraction       << G4endl;
  out << "Geometry:            "  << info.geometryTypeAndPath      << G4endl;
  return out;
}

#include "BDSFieldMagGradient.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "BDSInterpolatorType.hh"
#include "BDSFieldMag.hh"
#include "BDSMagnetStrength.hh"
#include "G4String.hh"
#include "G4Transform3D.hh"
#include "BDSFieldMag.hh"
#include "BDSFieldMagInterpolated2D.hh"
#include "BDSFieldMagSkew.hh"

#include <string>
#include <vector>

BDSFieldMagGradient::BDSFieldMagGradient()
{;}

G4double BDSFieldMagGradient::GetBy(BDSFieldMag* field, G4double x, G4double y) const
{
  G4ThreeVector position(x, y, 0);
  G4ThreeVector fieldAtXY = field->GetField(position);
  G4double by = fieldAtXY[1]/CLHEP::tesla;
  return by;
}

BDSMagnetStrength* BDSFieldMagGradient::CalculateMultipoles(BDSFieldMag* BField,
							    G4int        order,
							    G4double     Brho)
{
    G4cout << "running field gradient calculations" << G4endl;
    BDSMagnetStrength* outputstrengths = new BDSMagnetStrength();
    G4double h =0.5; //distance apart in CLHEP distance units (mm) to place query points.

    G4double brhoinv = 1./Brho;
    G4int centreIndex = 0;
    std::vector<G4double> d = PrepareValues(BField, 5, 0, h, centreIndex);
    G4int centreIndexSkew = 0;
    std::vector<std::vector<G4double>> dskew = PrepareSkewValues(BField,5,0,h,centreIndexSkew);
    // o+1 as we start from k1 upwards - ie, 0th order isn't calculated
    for (G4int o = 0; o < order; ++o)
      {
        (*outputstrengths)["k" + std::to_string(o+1)] = Derivative(d, o+1, centreIndex, h) * brhoinv * pow(CLHEP::meter,o+1);
        G4cout << "k" << o+1 << " = " << (*outputstrengths)["k" + std::to_string(o+1)] << G4endl;
        (*outputstrengths)["k" + std::to_string(o+1) + "s"] = Derivative(dskew[o], o+1, centreIndex, h)*brhoinv*pow(CLHEP::meter,o+1);
        G4cout << "k" << o+1 << "s"<< " = " << (*outputstrengths)["k" + std::to_string(o+1) + "s"]<< G4endl;
      }
    return outputstrengths;
}

std::vector<G4double> BDSFieldMagGradient::PrepareValues(BDSFieldMag* field,
							 G4int        order,
							 G4double     centreX,
							 G4double     h,
							 G4int&       centreIndex) const
{

  G4int maxN = 2*order + 1;
  centreIndex = maxN; // write out maxN to centre index
  std::vector<G4double> data(2*maxN+1); // must initialise vector as not using push_back
  
  for (G4int i = -maxN; i <= maxN; i++)
    {
    data[maxN + i] = GetBy(field, centreX+(G4double)i*h);
    }
  return data;
}

std::vector<std::vector<G4double>> BDSFieldMagGradient::PrepareSkewValues(BDSFieldMag* field,
                                                                          G4int        order,
                                                                          G4double     centreX,
                                                                          G4double     h,
                                                                          G4int&       centreIndex) const
{
    G4double rotation[5] = {CLHEP::pi/4, CLHEP::pi/6, CLHEP::pi/8, CLHEP::pi/10, CLHEP::pi/12};
    G4int maxN = 2*order + 1;
    centreIndex = maxN;
    std::vector<G4double> data(2*maxN+1);
    std::vector<std::vector<G4double>> SkewValues(order+1);
    for (G4int j=0; j<order; j++)
    {
        BDSFieldMagSkew* skewField = new BDSFieldMagSkew(field, rotation[j]);
        for (G4int i = -maxN; i < maxN; i++)
        {
            data[maxN + i] = GetBy(skewField, centreX + (G4double) i * h);
            SkewValues[j].push_back(data[maxN + i]);
        }
        delete skewField;
    }
    return SkewValues;
}



G4double BDSFieldMagGradient::Derivative(const std::vector<G4double>& data,
					 const G4int                  order,
					 const G4int                  startIndex,
					 const G4double               h) const
{
  if (order == 0)
    {return data.at(startIndex);}
  G4int subOrder = order-1;
  G4double result = -Derivative(data, subOrder,startIndex+2,h)
    + 8*Derivative(data, subOrder,startIndex+1, h)
    - 8*Derivative(data, subOrder,startIndex-1, h)
    + Derivative(data, subOrder,startIndex-2, h);
  result /= 12*(h);
  return result;
}
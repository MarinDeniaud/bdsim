/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2018.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "BDSDebug.hh"
#include "BDSFieldMagDipole.hh"
#include "BDSGlobalConstants.hh"
#include "BDSIntegratorDipoleRodrigues2.hh"
#include "BDSIntegratorDipoleQuadrupole.hh"
#include "BDSIntegratorQuadrupole.hh"
#include "BDSMagnetStrength.hh"
#include "BDSStep.hh"
#include "BDSUtilities.hh"

#include "globals.hh"
#include "G4AffineTransform.hh"
#include "G4Mag_EqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>

BDSIntegratorDipoleQuadrupole::BDSIntegratorDipoleQuadrupole(BDSMagnetStrength const* strengthIn,
							     G4double                 brhoIn,
							     G4Mag_EqRhs*             eqOfMIn,
							     G4double minimumRadiusOfCurvatureIn,
							     const G4double&          tiltIn):
  BDSIntegratorMag(eqOfMIn, 6),
  bRho(brhoIn),
  eq(static_cast<BDSMagUsualEqRhs*>(eqOfM)),
  bPrime(std::abs(brhoIn) * (*strengthIn)["k1"]),
  beta0((*strengthIn)["beta0"]),
  rho((*strengthIn)["length"]/(*strengthIn)["angle"]),
  fieldRatio((*strengthIn)["field"]/(*strengthIn)["scaling"] / (bRho/rho)),
  nominalEnergy((*strengthIn)["nominalEnergy"]),
  fieldArcLength((*strengthIn)["length"]),
  fieldAngle((*strengthIn)["angle"]),
  tilt(tiltIn),
  dipole(new BDSIntegratorDipoleRodrigues2(eqOfMIn, minimumRadiusOfCurvatureIn))
{
  primaryMass = BDSGlobalConstants::Instance()->BeamParticleDefinition()->Mass();
  primaryCharge = BDSGlobalConstants::Instance()->BeamParticleDefinition()->Charge();
  zeroStrength = !BDS::IsFinite((*strengthIn)["field"]);
  BDSFieldMagDipole* dipoleField = new BDSFieldMagDipole(strengthIn);
  unitField = (dipoleField->FieldValue()).unit();
  delete dipoleField;
  angleForCL = fieldRatio != 1 ? fieldAngle * fieldRatio : fieldAngle;
}

BDSIntegratorDipoleQuadrupole::~BDSIntegratorDipoleQuadrupole()
{
  delete dipole;
}

void BDSIntegratorDipoleQuadrupole::Stepper(const G4double yIn[6],
					    const G4double dydx[6],
					    const G4double h,
					    G4double       yOut[6],
					    G4double       yErr[6])
{
  // charge and unit normalisation
  const G4double fcof = eqOfM->FCof();
  
  // Protect against very small steps, neutral particles, and zero field: drift through.
  if (h < 1e-12 || !BDS::IsFinite(fcof) || zeroStrength)
    {
      AdvanceDriftMag(yIn,h,yOut,yErr);
      SetDistChord(0);
      return;
    }

  // Revert to the backup stepper if the particle charge or mass is not nominal
  G4double charge = eq->FCof()/(eplus*c_light); //same equation as in BDSMagUsualEqRhs
  G4double mass = std::sqrt(eq->Mass()); // Mass() returns mass squared
  if (!((mass == primaryMass) && (charge == primaryCharge)))
    {
      dipole->Stepper(yIn, dydx, h, yOut, yErr); // more accurate step with error
      SetDistChord(dipole->DistChord());
      return;
    }

  // Revert to the backup stepper if the particle momentum is greater than 5% off nominal
  G4double totalEnergy = BDSGlobalConstants::Instance()->BeamParticleDefinition()->TotalEnergy();
  G4double relEnergyDiff = (eq->TotalEnergy(yIn) - totalEnergy)/totalEnergy;
  if (std::abs(relEnergyDiff) > 0.05)
    {
      dipole->Stepper(yIn, dydx, h, yOut, yErr); // more accurate step with error
      SetDistChord(dipole->DistChord());
      return;
    }

  // try out a simple dipole step first - dumb stepper has no error but is 2x as fast as normal step
  dipole->SingleStep(yIn, h, yOut);

  G4double dipoleDC = dipole->DistChord();
  // We assume that the major effect is the dipole component and the quadrupole
  // component is low. Therefore, we can safely set distchord from the dipole.
  SetDistChord(dipoleDC);

  // return just dipole kick for small step size
  if (h < 1e-4) // 1e-4mm
    {
      dipole->Stepper(yIn, dydx, h, yOut, yErr); // more accurate step with error
      SetDistChord(dipole->DistChord());
      return;
    }

  // If the particle might spiral, we return and just use the dipole only component
  // Aimed at particles of much lower momentum than the design energy.
  G4double radiusOfCurvature = dipole->RadiusOfHelix();
  if (dipoleDC > 0.3*radiusOfCurvature)
    {
      dipole->Stepper(yIn, dydx, h, yOut, yErr); // more accurate step with error
      SetDistChord(dipole->DistChord());
      return;
    }
  
  // not going to spiral so proceed
  // convert to true curvilinear
  G4ThreeVector globalPos   = G4ThreeVector(yIn[0], yIn[1], yIn[2]);
  G4ThreeVector globalMom   = G4ThreeVector(yIn[3], yIn[4], yIn[5]);

  BDSStep       localCL     = GlobalToCurvilinear(fieldArcLength, unitField, angleForCL,
						  globalPos, globalMom, h, true, fcof, tilt);
  G4ThreeVector localCLPos  = localCL.PreStepPoint();
  G4ThreeVector localCLMom  = localCL.PostStepPoint();
  G4ThreeVector localCLMomU = localCLMom.unit();

  // only proceed with thick matrix if particle is paraxial
  // judged by forward momentum > 0.9 and |transverse| < 0.1
  if (localCLMomU.z() < 0.9 || std::abs(localCLMomU.x()) > 0.1 || std::abs(localCLMomU.y()) > 0.1)
    {
      dipole->Stepper(yIn, dydx, h, yOut, yErr); // more accurate step with error
      SetDistChord(dipole->DistChord());
      return;
    }
  
  // calculate new position
  G4ThreeVector localCLPosOut;
  G4ThreeVector localCLMomOut;
  OneStep(localCLPos, localCLMom, localCLMomU, h, fcof, localCLPosOut, localCLMomOut);

  // convert to global coordinates for output
  BDSStep globalOut = CurvilinearToGlobal(fieldArcLength, unitField, angleForCL,
					  localCLPosOut, localCLMomOut, true, fcof, tilt);

  G4ThreeVector globalPosOut = globalOut.PreStepPoint();
  G4ThreeVector globalMomOut = globalOut.PostStepPoint();

  // error along direction of travel really
  G4ThreeVector globalMomOutU = globalMomOut.unit();
  globalMomOutU *= 1e-10;
  
  // chord distance (simple quadratic approx)
  G4double dc = h*h/(8*radiusOfCurvature);
  if (std::isnan(dc))
    {SetDistChord(0);}
  else
    {SetDistChord(dc);}
  
  // write out values and errors
  for (G4int i = 0; i < 3; i++)
    {
      yOut[i]     = globalPosOut[i];
      yOut[i + 3] = globalMomOut[i];
      yErr[i]     = globalMomOutU[i];
      yErr[i + 3] = 1e-40;
    }
}

void BDSIntegratorDipoleQuadrupole::OneStep(const G4ThreeVector& posIn,
					    const G4ThreeVector& momIn,
					    const G4ThreeVector& momUIn,
					    const G4double&      h,
					    const G4double&      fcof,
					    G4ThreeVector&       posOut,
					    G4ThreeVector&       momOut) const
{
  G4double momInMag = momIn.mag();
  G4double nomMomentum = std::abs(bRho * fcof); // safe as brho is nominal and abs(fcof) is scaling factor
  G4double energy = eq->TotalEnergy(momIn);

  // get beta (v/c)
  //G4double beta = eq->Beta(momIn);

  // deltaE/P0 to match literature.
  G4double deltaEoverP = (energy - nominalEnergy) / (nomMomentum);

  // quad strength k normalised to charge and nominal momentum
  // eqOfM->FCof() gives us conversion to MeV,mm and rigidity in Tm correctly
  // as well as charge of the given particle
  G4double K1  = std::abs(fcof)*bPrime/nomMomentum;

  // separate focussing strengths for vertical and horizontal axes.
  // Used by matrix elements so must be derived from nominal values.
  G4double kx2 = std::pow(1.0 / rho, 2) + K1;
  G4double kx  = std::sqrt(std::abs(kx2));
  G4double ky2 = -K1;
  G4double ky  = std::sqrt(std::abs(ky2));
  G4double kxl = kx * h;
  G4double kyl = ky * h;

  G4bool focussing = K1 >= 0; // depends on charge as well (in eqOfM->FCof())

  G4double x0  = posIn.x();
  G4double y0  = posIn.y();
  G4double s0  = posIn.z();
  G4double xp  = momUIn.x();
  G4double yp  = momUIn.y();
  G4double zp  = momUIn.z();

  G4double x1  = x0;
  G4double y1  = y0;
  G4double xp1 = xp;
  G4double yp1 = yp;
  G4double zp1 = zp;

  G4double X11=0,X12=0,X16=0,X21=0,X22=0,X26=0;
  G4double Y11=0,Y12=0,Y21=0,Y22=0;

  // matrix elements. All must derived from nominal parameters.
  if (focussing)
    {//focussing
      X11= std::cos(kxl);
      X12= std::sin(kxl)/kx;
      X21=-std::abs(kx2)*X12;
      X22= X11;
      X16 = (1.0/beta0) * ((1.0/rho) / kx2) * (1 - std::cos(kxl));
      X26 = (1.0/beta0) * (1.0/rho) * X12;

      Y11= std::cosh(kyl);
      Y12= std::sinh(kyl)/ky;
      if (std::isnan(Y12))  //y12 can be nan from div by 0 in previous line
        {Y12 = h;}
      Y21= std::abs(ky2)*Y12;
      Y22= Y11;
    }
  else
    {// defocussing
      X11= std::cosh(kxl);
      X12= std::sinh(kxl)/kx;
      X21= std::abs(kx2)*X12;
      X22= X11;
      X16 = (1.0/beta0) * ((1.0/rho) / kx2) * (1 - std::cosh(kxl));
      X26 = (1.0/beta0) * (1.0/rho) * X12;
      
      Y11= std::cos(kyl);
      Y12= std::sin(kyl)/ky;
      if (std::isnan(Y12))  //y12 can be nan from div by 0 in previous line
        {Y12 = h;}
      Y21= -std::abs(ky2)*Y12;
      Y22= Y11;
    }

  x1  = X11*x0 + X12*xp + X16*deltaEoverP;
  xp1 = X21*x0 + X22*xp + X26*deltaEoverP;
  
  y1  = Y11*y0 + Y12*yp;    
  yp1 = Y21*y0 + Y22*yp;

  G4double s1 = s0 + h;
  
  // relies on normalised momenta otherwise this will be nan.
  zp1 = std::sqrt(1 - xp1*xp1 - yp1*yp1);
  if (std::isnan(zp1))
    {zp1 = zp;} // ensure not nan

  G4ThreeVector momOutUnit = G4ThreeVector(xp1, yp1, zp1);
  momOut = momOutUnit * momInMag;

  posOut = G4ThreeVector(x1, y1, s1);

}


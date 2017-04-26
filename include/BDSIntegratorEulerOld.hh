#ifndef BDSINTEGRATOREULEROLD_H
#define BDSINTEGRATOREULEROLD_H

#include "BDSIntegratorMag.hh"

#include "globals.hh"

class G4Mag_EqRhs;

/**
 * @brief Common functionality for Euler integrators.
 *
 * @author Laurie Nevay
 */

class BDSIntegratorEulerOld: public BDSIntegratorMag
{
public:
  explicit BDSIntegratorEulerOld(G4Mag_EqRhs* eqOfMIn);
  virtual ~BDSIntegratorEulerOld(){;}

  /// Call AdvanceHelix pure virtual function. If the local unit momentum
  /// z component is less than 0.9 or the absolute magnitude of the momentum
  /// is < 40MeV, use the backup stepper. Otherwise compare a single full step
  /// and two half steps for the error calculation. Use the value of two half
  /// steps as the output coordinates.
  virtual void Stepper(const G4double yIn[],
		       const G4double dydx[],
		       const G4double h,
		       G4double       yOut[],
		       G4double       yErr[]);
  
protected:
  /// Expected method in derived class to calculate the output coordinates
  /// for given input coordinates along step length h.
  virtual void AdvanceHelix(const G4double yIn[],
			    G4double       h,
			    G4double       yOut[]) = 0;

  /// Advance chord by quadratic approximation. Can throw a std::out_of_range exception
  /// for inappropriate use of paraxial approximation.
  void AdvanceChord(const G4double       h,
		    G4ThreeVector&       localPos,
		    G4ThreeVector&       localMom,
		    const G4ThreeVector& localA);
};

#endif

/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: BDSOctStepper.hh,v 1.1.1.1 2004/11/18 17:42:41 ilia Exp $
// GEANT4 tag $Name:  $
//
//
// class BDSOctStepper
//
// Class description:
//

// History:
// - Created. G.Blair 1/11/00

#ifndef BDSOCTSTEPPER_HH
#define BDSOCTSTEPPER_HH
#include "globals.hh"
#include "G4MagIntegratorStepper.hh"
//#include "G4MagHelicalStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"
#include "G4Navigator.hh"

class BDSOctStepper : public G4MagIntegratorStepper
{

  public:  // with description

    BDSOctStepper(G4Mag_EqRhs *EqRhs);

    ~BDSOctStepper();

      void Stepper( const G4double y[],
		  const G4double dydx[],
		  const G4double h,
		  G4double yout[],
		  G4double yerr[]  );
      // The stepper for the Runge Kutta integration.
      // The stepsize is fixed, equal to h.
      // Integrates ODE starting values y[0 to 6]
      // Outputs yout[] and its estimated error yerr[].

    G4double DistChord()   const;
      // Estimate maximum distance of curved solution and chord ... 
 
     void SetBTrpPrime(G4double aBTrpPrime);
     G4double GetBTrpPrime();

     void StepperName();

  /*
   void DumbStepper( const G4double yInput[],
		     G4ThreeVector Bfld,
		     G4double hstep,
		     G4double yOut[]);

  */
  public: // without description
  
    G4int IntegratorOrder()const { return 2; }

  protected:
    //  --- Methods used to implement all the derived classes -----

    void AdvanceHelix( const G4double  yIn[],
                       G4ThreeVector Bfld,
		       G4double  h,
		       G4double  yOct[]);    // output 
      // A first order Step along a oct inside the field.

  private:
  
    G4ThreeVector yInitial, yMidPoint, yFinal;
      // Data stored in order to find the chord.

    G4Mag_EqRhs*  fPtrMagEqOfMot;

    G4double itsBTrpPrime;
  //    G4ThreeVector itsInitialPoint, itsFinalPoint, itsMidPoint,itsDistVec;
    G4double itsDist;

};

inline  void BDSOctStepper::SetBTrpPrime(G4double aBTrpPrime)
{itsBTrpPrime=aBTrpPrime;
}
inline G4double BDSOctStepper::GetBTrpPrime()
{return itsBTrpPrime;}

inline void BDSOctStepper::StepperName()
{G4cout<<"BDSOctStepper"<<G4endl;}

#endif /* BDSOCTSTEPPER_HH */

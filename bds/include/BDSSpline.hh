/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/

// BDSSpline
// Generates individual particles from PLACET output file
// G.A.Blair (RHUL) 10.12.01
// adapted from code supplied by D.Schulte (CERN)

#ifndef BDSSpline_h
#define BDSSpline_h 1

#include "BDSSpline_tab_entry.hh"
#include <vector>

using std::vector;

typedef vector<BDSSpline_tab_entry*> vTab;
typedef vector<G4double> vDbl;

class BDSSpline
{
public:
  BDSSpline(G4int nIn);
  ~BDSSpline();

  void initialise(vDbl* xIn,G4int xscalIn, 
		  vDbl* yIn,G4int yscalIn);
  G4double integrate(G4double xIn);

  vTab tab;
  G4int n,xscal,yscal;

private:
  vDbl u;
};

#endif

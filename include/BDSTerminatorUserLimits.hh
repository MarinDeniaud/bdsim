#ifndef BDSTERMINATORUSERLIMITS_H
#define BDSTERMINATORUSERLIMITS_H

#include "G4UserLimits.hh"
#include "globals.hh"  // geant4 basic stl types

/**
 * @brief Dynamic user limits for a volume that changed based
 * on the a parameter turnstaken in BDSGlobalConstants(). 
 *
 * Used to control the number of turns a particle completes around
 * a circular accelerator. An instance of these user limits is
 * typically attached to a volume near the end of the lattice.
 * 
 * @author Laurie Nevay
 */

class BDSTerminatorUserLimits: public G4UserLimits
{
public:
  /// Basic inheritance - just use everything normally from G4UserLimits but 
  /// replace one function in inherited class. Default values are defined
  /// in G4UserLimits so all particles continue.
  BDSTerminatorUserLimits(G4double ustepMax = DBL_MAX,
			  G4double utrakMax = DBL_MAX,
			  G4double utimeMax = DBL_MAX,
			  G4double uekinMin = 0.,     
			  G4double urangMin = 0.);    
  
  BDSTerminatorUserLimits(const G4String& type,
			  G4double ustepMax = DBL_MAX,
			  G4double utrakMax = DBL_MAX,
			  G4double utimeMax = DBL_MAX,
			  G4double uekinMin = 0.,
			  G4double urangMin = 0.);
  
  virtual ~BDSTerminatorUserLimits(){};

  /// Replace this function alone in G4UserLimits
  /// where we only set the limit to E=0 eV on the
  /// last turn.
  virtual G4double GetUserMinEkine(const G4Track&);

protected:
  G4double keeprunningEK; ///< Minimum energy particle must have to keep going
  G4double stoprunningEK; ///< Same, so everything < DBL_MAX so everything stopped
  
private:
  G4bool verbose;
  G4int turnsToTake;
};

#endif

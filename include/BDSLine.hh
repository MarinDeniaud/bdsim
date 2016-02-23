#ifndef BDSLINE_H
#define BDSLINE_H 

#include "BDSAcceleratorComponent.hh"

#include <iterator>
#include <vector>

/**
 * @brief a class that hold multiple accelerator components
 * 
 * Not physically used in Geant4, only as placeholder.
 * Therefore virtual methods are not specified, but could be in the future.
 *
 * Memory of components is not owned
 */

class BDSLine: public BDSAcceleratorComponent
{
private:
  /// Private typedef first so it can be used throughout. Typedef for member
  /// vector.
  typedef std::vector<BDSAcceleratorComponent*> BDSLineVector;

  /// Member vector to store components
  BDSLineVector line;
  
public:
  BDSLine(G4String name);
  virtual ~BDSLine(){;}

  /// Add a component to the line.
  void AddComponent(BDSAcceleratorComponent* component);

  /// @{ Iterator mechanics
  typedef BDSLineVector::iterator       iterator;
  typedef BDSLineVector::const_iterator const_iterator;
  iterator       begin()       {return line.begin();}
  iterator       end()         {return line.end();}
  const_iterator begin() const {return line.begin();}
  const_iterator end()   const {return line.end();}
  G4bool         empty() const {return line.empty();}
  /// @}
  
  /// Override the BDSAccelerator::Initialise() function to loop over the
  /// line and call that function belonging to each member
  virtual void Initialise();

private:
  /// define pure virtual method
  virtual void BuildContainerLogicalVolume(){};

  /// assignment and copy constructor not implemented nor used
  BDSLine& operator=(const BDSLine&);
  BDSLine(BDSLine&);
};

#endif

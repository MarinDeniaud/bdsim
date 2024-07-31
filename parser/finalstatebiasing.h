/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2023.

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
#ifndef FINALSTATEBIASING_H
#define FINALSTATEBIASING_H

#include <string>
#include <list>
#include <vector>

namespace GMAD
{
  class Array;
}

namespace GMAD
{
  /**
   * @brief Final state biasing class for parser
   * 
   * @author Jochem Snuverink and Marin Deniaud
   */
  class FinalStateBiasing {
    
  public:
    /// bias name
    std::string name;
    /// particle to bias name
    std::string particle;
    /// geant4 process(es): single string, but can have multiple processes separated with a space
    std::string process;
    std::vector<std::string> processList;
    /// factor(s) corresponding to process
    std::vector<double> factor;
    /// threshold(s) corresponding to process
    std::vector<double> threshold;
    /// product particle(s) that will be affected by the splitting : single string, but can have multiple processes separated with a space
    std::string product;
    std::vector<std::string> productList;

    
    /// constructor
    FinalStateBiasing();
    /// reset
    void clear();
    /// print some properties
    void print()const;

    /// set methods by property name, numeric values
    void set_value(const std::string& property, double value);
    /// set methods by property name, list of numeric values
    void set_value(const std::string& property, Array* value);
    /// set methods by property name, string values
    void set_value(const std::string& property, std::string value);
  };

}

#endif

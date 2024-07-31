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
#include "finalstatebiasing.h"

#include <iomanip>
#include <iostream>
#include <sstream>

#include "array.h"

using namespace GMAD;

FinalStateBiasing::FinalStateBiasing(){}

void FinalStateBiasing::clear()
{
  name = "";
  particle = "";
  process = "";
  processList.clear();
  factor.clear();
  threshold.clear();
  product = "";
  productList.clear();
}

void FinalStateBiasing::print() const
{
  std::cout << "finalstatebiasing: " << std::endl
	    << "name "     << name     << " " << std::endl
	    << "particle " << particle << " " << std::endl
	    << "process "  << process  << " " << std::endl
        << "factor ";
  for (const auto& i : factor) std::cout << static_cast<double>(i) << " ";
  std::cout << std::endl << "threshold ";
  for (const auto& i : threshold) std::cout << static_cast<double>(i) << " ";
  std::cout << std::endl;
  std::cout << "product " << product << " " << std::endl;
}

void FinalStateBiasing::set_value(const std::string& property, double value )
{
#ifdef BDSDEBUG
  std::cout << "parser> Setting value " << std::setw(25) << std::left << property << value << std::endl;
#endif

  if (property=="fsfact") {factor.push_back(value); return;}
  if (property=="threshold") {threshold.push_back(value); return;}
  
  std::cerr << "Error: parser> unknown finalstatebiasing option \"" << property << "\" with value " << value << std::endl;
  exit(1);
}

void FinalStateBiasing::set_value(const std::string& property, Array* value)
{
  if (property=="fsfact")
    {value->set_vector(factor);}
  else if (property=="threshold")
    {value->set_vector(threshold);}
  else
    {
      std::cerr << "Error: parser> unknown finalstatebiasing option \"" << property << "\" with value ";
      for (const auto& i : value->GetData())
	{std::cout << i << " ";}
      std::cout << std::endl;
      exit(1);
    }
}

void FinalStateBiasing::set_value(const std::string& property, std::string value)
{
#ifdef BDSDEBUG
  std::cout << "parser> Setting value " << std::setw(25) << std::left << property << value << std::endl;
#endif

  if (property=="name")
    {name     = value;}
  else if (property=="particle")
    {particle = value;}
  else if ((property=="proc") || (property=="process"))
  {
    process = value;
    std::stringstream ss(process);
    std::string tok;
    while(ss >> tok)
    {processList.push_back(tok);}
    return;
  }
  else if ((property=="prod") || (property=="product"))
  {
    product = value;
    std::stringstream ss(product);
    std::string tok;
    while(ss >> tok)
    {productList.push_back(tok);}
    return;
  }
  else
    {
      std::cerr << "Error: parser> unknown finalstatebiasing option \"" << property
		<< "\" with value " << value << std::endl;
      exit(1);
    }
}

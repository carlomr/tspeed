/**
 * @file Parameters.cpp
 * @brief  Implementation of the Parameters methods
 *
 * @author Carlo Marcati
 * @date 2013-09-08
 */
/* This program is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version. 
 *  
 *  This program is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details. 
 *  
 *  You should have received a copy of the GNU General Public License 
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include"FESpace.hpp"
namespace Tspeed{
 void Parameters::setp(std::string const & p, const unsigned int lab, const double lambda)
    {
	unsigned int count = 0;
	for(auto ie: M_mesh->elements())
	{

	    if(ie.reg() == lab)
	   {
	       (*M_map[p])(ie.id()) = lambda;
	       ++count;
	   }
	}
	std::cout << "Param :: set "<<p << " = " << lambda << " on " << count << " elements." << std::endl;
    };
 double Parameters::avg_p(std::string const & p,int i, int j)const
 {
     double loc_p = (*M_map.at(p))(i);
     double out_p = (*M_map.at(p))(j);
     if(loc_p+out_p == 0)
	 return 0;
     else
	 return 2*loc_p*out_p/(loc_p+out_p);
 }

}

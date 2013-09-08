/**
 * @file Receivers.cpp
 * @brief Implementation of the Receivers methods
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
#include"Receivers.hpp"

namespace Tspeed{

    void Receivers::add(double const & x, double const &y, unsigned int const & ir, unsigned int const & step)
    {
	//std::cout << M_nel << " " << step+1 << std::endl;
	M_val.resize(M_nel, step+1);
	M_val(ir, step) = std::array<double,2>{{x,y}};

    }
    void Receivers::write(std::string const & fn)const
    {
	std::ofstream outf;
	for(unsigned int i = 0; i<M_nel; ++i )
	{
	    outf.open(fn+"_rcv_"+std::to_string(i)+".out");
	    for(unsigned int j = 0; j<M_val.cols(); ++j )
	    {
		outf << M_val(i,j)[0] << " " << M_val(i,j)[1] << std::endl;
	    }
	    outf.close();
	}
    }

}

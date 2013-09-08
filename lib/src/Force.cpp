/**
 * @file Force.cpp
 * @brief Implementation of the Force method
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
#include"Force.hpp"
namespace Tspeed{
    Eigen::VectorXd PointWiseForce::eval(const double & t)const
    {
		
	std::array<double, 2> F = M_f(t);
	Vec rhs=Vec::Zero(M_totSize);
	rhs.segment(M_startIndex1,M_size) = F[0] * M_shape[0];
	rhs.segment(M_startIndex2,M_size) = F[1] * M_shape[0];
	return rhs;
    }

}

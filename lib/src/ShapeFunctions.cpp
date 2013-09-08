/**
 * @file ShapeFunctions.cpp
 * @brief Implementation of the jacobi polynomials
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
#include"ShapeFunctions.hpp"
typedef Eigen::ArrayXd Arr;


Eigen::ArrayXd Tspeed::jacobi_polynomial(int N, int alpha, int beta, Eigen::ArrayXd const & z)
{
    if(N==0)
	return Arr::Ones(z.size());
    if(N==1)
	return 0.5*(alpha - beta + (alpha + beta + 2.0)*z);

    unsigned int dims = z.size();
    //Arr jf = Arr::Zero(dims);
    const float one = 1.0;
    const float two = 2.0;

    int apb = alpha + beta;

    Arr poly   = Arr::Zero(dims);
    Arr polyn2 = Arr::Ones(dims);
    Arr polyn1 = 0.5*(alpha - beta + (alpha + beta + two)*z);

    double a1, a2, a3, a4;
    for(int k=2; k<=N;++k)
    {
	a1 =  two*k*(k + apb)*(two*k + apb - two);
	a2 = (two*k + apb - one)*(alpha*alpha - beta*beta);
	a3 = (two*k + apb - two)*(two*k + apb - one)*(two*k + apb);
	a4 =  two*(k + alpha - one)*(k + beta - one)*(two*k + apb);

	a2 = a2/a1;
	a3 = a3/a1;
	a4 = a4/a1;

	poly   = (a2 + a3*z)*polyn1 - a4*polyn2;
	polyn2 = polyn1;
	polyn1 = poly;
    }
    return poly;
}


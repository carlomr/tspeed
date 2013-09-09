/**
 * @file TimeAdvance_imp.hpp
 * @brief Implementation of the time advancing template methods
 *
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
#include"TimeAdvance.hpp"

namespace Tspeed
{
    template <int N, typename Q, typename S>
	TimeAdvance::TimeAdvance(FESpace_ptr<N,Q,S> Xh, Parameters const & p, Receivers const & r):M_penalty(2),M_dt(0),M_tmax(0),uh(), uhold(), uholdold(),M_recv(r),M_mat(Xh,p),M_completed(false),M_last_time(0),M_last_step(0),M_recv_written(0),M_nln(Xh->nln()),M_ne(Xh->ne()){
	}
    template <int N, typename Q, typename S>
	LeapFrog::LeapFrog(FESpace_ptr<N,Q,S> Xh, Parameters const & p, Receivers const & r):TimeAdvance(Xh, p, r)
	{
	    B = std::move(M_penalty*N*N*M_mat.getS() + M_mat.getI() + M_mat.getA());
	    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	    //std::cout << M_mat.getI().component(0,0).block(0,0) << std::endl;
	    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	    //std::cout << M_mat.getA().component(0,0).block(0,0) << std::endl;
	    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	    //std::cout << B.component(0,0).block(0,0) << std::endl;
	    //M_mat.getA().sum(M_mat.getI(),B);
	}
    template<int N, typename Q, typename S>
	void TimeAdvance::set_initial_v(FESpace_ptr<N,Q,S> Xh, std::function<std::array<double,2>(double,double)> fun)
	{
	    initial_v = Xh->transform(fun);
	}
    template<int N, typename Q, typename S>
	void TimeAdvance::set_initial_u(FESpace_ptr<N,Q,S> Xh, std::function<std::array<double,2>(double,double)> fun)
	{
	    uholdold = Xh->transform(fun);
	}
}

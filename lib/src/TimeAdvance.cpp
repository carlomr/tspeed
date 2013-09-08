/**
 * @file TimeAdvance.cpp
 * @brief Implmentation of the TimeAdvance methods
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
    void LeapFrog::first_step()
    {
	if(M_last_step != 0)
	    std::cerr << "Warning: attempt to do first step after others" << std::endl; 
	fold = this->M_f->eval(0.);
	//std::cout << M_dt*M_dt/2 *(B*uholdold - f) << std::endl;
	uhold  = uholdold + M_dt*initial_v - M_mat.getinvM()*(M_dt*M_dt/2 *(B*uholdold - fold));
	uh = uhold;
	//std::cout << "+++++++++++++++++++"<<std::endl;
	////std::cout << uholdold << std::endl;
	//std::cout << (B*uholdold ) << std::endl;
	//std::cout << "+++++++++++++++++++"<<std::endl;
	//std::cout << uhold<< std::endl;
    }

    void LeapFrog::step(double t)
    {
	if(t != M_last_step + M_dt)
	{
	    std::cerr << "Warning: provided time is not sequential to last one. Resetting t = " << M_last_step + M_dt << std::endl;
	    t = M_last_step + M_dt;
	}
	std::cout << "Solving timestep at t = " << t << std::endl;

	f = this->M_f->eval(t);;
	uh = 2*uhold - uholdold + M_mat.getinvM() * (- M_dt*M_dt*(B*uhold - f));
	
	update_variables(t);

	if(t >= M_tmax-M_dt )
	    M_completed = true;
	
    }
    Eigen::Matrix2d CTensorProduct(Eigen::Matrix2d const & A, double lambda, double mu)
    {
	Eigen::Matrix2d B;// = new Eigen::Matrix2d;
       	B << (lambda + 2*mu) * A(0,0) + lambda * A(1,1) , 2*mu*A(0,1) , 2*mu*A(0,1) , (lambda + 2*mu) * A(1,1) + lambda * A(0,0);
       	return B;
    }
    double mat_dot(Eigen::Matrix2d const & a, Eigen::Matrix2d const & b){return a(0,0)*b(0,0)+a(1,0)*b(1,0)+a(0,1)*b(0,1)+a(1,1)*b(1,1);};

    void TimeAdvance::eval_receivers()
    {
	double x,y;
	for(int i = 0; i<M_recv.size(); ++i)
	{
	    x = M_recv.shape(i).matrix().dot(uh.segment(M_recv.elem(i)*M_nln, M_nln));
	    y = M_recv.shape(i).matrix().dot(uh.segment(M_recv.elem(i)*M_nln+M_nln*M_ne, M_nln));
	    //std::cout << x << ", " << y << std::endl;
	    M_recv.add(x,y,i,M_recv_written);
	}
	++M_recv_written;
    }
}

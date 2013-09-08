/**
 * @file Force.hpp
 * @brief Header file for the force 
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
#ifndef __FORCE_HPP__
#define __FORCE_HPP__ 1
#include<functional>
#include<Eigen/SparseCore>
#include"Receivers.hpp" //TODO: move pointwiseentity somewhere else, and include it here and in receivers.hpp
#include"FESpace.hpp"
#include<array>

namespace Tspeed{
    //template <typename T>
    /**
     * @brief virtual base class for forces
     */
    class Force
    {
	public:
	    typedef Eigen::SparseVector<double> SPVec;
	    typedef Eigen::VectorXd Vec;
	    Force(){};
	    /**
	     * @brief Constructor, taking the function (time dependent)
	     *
	     * @param fun The function
	     */
	    Force(std::function<std::array<double,2>(const double &)> const & fun);
	    virtual ~Force(){};
	    virtual Vec eval(const double &)const=0;
	protected:
	    std::function<std::array<double,2>(const double &)> M_f;
	    //T M_f;
    };
    typedef std::shared_ptr<Force> Force_ptr;
    //template <typename T>
    /**
     * @brief Time dependent force acting on a point
     */
    class PointWiseForce : public Force, public PointWiseEntity
    {
	public:
	    /**
	     * @brief Costructor taking the function, the point where the force is applied and the function space
	     *
	     * @tparam N,Q,S the template parameters of the function space
	     * @param f the force
	     * @param p the point
	     * @param Xh the space
	     */
	    template<int N,typename Q,typename S>
		PointWiseForce(std::function<std::array<double,2>(const double &)> const& f, Geo::Point p, FESpace_ptr<N,Q,S> Xh );
		//PointWiseForce(T const& , Geo::Point , FESpace<N,Q,S> const & );
	    virtual ~PointWiseForce(){};
	    /**
	     * @brief Get value of force vector, i.e. 
	     * \f[ r_i = \int_K f \psi_i \f]
	     * at time  t, where \f$f\f$ is non null in \f$K\f$
	     *
	     * @param t the time
	     *
	     * @return a vector with elements \f$ r_i \f$
	     */
	    Vec eval(const double & t)const;
	private:
	    unsigned int M_startIndex1;
	    unsigned int M_startIndex2;
	    unsigned int M_size;
	    unsigned int M_totSize;
    };
    //template<typename T>
/*    class DiffuseForce : public Force{*/
	//public: 
	    //template<int N,typename Q,typename S>
		////DiffuseForce(T const& , FESpace<N,Q,S> const & );
		//DiffuseForce(std::function<std::array<double,2>(const double &)> const& , FESpace<N,Q,S> const & );
	    //virtual ~PointWiseForce(){};
	    //Vec eval(const double &)const;

    /*};*/
}

#include"Force_imp.hpp"
#endif

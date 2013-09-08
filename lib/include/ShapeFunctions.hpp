/**
 * @file ShapeFunctions.hpp
 * @brief Header file for the definition of the shape functions
 *
 * A base class is used, and the Dubiner and Boundary adapted derived classes are implemented.
 * See 
 *
 * [1] M. Dubiner, Spectral methods on triangles and other domains, Journal
 * of Scientific Computing 6 (1991), no. 4, 345–390
 * 
 * [2] G. E. Karniadakis and S. J. Sherwin, Spectral/hp element methods
 * for computational fluid dynamics, second ed., Numerical Mathemat
 * ics and Scientific Computation, Oxford University Press, New York.
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
#ifndef __SHAPEFUNCTIONS_HPP__
#define __SHAPEFUNCTIONS_HPP__ 1

#include<functional>
#include<vector>
#include<Eigen/Dense>
namespace Tspeed
{
    /**
     * @brief Base class for the shared functions
     *
     * @tparam N degree of the space \f$\mathbb{P}_N\f$
     */
    template<int N>
	class ShapeFunction
    {
	typedef Eigen::Array<double,Eigen::Dynamic, 2> ArrG;
	typedef Eigen::ArrayXd Arr;
	public:
	/**
	 * @brief Number of degrees of freedom
	 *
	 */
	    enum {gdl = (N+1)*(N+2)/2};
	    /**
	     * @brief Orthonormality of the basis
	     */
	    enum {is_orthonormal = false};
	    virtual ~ShapeFunction(){};
	    ShapeFunction(){};
	    /**
	     * @brief Get value of base function with index s, on points (v,w) 
	     *
	     * @param s index of the basis function
	     * @param v x-coordinates of the points
	     * @param w y-coordinates of the points
	     *
	     * @return An Eigen array with the values
	     */
	    Eigen::ArrayXd phi(unsigned int s, Arr const & v, Arr const &w)const{return M_phi[s](v, w);}; 
	    /**
	     * @brief Get value of base function with index s, on point (x,y) 
	     *
	     * @param s index of the basis function
	     * @param x x-coordinates of the point
	     * @param y y-coordinates of the point
	     *
	     * @return the value
	     */
	    double phi(unsigned int s, double x, double y)const{return *(M_phi[s](x*Arr::Ones(1), y*Arr::Ones(1)).data());};
	    /**
	     * @brief Get values of gradient of basis function s, on points (v,w)
	     *
	     * @param s Index f the function
	     * @param v x-coordinates of the points
	     * @param w y-coordinates of the points
	     *
	     * @return An array of dimension length(v) , 2 with the values
	     */
	    ArrG grad(unsigned int s, Arr const & v, Arr const &w){return M_grad[s](v,w);}; 

	protected:
	    //std::vector<std::function<Eigen::ArrayXd(Eigen::ArrayXd const &,Eigen::ArrayXd const &)>> M_phi;
	    std::vector<std::function<Arr(Arr const &, Arr const &)>> M_phi;
	    //std::vector<std::function<Eigen::Array<double,Dynamic,2>(Eigen::ArrayXd const &, Eigen::ArrayXd const &)>> M_grad;
	    std::vector<std::function<ArrG(Arr const &, Arr const &)>> M_grad;
    };


    /**
     * @brief Dubiner [1] basis
     *
     * @tparam N degree of the space \f$\mathbb{N}\f$
     */
    template<int N>
	class Dubiner:public ShapeFunction<N>
    {
	typedef Eigen::Array<double,Eigen::Dynamic, 2> ArrG;
	typedef Eigen::ArrayXd Arr;
	public:
	/**
	 * @brief Dubiner basis is orthonormal
	 */
	    enum {is_orthonormal = true};
	    virtual ~Dubiner(){};
	    Dubiner();
    };
    /**
     * @brief Boundary adapted [2] basis  
     *
     * @tparam N degree of the space \f$\mathbb{N}\f$
     */
    template<int N>
	class BoundaryAdapted:public ShapeFunction<N>
    {
	typedef Eigen::Array<double,Eigen::Dynamic, 2> ArrG;
	typedef Eigen::ArrayXd Arr;
	public:
	/**
	 * @brief The basis is not orthonormal
	 */
	    enum {is_orthonormal = false};
	    virtual ~BoundaryAdapted(){};
	    BoundaryAdapted();
    };
;
}

#include"ShapeFunctions_imp.hpp"
#endif

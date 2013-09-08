/**
 * @file Receivers.hpp
 * @brief Header file containg the class for receivers and a base class for pointwise entities
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
#ifndef __RECEIVERS_HPP__
#define __RECEIVERS_HPP__ 1
#include<string>
#include"Geometry.hpp"
#include"FESpace.hpp"
#include<fstream>
#include<vector>

namespace Tspeed{
    //template<int N, typename Q, typename S>
    /**
     * @brief A base class for pointwise entities, with the points and the basis function in that points
     */
    class PointWiseEntity{
	public:
	    virtual ~PointWiseEntity(){};
	    /**
	     * @brief All shape functions at point i 
	     *
	     * @param i the index of the point
	     *
	     * @return an array of all functions
	     */
	    Eigen::ArrayXd const & shape(int i)const{return M_shape[i];};
	    /**
	     * @brief Point i, with coordinates in the reference triangle
	     *
	     * @param i the index of the  point
	     *
	     * @return The point
	     */
	    Geo::Point const & point(int i)const{return M_relp[i];};
	    /**
	     * @brief The index of the element where point i resides
	     *
	     * @param i the index of the point
	     *
	     * @return the index of the triangle
	     */
	    unsigned int const & elem(int i)const{return M_ie[i];};
	    /**
	     * @brief The number of points
	     *
	     */
	    unsigned int size()const{return M_ie.size();};
	protected:
	    std::vector<unsigned int> M_ie;
	    std::vector<Geo::Point> M_relp;
	    std::vector<Eigen::ArrayXd> M_shape;
	    unsigned int M_nel;
	    template<int N, typename Q, typename S>
		void M_add(FESpace_ptr<N,Q,S>, Geo::Point const & );

    };
    /**
     * @brief A class for seismic receivers, i.e., receivers recording the movement at a point
     */
    class Receivers: public PointWiseEntity
    {
	public:
	    /**
	     * @brief Constructor taking the function space and a file with the coordinates of the receivers listed (x-coord and y-coord on every row)
	     *
	     * @tparam N,Q,S the template parameters of the function space
	     * @param Xh the space
	     * @param fname the name of the file with the receivers
	     */
	    template<int N, typename Q, typename S>
		Receivers(FESpace_ptr<N,Q,S> Xh,std::string const & fname);
	    /**
	     * @brief Constructor taking the function space and a point
	     *
	     * @tparam N,Q,S the template parameters of the function space
	     * @param Xh the space
	     * @param p the point where the receiver is
	     */
	    template<int N, typename Q, typename S>
		Receivers(FESpace_ptr<N,Q,S> Xh, Geo::Point const & p);
	    /**
	     * @brief Add the the value (x,y) of receiver ir at time step step
	     *
	     * @param x the x displacement recorded
	     * @param y the y displacement recorded
	     * @param ir the index of the receiver
	     * @param step the time step of the simulation
	     */
	    void add(double const & x, double const &y, unsigned int const & ir, unsigned int const & step);
	    /**
	     * @brief Write all recorded values to file
	     *
	     * @param fn the name of the output file
	     */
	    void write(std::string const & fn)const;
	private:
	    //std::vector<std::vector<std::array<double,2>>> M_val;
	    Eigen::Matrix<std::array<double,2>,Eigen::Dynamic, Eigen::Dynamic> M_val;
	    //std::vector<unsigned int> M_ie;
	    //std::vector<Geo::Point> M_relp;
	    //std::vector<Eigen::ArrayXd> M_shape;
	    //template<int N, typename Q, typename S>
		//void M_add(FESpace<N,Q,S> const  &, Geo::Point const & );
    };

}
#include"Receivers_imp.hpp"
#endif

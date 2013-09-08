/**
 * @file QuadratureRule.hpp
 * @brief Header file for the quadrature rules
 *
 * A base class is implemented, with derived classes which implement Gauss quadrature
 * on the triangle and Dunavant quadrature
 *
 * Reference:
 * [1] D. A. Dunavant, High degree efficient symmetrical Gaussian quadra-
 * ture rules for the triangle, Internat. J. Numer. Methods Engrg. 21
 * (1985), no. 6, 1129–1148. 
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
#ifndef __QUADRATURERULE_HPP__
#define __QUADRATURERULE_HPP__ 1
#include<Eigen/Dense> 
#include<limits>
#include<iostream>
#include"Geometry.hpp"
#include"Dunavant.hpp"
//namespace
namespace Tspeed
{
    template<int N>
	constexpr int dunavant_num_points();

    template<> constexpr int dunavant_num_points<1>() {   return 1;}
    template<> constexpr int dunavant_num_points<2>() {   return 3;}
    template<> constexpr int dunavant_num_points<3>() {   return 4;}
    template<> constexpr int dunavant_num_points<4>() {   return 6;}
    template<> constexpr int dunavant_num_points<5>() {   return 7;}
    template<> constexpr int dunavant_num_points<6>() {   return 12;}
    template<> constexpr int dunavant_num_points<7>() {   return 13;}
    template<> constexpr int dunavant_num_points<8>() {   return 16;}
    template<> constexpr int dunavant_num_points<9>() {   return 19;}
    template<> constexpr int dunavant_num_points<10>(){  return 25;}
    template<> constexpr int dunavant_num_points<11>(){  return 27;}
    template<> constexpr int dunavant_num_points<12>(){  return 33;}
    template<> constexpr int dunavant_num_points<13>(){  return 37;}
    template<> constexpr int dunavant_num_points<14>(){  return 42;}
    template<> constexpr int dunavant_num_points<15>(){  return 48;}
    template<> constexpr int dunavant_num_points<16>(){  return 52;}
    template<> constexpr int dunavant_num_points<17>(){  return 61;}
    template<> constexpr int dunavant_num_points<18>(){  return 70;}
    template<> constexpr int dunavant_num_points<19>(){  return 73;}
    template<> constexpr int dunavant_num_points<20>(){  return 79;}
//}
//namespace Tspeed
//{
    /**
     * @brief Base class for quadrature rules
     *
     * @tparam N order of the rule
     */
    template<int N>
	class QuadratureRule
	{
	    protected:
		typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec2;
		typedef Eigen::Matrix<double, Eigen::Dynamic, 2> Mat;
		Vec2 M_w_2D;
		Mat M_node_2D;
		typedef Eigen::Matrix<double, N, 1> Vec;
		Vec M_node_1D;
		Vec M_w_1D;
		size_t   M_nqn_1D;
		size_t   M_nqn_2D;


	    public:
		/**
		 * @brief Required by Eigen
		 */
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
		QuadratureRule() = default;
		/**
		 * @brief Get weights on the edge
		 *
		 * @return The vector of the weights
		 */
		Vec edge_weights()const{return M_w_1D;};
		/**
		 * @brief Get internal weights
		 *
		 * @return  The vector of the weights
		 */
		Vec2 int_weights()const{return M_w_2D;};
		/**
		 * @brief Get i-th qudrature node
		 *
		 * @param i The index of the node
		 *
		 * @return A Point
		 */
		Geo::Point inode(unsigned int i)const{return Geo::Point(M_node_2D(i,0),M_node_2D(i,1));};
		/**
		 * @brief Get i-th weight on the edge
		 *
		 * @param i the index of the weight
		 *
		 * @return the weight
		 */
		double eweight(unsigned int i)const{return M_w_1D[i];};
		/**
		 * @brief Get i-th interior weight 
		 *
		 * @param i the index of the weight
		 *
		 * @return the weight
		 */
		double iweight(unsigned int i)const{return M_w_2D[i];};
		/**
		 * @brief Number of nodes or weights on the edge
		 *
		 */
		size_t edge_n()const{return M_nqn_1D;};
		/**
		 * @brief Number of interior nodes or weights
		 *
		 */
		size_t int_n()const{return M_nqn_2D;};
		/**
		 * @brief Get all interior nodes
		 *
		 * @return A matrix of size number of nodes x 2
		 */
		Mat int_nodes()const{return M_node_2D;};
		/**
		 * @brief Get all nodes on edge i of the reference triangle
		 *
		 * Note that ede nodes always have a Gauss-Legendre distribution,
		 * whatever rule is used in the interior 
		 *
		 * @param i The edge of the reference triangle
		 *
		 * @return A matrix of size number of nodes x 2
		 */
		inline Eigen::Matrix<double, N, 2> edge_nodes(int i)const;
		virtual ~QuadratureRule(){};
		/**
		 * @brief Number of internal nodes/weights
		 *
		 */
		unsigned int nqn2d()const{return M_w_2D.size();};
	};


/**
 * @brief Gauss quadrature rule on the triangle
 *
 * @tparam N the order of the rule
 *
 * The internal nodes will be N^2; N+1 is the required order to integrate 
 * polynomials of order 2N.
 */
    template<int N>
	class Gauss : public QuadratureRule<N>
    {
    public:
	typedef typename QuadratureRule<N>::Vec Vec;
	typedef typename QuadratureRule<N>::Vec Mat;
	typedef typename QuadratureRule<N>::Vec Vec2;
	enum{ nqn2d = N*N};
	enum{ nqn1d = N};
	Gauss();
    };

/**
 * @brief Dunavant [1] quadrature rule.
 *
 *
 * @tparam N the order of the rule.
 *
 * A rule of order N integrates exactly polynomials of order N
 */
    template<int N>
	class Dunavant: public QuadratureRule<N>
    {
    public:
	typedef typename QuadratureRule<N>::Vec Vec;
	typedef typename QuadratureRule<N>::Vec Mat;
	typedef typename QuadratureRule<N>::Vec Vec2;
	enum{ nqn2d = dunavant_num_points<N>()};
	enum{ nqn1d = N};
	Dunavant();
    };

}
#include"QuadratureRule_imp.hpp"
#endif

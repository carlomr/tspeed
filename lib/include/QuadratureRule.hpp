#ifndef __QUADRATURERULE_HPP__
#define __QUADRATURERULE_HPP__ 1
#include<Eigen/Dense> 
#include<limits>
#include<iostream>
#include"Geometry.hpp"
#include"Dunavant.hpp"
namespace Tspeed
{
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
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
		QuadratureRule() = default;
		Vec edge_weights()const{return M_w_1D;};
		Vec2 int_weights()const{return M_w_2D;};
		Geo::Point inode(unsigned int i)const{return Geo::Point(M_node_2D(i,0),M_node_2D(i,1));};
		double eweight(unsigned int i)const{return M_w_1D[i];};
		double iweight(unsigned int i)const{return M_w_2D[i];};
		size_t edge_n()const{return M_nqn_1D;};
		size_t int_n()const{return M_nqn_2D;};
		Mat int_nodes()const{return M_node_2D;};
		inline Eigen::Matrix<double, N, 2> edge_nodes(int i)const;
		virtual ~QuadratureRule(){};
		unsigned int nqn2d()const{return M_w_2D.size();};
	};


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
;


}
#include"QuadratureRule_imp.hpp"
#endif

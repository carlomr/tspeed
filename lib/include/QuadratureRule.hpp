#ifndef __QUADRATURERULE_HPP__
#define __QUADRATURERULE_HPP__ 1
#include<Eigen/Dense> 
#include<limits>
#include<iostream>
#include"Geometry.hpp"
namespace Tspeed
{
    template<int N>
	class QuadratureRule
	{
	    typedef Eigen::Matrix<double, N, 1> Vec;
	    typedef Eigen::Matrix<double, N*N, 1> Vec2;
	    typedef Eigen::Matrix<double, N*N, 2> Mat;
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
	protected:
	    Vec2 M_w_2D;
	    Vec M_w_1D;
	    size_t   M_nqn_1D;
	    size_t   M_nqn_2D;
	    Mat M_node_2D;
	    Vec M_node_1D;

	};


    template<int N>
	class Gauss : public QuadratureRule<N>
    {
	typedef Eigen::Matrix<double, N, 1> Vec;
	typedef Eigen::Matrix<double, N*N, 1> Vec2;
	typedef Eigen::Matrix<double, N*N, 2> Mat;
    public:
	enum{ nqn2d = N*N};
	enum{ nqn1d = N};
	Gauss();
    };


}
#include"QuadratureRule_imp.hpp"
#endif

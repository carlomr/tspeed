#ifndef __SHAPEFUNCTIONS_HPP__
#define __SHAPEFUNCTIONS_HPP__ 1

#include<functional>
#include<vector>
#include<Eigen/Dense>
namespace Tspeed
{
    template<int N>
	class ShapeFunction
    {
	typedef Eigen::Array<double,Eigen::Dynamic, 2> ArrG;
	typedef Eigen::ArrayXd Arr;
	public:
	    enum {gdl = (N+1)*(N+2)/2};
	    enum {is_orthonormal = false};
	    virtual ~ShapeFunction(){};
	    ShapeFunction(){};
	    Eigen::ArrayXd phi(unsigned int s, Arr const & v, Arr const &w)const{return M_phi[s](v, w);}; 
	    double phi(unsigned int s, double x, double y)const{return *(M_phi[s](x*Arr::Ones(1), y*Arr::Ones(1)).data());};
	    ArrG grad(unsigned int s, Arr const & v, Arr const &w){return M_grad[s](v,w);}; 

	protected:
	    //std::vector<std::function<Eigen::ArrayXd(Eigen::ArrayXd const &,Eigen::ArrayXd const &)>> M_phi;
	    std::vector<std::function<Arr(Arr const &, Arr const &)>> M_phi;
	    //std::vector<std::function<Eigen::Array<double,Dynamic,2>(Eigen::ArrayXd const &, Eigen::ArrayXd const &)>> M_grad;
	    std::vector<std::function<ArrG(Arr const &, Arr const &)>> M_grad;
    };


    template<int N>
	class Dubiner:public ShapeFunction<N>
    {
	typedef Eigen::Array<double,Eigen::Dynamic, 2> ArrG;
	typedef Eigen::ArrayXd Arr;
	public:
	    enum {is_orthonormal = true};
	    virtual ~Dubiner(){};
	    Dubiner();
    };
}

#include"ShapeFunctions_imp.hpp"
#endif

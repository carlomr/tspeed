#ifndef __FESPACE_HPP__
#define __FESPACE_HPP__ 1
#include"QuadratureRule.hpp"
#include"ShapeFunctions.hpp"
#include"Mesh.hpp"
#include<Eigen/Dense>
#include<Eigen/StdVector>
#include<functional>
namespace Tspeed
{
    template<int N, typename Q = Gauss<N+1>, typename S = Dubiner<N>>
	class FESpace
	{
	    typedef Eigen::Matrix<double, Q::nqn2d, (N+1)*(N+2)/2> B_mat;
	    typedef Eigen::Matrix<double, 2, Q::nqn2d> G_mat;
	    typedef Eigen::Matrix<double,Q::nqn1d, (N+1)*(N+2)/2>  BEdg_mat;
	    typedef Eigen::Matrix<double,2 ,Q::nqn1d>  GEdg_mat;
	    //typedef std::shared_ptr<Mesh> Mesh_ptr;

	    public:
		explicit FESpace(Mesh_ptr );
		virtual ~FESpace(){};
		Mesh_ptr  mesh()const{return M_mesh;};
		Q const & quad()const{return M_quad;};
		ShapeFunction<N> const & shape()const{return M_sh;};
		unsigned int nln()const{return M_nln;};
		unsigned int ne()const{return M_ne;};
		Eigen::Vector2d  grad(unsigned int k, unsigned int i)const{return M_G[i].col(k);};
		Eigen::VectorXd  b_edge(unsigned int k, unsigned int iedg)const{return M_Bedge[iedg].row(k);};
		Eigen::Vector2d  g_edge(unsigned int k, unsigned int i, unsigned short int edg)const{return M_Gedge[i*3+edg].col(k);};
		Eigen::VectorXd inverse_transform(std::function<std::array<double,2>(double,double)> const &)const;
		double L2error(std::function<std::array<double,2>(double,double)> const &, Eigen::VectorXd const &)const;
		Eigen::VectorXd loc_rhs(Geo::Triangle const & ie, std::function <std::array<double,2>(double,double)> const & fun)const{return M_loc_rhs(ie, fun);};
		void points_out(std::string const &fname)const;
		void field_out(std::string const &fname, Eigen::VectorXd const &uh, unsigned int step)const;
	    private:
		Mesh_ptr M_mesh;
		Q M_quad;
		S M_sh;
		unsigned int M_nln;
		unsigned int M_nqn2d;
		unsigned int M_ne;
		B_mat M_B;
		Eigen::Matrix<double, (N+1)*(N+2)/2,4> M_vert_basis;
		std::array<G_mat, (N+1)*(N+2)/2> M_G; 
		std::array<BEdg_mat, 3> M_Bedge;
		std::array<GEdg_mat, 3*(N+1)*(N+2)/2> M_Gedge;
		Eigen::VectorXd M_loc_rhs(Geo::Triangle const &, std::function <std::array<double,2>(double,double)> const &)const;
	};
    template <int N, typename Q=Gauss<N+1>, typename S=Dubiner<N>>
	using FESpace_ptr = std::shared_ptr<FESpace<N,Q,S>>;
    //template <int N, typename Q=Gauss<N+1>, typename S=Dubiner<N>>
	//using FESpace_ptr_const = std::shared_ptr<const FESpace<N,Q,S>>;

    class Parameters
    {
	typedef Eigen::ArrayXd Arr;
	public:
	    virtual ~Parameters(){};
	    //Parameters(std::shared_ptr<Mesh> const & m):M_mesh(m),M_lambda(new Arr(m->ne())),M_rho(new Arr(m->ne())),M_mu(new Arr(m->ne())),M_map{{{"lambda",M_lambda},{"mu",M_mu},{"rho",M_rho}}}{};
	    Parameters(Mesh_ptr m):M_mesh(m),M_lambda(new Arr(m->ne())),M_rho(new Arr(m->ne())),M_mu(new Arr(m->ne())),M_map{{{"lambda",M_lambda},{"mu",M_mu},{"rho",M_rho}}}{};
	    //void set_lambda(int const , double const);
	    //void set_rho(int const , double const);
	    void setp(std::string const &, int const  , double const );
	    double const & lambda(int i)const{return (*M_lambda)[ i ];}
	    double const & mu(int i)const{return (*M_mu)[ i ];}
	    double const & rho(int i)const{return (*M_rho)[ i ];}
	    double avg_p(std::string const &,int i, int j)const;

	private:
	    //const std::shared_ptr<Mesh> & M_mesh;
	    Mesh_ptr M_mesh;
	    std::shared_ptr<Arr> M_lambda;
	    std::shared_ptr<Arr> M_rho;
	    std::shared_ptr<Arr> M_mu;
	    std::map<std::string, std::shared_ptr<Arr>> M_map;
    };
}
#include"FESpace_imp.hpp"

#endif

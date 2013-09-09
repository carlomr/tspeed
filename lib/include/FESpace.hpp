/**
 * @file FESpace.hpp
 * @brief Header file for the Galerkin space and for the parameters of the elastodynamics equation
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
    /**
     * @brief Functional space
     *
     * @tparam N order of the polynomials
     * @tparam Q quadrature rule
     * @tparam S basis functions
     */
    template<int N, typename Q = Gauss<N+1>, typename S = Dubiner<N>>
	class FESpace
	{
	    typedef Eigen::Matrix<double, Q::nqn2d, (N+1)*(N+2)/2> B_mat;
	    typedef Eigen::Matrix<double, 2, Q::nqn2d> G_mat;
	    typedef Eigen::Matrix<double,Q::nqn1d, (N+1)*(N+2)/2>  BEdg_mat;
	    typedef Eigen::Matrix<double,2 ,Q::nqn1d>  GEdg_mat;
	    //typedef std::shared_ptr<Mesh> Mesh_ptr;

	    public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
		/**
		 * @brief Constructor from the mesh
		 *
		 * @param m shared pointer to Mesh
		 */
		explicit FESpace(Mesh_ptr m);
		virtual ~FESpace(){};
		/**
		 * @brief Get pointer ot mesh
		 *
		 * @return pointer to the mesh
		 */
		Mesh_ptr  mesh()const{return M_mesh;};
		/**
		 * @brief Get quadrature rule
		 *
		 */
		Q const & quad()const{return M_quad;};
		/**
		 * @brief Get Shapefunction
		 *
		 */
		ShapeFunction<N> const & shape()const{return M_sh;};
		/**
		 * @brief Number of degrees of freedom per element
		 *
		 */
		unsigned int nln()const{return M_nln;};
		/**
		 * @brief Number of elements in the mesh
		 *
		 */
		unsigned int ne()const{return M_ne;};
		/**
		 * @brief Get value of the gradient of the basis function i on quadrature point k.
		 *
		 * @param k index of the quadrature node
		 * @param i index of the basis function
		 *
		 * @return A vector with the x and y derivatives
		 */
		Eigen::Vector2d  grad(unsigned int k, unsigned int i)const{return M_G[i].col(k);};
		/**
		 * @brief Get value of all basis function on quadrature node k, on edge iedg
		 *
		 * @param k index of the quadrature point
		 * @param iedg index of the edge (=0,1,2)
		 *
		 * @return  vector of all the values
		 */
		Eigen::VectorXd  b_edge(unsigned int k, unsigned int iedg)const{return M_Bedge[iedg].row(k);};
		/**
		 * @brief value of the gradient of basis function i, on edge edg, quadrature node k
		 *
		 * @param k index of the quadrature node
		 * @param i indec of the basis function
		 * @param edg edg number
		 *
		 * @return a vector with the x and y derivatives
		 */
		Eigen::Vector2d  g_edge(unsigned int k, unsigned int i, unsigned short int edg)const{return M_Gedge[i*3+edg].col(k);};
		/**
		 * @brief Transform a function \f$u(x,y)\f$ into its expantion modes \f$\hat{u}_i\f$ s.t.
		 *  \f[ \sum_i \hat{u}_i\psi_i(x,y) = u(x,y) \f]
		 *
		 * @param fun the function
		 *
		 * @return the vector of \f$\hat{u}_i\f$
		 */
		Eigen::VectorXd transform(std::function<std::array<double,2>(double,double)> const & fun)const;
		/**
		 * @brief L2 norm of the difference uex-uh
		 *
		 * @param uex A function
		 * @param uh A vector of expansion modes
		 *
		 * @return The norm
		 */
		double L2error(std::function<std::array<double,2>(double,double)> const & uex, Eigen::VectorXd const & uh)const;
		/**
		 * @brief Integration of a function against all basis function, in triangle ie, i.e.
		 * \f[ l_i = \int_\mathrm{ie} f \psi_i  \f]
		 *
		 * @param ie the index of the triangle
		 * @param fun the function to be integrated
		 *
		 * @return the vector made by \f$ l_i \f$
		 */
		Eigen::VectorXd loc_rhs(Geo::Triangle const & ie, std::function <std::array<double,2>(double,double)> const & fun)const{return M_loc_rhs(ie, fun);};
		/**
		 * @brief Write a set of points of the mesh to file
		 *
		 * @param fname name of the outpu file
		 */
		void points_out(std::string const &fname)const;
		/**
		 * @brief Write a field to file, on the points given by points_out
		 *
		 * @param fname name of the output file
		 * @param uh vector of the coefficients of a FE function
		 * @param step the time step (gets appended to the file name)
		 */
		void field_out(std::string const &fname, Eigen::VectorXd const &uh, unsigned int step)const;
		/**
		 * @brief Evaluate the mass matrix \f$\mathbf{M}\f$ s.t.
		 * \f[ M_{ij} = \int_{\mathcal{T}^2} \psi_i \psi_j\f]
		 *
		 * @return the matrix \f$ \mathbf{M} \f$
		 */
		Eigen::MatrixXd base_mass()const{return M_base_mass;}
		/**
		 * @brief Evaluate the inverse of the mass matrix \f$ \mathbf{M}^{-1} \f$
		 *
		 * @return \f$ \mathbf{M}^{-1} \f$
		 */
		Eigen::MatrixXd base_invmass()const{return M_base_invmass;}
	    private:
		Mesh_ptr M_mesh;
		Q M_quad;
		S M_sh;
		unsigned int M_nln;
		unsigned int M_nqn2d;
		unsigned int M_ne;
		Eigen::MatrixXd M_base_mass;
		Eigen::MatrixXd M_base_invmass;
		void base_mass_and_inv();
		B_mat M_B;
		Eigen::Matrix<double, (N+1)*(N+2)/2,4> M_vert_basis;
		std::array<G_mat, (N+1)*(N+2)/2> M_G; 
		std::array<BEdg_mat, 3> M_Bedge;
		std::array<GEdg_mat, 3*(N+1)*(N+2)/2> M_Gedge;
		Eigen::VectorXd M_loc_rhs(Geo::Triangle const &, std::function <std::array<double,2>(double,double)> const &)const;
	};
    /**
     * @brief template pointer to FESpace
     *
     * @tparam N,Q,S the same template parameters as FESpace
     */
    template <int N, typename Q=Gauss<N+1>, typename S=Dubiner<N>>
	using FESpace_ptr = std::shared_ptr<FESpace<N,Q,S>>;
    //template <int N, typename Q=Gauss<N+1>, typename S=Dubiner<N>>
	//using FESpace_ptr_const = std::shared_ptr<const FESpace<N,Q,S>>;

    /**
     * @brief Class for the parameters \f$\lambda, \rho, \mu\f$ of the elastodynamics equations
     */
    class Parameters
    {
	typedef Eigen::ArrayXd Arr;
	public:
	    virtual ~Parameters(){};
	    //Parameters(std::shared_ptr<Mesh> const & m):M_mesh(m),M_lambda(new Arr(m->ne())),M_rho(new Arr(m->ne())),M_mu(new Arr(m->ne())),M_map{{{"lambda",M_lambda},{"mu",M_mu},{"rho",M_rho}}}{};
	    Parameters(Mesh_ptr m):M_mesh(m),M_lambda(new Arr(m->ne())),M_rho(new Arr(m->ne())),M_mu(new Arr(m->ne())),M_map{{{"lambda",M_lambda},{"mu",M_mu},{"rho",M_rho}}}{};
	    //void set_lambda(int const , double const);
	    //void set_rho(int const , double const);
	    /**
	     * @brief Set a parameter
	     *
	     * @param p string with the name of the parameter ("lambda", "mu", "rho")
	     * @param lab attribute of the mesh partition on which the parameter is set
	     * @param lambda value of the parameter
	     */
	    void setp(std::string const & p, unsigned int const lab , double const lambda);
	    /**
	     * @brief Get lambda on element i
	     *
	     * @param i the index of the element
	     */
	    double const & lambda(int i)const{return (*M_lambda)[ i ];}
	    /**
	     * @brief Get mu on element i
	     *
	     * @param i the index of the element
	     */
	    double const & mu(int i)const{return (*M_mu)[ i ];}
	    /**
	     * @brief Get rho on element i
	     *
	     * @param i the index of the element
	     */
	    double const & rho(int i)const{return (*M_rho)[ i ];}
	    /**
	     * @brief Get the value of the harmonic average of a parameter between two elements 
	     *
	     * @param p string with the name of the parameter ("lambda", "mu", "rho")
	     * @param i index of one element
	     * @param j index of the seocnd element
	     *
	     * @return harmonic average
	     */
	    double avg_p(std::string const & p ,int i, int j)const;

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

/**
 * @file TimeAdvance.hpp
 * @brief Header file for the implementation of the time stepping and of the matrices for the method
 * @author Carlo Marcati
 * @version 
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
#ifndef __TIMEADVANCE_HPP__
#define __TIMEADVANCE_HPP__ 1

#include<Eigen/SparseCore>
#include<Eigen/Dense>
#include"FESpace.hpp"
#include"Receivers.hpp"
#include"Geometry.hpp"
#include"Force.hpp"
#include"MyMat.hpp"
#include<memory>
#include<limits>

namespace Tspeed
{
    /**
     * @brief Dot product between two 2x2 matrices
     *
     * @param a,b the matrices \f$\mathbf{A}, \mathbf{B}\f$
     *
     * @return \f$\sum_{ij} A_{ij} B_{ij} \f$
     */
    double mat_dot(Eigen::Matrix2d const & a, Eigen::Matrix2d const & b);
    /**
     * @brief tensor product between Hooke's tensor and matrix A
     *
     * @param A the matrix
     * @param lambda, mu Hooke's constants
     *
     * @return 
     */
    Eigen::Matrix2d CTensorProduct(Eigen::Matrix2d const & A, double lambda, double mu);
    /**
     * @brief The class containing the matrices resulting from the spatial discretization
     */
    class Matrices{
	public:
	    typedef Eigen::SparseMatrix<double> SpMat;
	    //Matrices() = default;
	    /**
	     * @brief Constructor taking the space and the parameters
	     *
	     */
	    template <int N, typename Q, typename S>
		Matrices(FESpace_ptr<N,Q,S>, Parameters const &);
	    /**
	     * @brief Get the stiffness matrix 
	     * \f[ A_{ij} = \sum_K \int_K \sigma(\phi_j):\epsilon(\phi_i)  \f]
	     *
	     * @return \f$\mathbf{A}\f$
	     */
	    MyMatMultiDim<MyMatBlockDiag> const & getA(){return A;};
	    /**
	     * @brief Get the stability matrix 
	     * \f[ S_{ij} = \sum_e \int_e [[\phi_j]] : [[\phi_i]]  \f]
	     *
	     * @return \f$\mathbf{S}\f$
	     */
	    MyMatMultiDim<MyMat> const & getS(){return S;};
	    /**
	     * @brief Get the interelement matrix 
	     * \f[ S_{ij} = \sum_e \int_e \lbrace\lbrace\sigma(\phi_j)\rbrace\rbrace : [[\phi_i]]  + \lbrace\lbrace\sigma(\phi_i)\rbrace\rbrace : [[\phi_j]] \f]
	     *
	     * @return \f$\mathbf{S}\f$
	     */
	    MyMatMultiDim<MyMat> const & getI(){return I;};
	    /**
	     * @brief Get inverse of the global mass matrix
	     *
	     * @return \f$ \mathbf{M}^{-1}\f$
	     */
	    MyMatMultiDimBlockDiag<MyMatBlockDiag> const & getinvM(){return invM;};


	private:
	    unsigned int  kien(unsigned int gdl, unsigned int k, unsigned int i, unsigned int e, unsigned int n){return k*gdl*6+i*6+e*2+n;};
	    unsigned int  kin(unsigned int gdl, unsigned int k, unsigned int i, unsigned int n){return k*gdl*2+i*2+n;};
	    MyMatMultiDim<MyMatBlockDiag> A;
	    MyMatMultiDim<MyMat> S;
	    MyMatMultiDim<MyMat> I;
	    MyMatMultiDimBlockDiag<MyMatBlockDiag> M;
	    MyMatMultiDimBlockDiag<MyMatBlockDiag> invM;
	    //Eigen::MatrixXd baseMass;
	    //Eigen::MatrixXd baseInvMass;
	    //void base_mass_and_inv(FESpace_ptr<N,Q,T>);
	    template <int N, typename Q, typename T >
		void sigmaeps_edge(Geo::Triangle const & ie, FESpace_ptr<N,Q,T> Xh, std::array<Eigen::Matrix2d,6*Q::nqn1d*T::gdl> &, std::array<Eigen::Matrix2d,6*Q::nqn1d*T::gdl> &, double lambda, double mu);
	    template <int N, typename Q, typename T >
		void sigmaeps(Geo::Triangle const & ie, FESpace_ptr<N,Q,T> Xh, std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> &, std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> &, double lambda, double mu);
	    template <int N, typename Q, typename T >
		void stress_loc(Geo::Triangle const & , std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> &,std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> &, FESpace_ptr<N,Q,T> );
	    template<int N, typename Q, typename T> 
		void stability(FESpace_ptr<N,Q,T> Xh, Geo::Triangle const & ie,Parameters const & p);
	    template<int N, typename Q, typename T> 
		void interelement(FESpace_ptr<N,Q,T> Xh, Geo::Triangle const & ie,std::array<Eigen::Matrix2d,3*2*Q::nqn1d*T::gdl> &s);
	    template<int N, typename Q, typename T>
		void mass(FESpace_ptr<N,Q,T> Xh, Geo::Triangle const & ie, double rho);
	    
    };
    /**
     * @brief Base class for time stepping methods
     */
    class TimeAdvance
    {

	typedef Eigen::VectorXd Vec;
	public:
	/**
	 * @brief constructor from the space, parameters and receivers
	 *
	 * @tparam N,Q,S template parameters of the space
	 * @param Xh the function space
	 * @param p the parameters of the materials
	 * @param r the receivers
	 */
	    template <int N, typename Q, typename S>
		TimeAdvance(FESpace_ptr<N,Q,S> Xh, Parameters const & p, Receivers const & r);
	/**
	 * @brief The first step of the method (which is different for 2nd order methods)
	 */
	    void first_step();
	    /**
	     * @brief step at time t
	     *
	     * @param t the time 
	     */
	    void step(double t);
	    virtual ~TimeAdvance(){};
	    /**
	     * @brief Set time step \f$\delta t \f$
	     *
	     * @param dt the time step
	     */
	    void set_dt(double dt){M_dt = dt;};
	    /**
	     * @brief Set end time of the simulation
	     *
	     * @param tmax the end time
	     */
	    void set_tmax(double tmax){M_tmax = tmax;};
	    /**
	     * @brief set penalty for the stability matrix
	     *
	     * @param p the penalty value 
	     */
	    void set_penalty(double p){M_penalty = p;};
	    /**
	     * @brief Add the force
	     *
	     * @param f force
	     */
	    void add_force(std::shared_ptr<Force> f){M_f = f;};
	    
	    /**
	     * @brief Set initial speed \f$ \dot{u} \f$
	     *
	     * @param Xh the function space
	     * @param fun \f$ \dot{u}(x,y)\f$
	     */
	    template<int N, typename Q, typename S>
		void set_initial_v(FESpace_ptr<N,Q,S> Xh, std::function<std::array<double,2>(double,double)> fun);
	    /**
	     * @brief Set initial displacement \f$ {u} \f$
	     *
	     * @param Xh the function space
	     * @param fun \f$ {u}(x,y)\f$
	     */
	    template<int N, typename Q, typename S>
		void set_initial_u(FESpace_ptr<N,Q,S> Xh, std::function<std::array<double,2>(double,double)> fun);
	    //template<typename T>
	    //void setf(T f, Geo::Point ){M_f = new(T);};
	    /**
	     * @brief Check if the method has arrived at the final time
	     *
	     * @return TRUE if it is still running
	     */
	    bool is_running(){return !M_completed;};
	    /**
	     * @brief Get the coefficients of the numerical solution \f$u_h\f$
	     *
	     * @return \f$ u_h \f$
	     */
	    Vec const & get_uh()const{return uh;}
	    Vec const & u()const{return uh;}
	    /**
	     * @brief Compute and store the value of the displacement at the receivers
	     */
	    void eval_receivers();
	    /**
	     * @brief Write the time series of the displacement at the receivers
	     *
	     * @param fn Base output file name
	     */
	    void write_receivers(std::string const & fn)const{M_recv.write(fn);};
	protected:
	    double M_penalty;
	    double M_dt;
	    double M_tmax;
	    Vec f;
	    Vec fold;
	    Vec foldold;
	    Vec uh;
	    Vec uhold;
	    Vec uholdold;
	    Vec initial_v;
	    Receivers M_recv;
	    //Force M_f;
	    Matrices M_mat;
	    MyMatMultiDim<MyMat> B;
	    std::shared_ptr<Force> M_f;
	    void update_variables(double t){M_last_step = t;uholdold = uhold; uhold = uh; foldold=fold; fold=f;};
	    bool M_completed;
	    double M_last_step;
	    unsigned int M_recv_written;
	    unsigned int M_nln;
	    unsigned int M_ne;
    };

    /**
     * @brief Implementation of the second order Leap-Frog explicit time stepping scheme
     */
    class LeapFrog : public TimeAdvance
    {
	public:
	template <int N, typename Q, typename S>
	    LeapFrog(FESpace_ptr<N,Q,S>, Parameters const & ,Receivers const &);
	/**
	 * @brief First step for the Leap-Frog method
	 */
	void first_step();
	/**
	 * @brief Step at time t for the Leap-Frog method
	 *
	 * @param t time
	 */
	void step(double);
    };


}
#include"Matrices_imp.hpp"
#include"TimeAdvance_imp.hpp"

#endif

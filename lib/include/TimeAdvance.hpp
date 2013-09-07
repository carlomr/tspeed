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
    double mat_dot(Eigen::Matrix2d const & a, Eigen::Matrix2d const & b);
    Eigen::Matrix2d CTensorProduct(Eigen::Matrix2d const & A, double lambda, double mu);
    class Matrices{
	public:
	    typedef Eigen::SparseMatrix<double> SpMat;
	    //Matrices() = default;
	    template <int N, typename Q, typename S>
		Matrices(FESpace_ptr<N,Q,S>, Parameters const &);
	    MyMatMultiDim<MyMatBlockDiag> const & getA(){return A;};
	    MyMatMultiDim<MyMat> const & getS(){return S;};
	    MyMatMultiDim<MyMat> const & getI(){return I;};
	    MyMatMultiDimBlockDiag<MyMatBlockDiag> const & getinvM(){return invM;};


	private:
	    unsigned int  kien(unsigned int gdl, unsigned int k, unsigned int i, unsigned int e, unsigned int n){return k*gdl*6+i*6+e*2+n;};
	    unsigned int  kin(unsigned int gdl, unsigned int k, unsigned int i, unsigned int n){return k*gdl*2+i*2+n;};
	    MyMatMultiDim<MyMatBlockDiag> A;
	    MyMatMultiDim<MyMat> S;
	    MyMatMultiDim<MyMat> I;
	    MyMatMultiDimBlockDiag<MyMatBlockDiag> M;
	    MyMatMultiDimBlockDiag<MyMatBlockDiag> invM;
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
    class TimeAdvance
    {

	typedef Eigen::VectorXd Vec;
	public:
	    void first_step();
	    void step(double);
	    virtual ~TimeAdvance(){};
	    template <int N, typename Q, typename S>
		TimeAdvance(FESpace_ptr<N,Q,S> , Parameters const &, Receivers const &);
	    void set_dt(double dt){M_dt = dt;};
	    void set_tmax(double tmax){M_tmax = tmax;};
	    void set_penalty(double p){M_penalty = p;};
	    void add_force(std::shared_ptr<Force> f){M_f = f;};
	    
	    template<int N, typename Q, typename S>
		void set_initial_v(FESpace_ptr<N,Q,S> Xh, std::function<std::array<double,2>(double,double)> fun);
	    template<int N, typename Q, typename S>
		void set_initial_u(FESpace_ptr<N,Q,S> Xh, std::function<std::array<double,2>(double,double)> fun);
	    //template<typename T>
	    //void setf(T f, Geo::Point ){M_f = new(T);};
	    bool is_running(){return !M_completed;};
	    Vec const & get_uh()const{return uh;}
	    void eval_receivers();
	    void write_receivers(std::string const & fn)const{M_recv.write(fn);};
	    Vec const & u()const{return uh;}
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

    class LeapFrog : public TimeAdvance
    {
	public:
	template <int N, typename Q, typename S>
	    LeapFrog(FESpace_ptr<N,Q,S>, Parameters const & ,Receivers const &);
	void first_step();
	void step(double);
    };


}
#include"Matrices_imp.hpp"
#include"TimeAdvance_imp.hpp"

#endif

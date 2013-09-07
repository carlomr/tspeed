#ifndef __FORCE_HPP__
#define __FORCE_HPP__ 1
#include<functional>
#include<Eigen/SparseCore>
#include"Receivers.hpp" //TODO: move pointwiseentity somewhere else, and include it here and in receivers.hpp
#include"FESpace.hpp"
#include<array>

namespace Tspeed{
    //template <typename T>
    class Force
    {
	public:
	    typedef Eigen::SparseVector<double> SPVec;
	    typedef Eigen::VectorXd Vec;
	    Force(){};
	    Force(std::function<std::array<double,2>(const double &)> const &);
	    virtual ~Force(){};
	    virtual Vec eval(const double &)const=0;
	protected:
	    std::function<std::array<double,2>(const double &)> M_f;
	    //T M_f;
    };
    //template <typename T>
    class PointWiseForce : public Force, public PointWiseEntity
    {
	public:
	    template<int N,typename Q,typename S>
		PointWiseForce(std::function<std::array<double,2>(const double &)> const& , Geo::Point , FESpace_ptr<N,Q,S>  );
		//PointWiseForce(T const& , Geo::Point , FESpace<N,Q,S> const & );
	    virtual ~PointWiseForce(){};
	    Vec eval(const double &)const;
	private:
	    unsigned int M_startIndex1;
	    unsigned int M_startIndex2;
	    unsigned int M_size;
	    unsigned int M_totSize;
    };
    //template<typename T>
/*    class DiffuseForce : public Force{*/
	//public: 
	    //template<int N,typename Q,typename S>
		////DiffuseForce(T const& , FESpace<N,Q,S> const & );
		//DiffuseForce(std::function<std::array<double,2>(const double &)> const& , FESpace<N,Q,S> const & );
	    //virtual ~PointWiseForce(){};
	    //Vec eval(const double &)const;

    /*};*/
}

#include"Force_imp.hpp"
#endif

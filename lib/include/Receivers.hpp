#ifndef __RECEIVERS_HPP__
#define __RECEIVERS_HPP__ 1
#include<string>
#include"Geometry.hpp"
#include"FESpace.hpp"
#include<fstream>
#include<vector>

namespace Tspeed{
    //template<int N, typename Q, typename S>
    class PointWiseEntity{
	public:
	    virtual ~PointWiseEntity(){};
	    Eigen::ArrayXd const & shape(int i)const{return M_shape[i];};
	    Geo::Point const & point(int i)const{return M_relp[i];};
	    unsigned int const & elem(int i)const{return M_ie[i];};
	    unsigned int size()const{return M_ie.size();};
	protected:
	    std::vector<unsigned int> M_ie;
	    std::vector<Geo::Point> M_relp;
	    std::vector<Eigen::ArrayXd> M_shape;
	    unsigned int M_nel;
	    template<int N, typename Q, typename S>
		void M_add(FESpace_ptr<N,Q,S>, Geo::Point const & );

    };
    class Receivers: public PointWiseEntity
    {
	public:
	    template<int N, typename Q, typename S>
		Receivers(FESpace_ptr<N,Q,S> ,std::string const &);
	    template<int N, typename Q, typename S>
		Receivers(FESpace_ptr<N,Q,S>, Geo::Point const & );
	    void add(double const & x, double const &y, unsigned int const &, unsigned int const &);
	    void write(std::string const &)const;
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

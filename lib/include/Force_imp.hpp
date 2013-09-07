#include"Force.hpp"
namespace Tspeed{
    template<int N,typename Q,typename S>
    PointWiseForce::PointWiseForce(std::function<std::array<double,2>(const double &)>const & f, Geo::Point p, FESpace_ptr<N,Q,S>  Xh):M_totSize(2*Xh->ne()*Xh->nln())
    {
	M_f = f;
	M_add(Xh, p);
	M_startIndex1 = Xh->nln()*(M_ie[0]);
	M_startIndex2 = Xh->nln()*(M_ie[0])+Xh->ne()*Xh->nln();
	M_size = Xh->nln();
    }
    //template<int N,typename Q,typename S>
    //DiffuseForce::DiffuseForce(std::function<std::array<double,2>(const double &, const double &, const double &)>const & f,  FESpace<N,Q,S> const & Xh):this->M_totSize(2*Xh.ne()*Xh.nln())
    //{
	//M_f = f;
    //}	
    
    /*template<int N,typename Q,typename S>*/
    //Eigen::VectorXd PointWiseForce::eval(const double & t, FESpace<N,Q,S> const & Xh)const
    //{
		
	//std::array<double, 2> F = M_f(t);
	//Vec rhs=Vec::Zero(M_totSize);
	//rhs.segment(M_startIndex1,M_size) = F[0] * M_shape[0];
	//rhs.segment(M_startIndex2,M_size) = F[1] * M_shape[0];
	//return rhs;
    /*}*/
   /* template<int N,typename Q,typename S>*/
	//Eigen::VectorXd DiffuseForce::eval(const double & t, FeSpace<N,Q,S> const & Xh)const
	//{
	    //Vec* rhs = new Vec(M_totSize);
	    //for(auto ie: Xh.elements())
	    //{
		//Xh.loc_rhs(ie, )
	    //}	

	/*}*/


}

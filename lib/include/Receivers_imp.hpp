#include"Receivers.hpp"
/*namespace{*/
    //bool intriangle(Tspeed::Geo::Point p, std::array<Tspeed::Geo::Point,3> v)
    //{
	//double s = v[0].y()*v[2].x() - v[0].x()*v[2].y() + (v[2].y() - v[0].y())*p.x() + (v[0].x() - v[2].x())*p.y();
	//double t = v[0].x()*v[1].y() - v[0].y()*v[1].x() + (v[0].y() - v[1].y())*p.x() + (v[1].x() - v[0].x())*p.y();
	//double A = 1/2*(-v[1].y()*v[2].x() + v[0].y()*(-v[1].x() + v[2].x()) + v[0].x()*(v[1].y() - v[2].y()) + v[1].x()*v[2].y());
	//return (s>-1e-6 && t>-1e-6 && (s+t)<(2+1e-6)*A); //tolerance for edge (sadly needed)
    //}

/*}*/
namespace Tspeed{

    template<int N, typename Q, typename S>
    void PointWiseEntity::M_add(FESpace_ptr<N,Q,S>  Xh, Geo::Point const & p)
    {
	Eigen::ArrayXd temp_shape(Xh->nln());
	for(auto ie : Xh->mesh()->elements())
	{
	    if(ie.intriangle(p))
	    {
		M_ie.push_back(ie.id());
		M_relp.push_back(ie.invmap(p));
		for(unsigned int s=0; s<Xh->nln(); ++s)
		    temp_shape[s] = Xh->shape().phi(s,M_relp.back().x(), M_relp.back().y());
		M_shape.push_back(temp_shape);
		return;
	    }
	}
	std::cerr << "Point " << p << " not found. Exiting" << std::endl; exit(1);
    }


    template<int N, typename Q, typename S>
    Receivers::Receivers(FESpace_ptr<N,Q,S>  Xh, Geo::Point const & p)
    {
	M_add(Xh, p);
	M_nel=1;

    }
    template<int N, typename Q, typename S>
    Receivers::Receivers(FESpace_ptr<N,Q,S>  Xh, std::string const & fname)
    {
	std::ifstream f;
	f.open(fname.c_str());
	double x,y;
	unsigned int count=0;
	//Eigen::ArrayXd temp_shape(Xh->nln());
	while(f>>x>>y)
	{
	    Geo::Point p(x,y);
	    M_add(Xh, p);
	    ++count;
	}
	f.close();
	std::cout << "Receivers :: added " << count << " receivers." << std::endl;
	M_nel = count;
    }



}

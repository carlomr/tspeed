#include "Geometry.hpp"
namespace Tspeed
{
    namespace Geo
    {

	std::ostream & operator<<(std::ostream& io, Point const & p)
	{
	    io << "(" << p.x()<<","<<p.y()<<")";
	    return io;
	}

	Point operator -(const Point & a, const Point & b)
	{
	    return Point(a.M_coord[0]-b.M_coord[0],
		    a.M_coord[1]-b.M_coord[1]);
	}

	Point operator -(const Eigen::Vector2d& a, const Point & b)
	{
	    return Point(a(0)-b.M_coord[0],
		    a(1)-b.M_coord[1]);
	}
	Point operator -(const Point & a, const Eigen::Vector2d& b)
	{
	    return Point(a.M_coord[0]-b(0),
		    a.M_coord[1]-b(1));
	}
	Point operator +(const Eigen::Vector2d& a, const Point & b)
	{
	    return Point(a(0)+b.M_coord[0],
		    a(1)+b.M_coord[1]);
	}
	Point operator +(const Point & a, const Eigen::Vector2d& b)
	{
	    return Point(a.M_coord[0]+b(0),
		    a.M_coord[1]+b(1));
	}

	Point operator +(const Point & a, const Point & b)
	{
	    return Point(a.M_coord[0]+b.M_coord[0],
		    a.M_coord[1]+b.M_coord[1]);
	}

	Point Point::operator *(const double & d) const
	{
	    return Point(d*M_coord[0],d*M_coord[1]);
	}

	Point & Point::operator =(const Point & p)
	{
	    if(this!=&p)
	    {
		M_coord[0]=p.M_coord[0];
		M_coord[1]=p.M_coord[1];
		M_id = p.M_id;
		M_reg = p.M_reg;
		M_bcId = p.M_bcId;
	    }
	    return *this;
	}

	Point operator *(const double & d, const Point & p)
	{
	    return p*d;
	}
	////////////////////////////////////////////////////////////
	Edge & Edge::operator=(Edge const & e)
	{
	    if(this!=&e)
	    {
		M_extr[0] = e.M_extr[0];
		M_extr[1] = e.M_extr[1];
		M_bcId = e.M_bcId;
		M_id = e.M_id;
	    }
	    return *this;
	}
	Eigen::Vector2d Edge::normal()const{
	    Eigen::Vector2d nn;
	    nn << -M_extr[0].y()+M_extr[1].y(), M_extr[0].x() - M_extr[1].x();
	    double norm = nn.norm();
	    return nn/norm;
	}
	////////////////////////////////////////////////////////////
	//Triangle::Triangle() 
	//{
	    //M_points[0] = M_points[1] = M_points[2] = nullptr;
	    //M_edges[0] = M_edges[1] = M_edges[2] = nullptr;
	//}
	bool Triangle::intriangle(const Point & p)const
	{
	    double s = M_points[0].y()*M_points[2].x() - M_points[0].x()*M_points[2].y() + (M_points[2].y() - M_points[0].y())*p.x() + (M_points[0].x() - M_points[2].x())*p.y();
	    double t = M_points[0].x()*M_points[1].y() - M_points[0].y()*M_points[1].x() + (M_points[0].y() - M_points[1].y())*p.x() + (M_points[1].x() - M_points[0].x())*p.y();
	    double A = 1./2.*(-M_points[1].y()*M_points[2].x() + M_points[0].y()*(-M_points[1].x() + M_points[2].x()) + M_points[0].x()*(M_points[1].y() - M_points[2].y()) + M_points[1].x()*M_points[2].y());
	    return (s>-1e-6 && t>-1e-6 && (s+t)<(2+1e-6)*A); //tolerance for edge (sadly needed)
	}
	Triangle::Triangle() 
	{
	    M_points[0] = M_points[1] = M_points[2] = Point(); 
	    M_edges[0] = M_edges[1] = M_edges[2] = Edge();
	}

	void Triangle::M_set_edges()
	{
	    M_edges[0] = Edge(M_points[1],M_points[2]);
	    M_edges[1] = Edge(M_points[2],M_points[0]);
	    M_edges[2] = Edge(M_points[0],M_points[1]);
	    for (auto i: {0,1,2})
		M_edges[i].id() = i;
	}


	Triangle::Triangle( const Point & a, const Point & b, const Point & c ):Entity()
	{
	    M_points[0] = a;
	    M_points[1] = b;
	    M_points[2] = c;
	    M_set_edges();
	    M_neigh = {{-1,-1,-1}};
	    M_neighedges = {{-1,-1,-1}};
	    setJac();
	    setinvJac();
	    setdetJ();
	}

	/*Triangle::Triangle( const Triangle & t )*/
	//{
	    //M_points[0] = t.M_points[0];
	    //M_points[1] = t.M_points[1];
	    //M_points[2] = t.M_points[2];
	    //M_edges[0]  = t.M_edges[0];
	    //M_edges[1]  = t.M_edges[1];
	    //M_edges[2]  = t.M_edges[2];

	    //M_id     = t.M_id;
	    //M_reg    = t.M_reg;
	    //M_bcId   = t.M_bcId;
	    //M_jac    = t.M_jac;
	    //M_invjac = t.M_invjac;
	    //M_detJ   = t.M_detJ;
	/*}*/

	Triangle& Triangle::operator =(const Triangle & t)
	{
	    if(this!=&t)
	    {
		M_points[0] = t.M_points[0];
		M_points[1] = t.M_points[1];
		M_points[2] = t.M_points[2];
		M_edges[0] = t.M_edges[0];
		M_edges[1] = t.M_edges[1];
		M_edges[2] = t.M_edges[2];
	    }
	    return *this;
	}
	void Triangle::setJac(){
	    Eigen::Matrix2d J;
	    J.col(0) = (M_points[1]-M_points[0]).toEig();
	    J.col(1) = (M_points[2]-M_points[0]).toEig();
	    M_jac = J;
	}
	Eigen::Matrix2d Triangle::Jac()const{
	    return M_jac;
	}
	double Triangle::detJ()const{
	    return M_detJ;
	}
	void Triangle::setdetJ(){
	    M_detJ = Jac().determinant();
	}
	Eigen::Matrix2d Triangle::invJac()const{
	    return M_invjac;
	}
	void Triangle::setinvJac(){
	    M_invjac = M_jac.inverse();
	}
	Point Triangle::invmap(Point const & p)const
	{
	    return Point(M_invjac*(p -M_points[0]).toEig());
	}
	Point Triangle::map(Point const & p)const
	{
	    Point ph;
	    return ph = M_points[0] + M_jac*p.toEig();
	}
	std::ostream & operator<<(std::ostream & io, const Triangle & t)
	{
	    for(auto i : t.all_pts())
		io << i << " ";
	    return io;
	}



    }
}

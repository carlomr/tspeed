#ifndef __GEOMETRY_HPP__
#define  __GEOMETRY_HPP__ 1

#include<array>
#include<cmath>
#include<Eigen/Dense>
#include<memory>
#include<limits>
#include<iostream>
const unsigned int NVAL=std::numeric_limits<unsigned int>::max();

namespace Tspeed
{
    enum class Bc {Dirichlet, Neumann, Internal, Unassigned};
    class Entity
    {
	public:
	    typedef unsigned int Id;
	    Entity():M_reg(NVAL),M_id(NVAL),M_bcId(Bc::Unassigned){}
	    bool unassignedId()const {return M_id==NVAL;}
	    bool unassignedBc()const {return M_bcId==Bc::Unassigned;}
	    bool unassignedReg()const {return M_reg==NVAL;}
	    //void setid(Id const i){M_id=i;};
	    //void setbcId(Bc const i){M_bcId=i;};
	    //void setreg(Id reg) {M_reg= reg;};
	    Id & reg(){return M_reg;}
	    Id const & reg()const {return M_reg;}
	    Id & id(){return M_id;}
	    Id const & id()const {return M_id;}
	    Bc & bcId(){return M_bcId;}
	    Bc const & bcId()const {return M_bcId;}
	protected:
	    Id M_reg;
	    Id M_id;
	    Bc M_bcId;
    };

    namespace Geo
    {
	using std::shared_ptr;
	static const int dim = 2;

	class Point: public Entity{
	    public:
		Point(const double& x=0, const double& y=0):Entity(){M_coord[0]=x;M_coord[1]=y;};
		Point(const Point & p){M_coord[0] = p.M_coord[0]; M_coord[1] = p.M_coord[1]; M_id = p.M_id; M_reg = p.M_reg; M_bcId = p.M_bcId;};
		Point(const Eigen::Vector2d& v){M_coord[0] = v(0); M_coord[1] = v(1);};
		virtual ~Point(){};
		double x()const{return M_coord[0];};
		double y()const{return M_coord[1];};
		double& x(){return M_coord[0];};
		double& y(){return M_coord[1];};
		friend Point operator+(const Point&, const Point&);
		friend Point operator+(const Eigen::Vector2d&, const Point&);
		friend Point operator+(const Point&,const Eigen::Vector2d&);
		friend Point operator-(const Point&, const Point&);
		friend Point operator-(const Eigen::Vector2d&, const Point&);
		friend Point operator-(const Point&,const Eigen::Vector2d&);
		friend Point operator*(const double &, const Point &);
		Point & operator=(const Point&);
		Point operator *(const double &)const;
		friend double dot(const Point & a, const Point &b){return a.M_coord[0]*b.M_coord[0] + a.M_coord[1]*b.M_coord[1]; };
		double norm()const{return std::sqrt(M_coord[0]*M_coord[0]+ M_coord[1]*M_coord[1]);};
		Eigen::Vector2d data()const{return M_coord;}
		//Eigen::Vector2d eig_v()const{Eigen::Vector2d v; v(0) = M_coord[0]; v(1) = M_coord[1]; return v;};
	    private:
		//std::array<double,dim> M_coord;
		Eigen::Vector2d M_coord;
	};

	class Edge : public Entity
	{
	    public:
		Edge():Entity(){M_extr[0] = M_extr[1] = Point();};
		//Edge(){M_extr[0] = M_extr[1] = nullptr;};
		Edge(const Point& a, const Point& b):Entity(){M_extr[0] = a; M_extr[1]= b;};
		//Edge(Point& a, Point& b){M_extr[0] = shared_ptr<Point>(&a); M_extr[1]= shared_ptr<Point>(&b);};
		Edge(const Edge& e){ M_extr[1] = e.M_extr[1]; M_extr[0] = e.M_extr[0]; M_bcId = e.M_bcId; M_id = e.M_id;};
		virtual ~Edge(){};
		double length()const{return (M_extr[0]-M_extr[1]).norm();};
		Eigen::Vector2d normal()const;
		//double length()const{return (M_extr[0]-M_extr[1]).norm();};
		Edge& operator=(const Edge&);
	    private:
		std::array<Point, 2> M_extr;
		//std::array<shared_ptr<Point>, 2> M_extr;

	};


	class Triangle : public Entity
	{
	    public:
		static const int numVertices=3;
		Triangle(); //Constructs an empty triangle
		Triangle(  const Point&,  const Point&,  const Point& ); //Points are given (by reference)
		//Triangle( const Point&,  const Point&,  const Point&, const Edge&, const Edge&,const  Edge& ); //Points are given (by reference)
		Triangle( const Triangle& ) = default;
		Triangle & operator=( const Triangle& );
		virtual ~Triangle(){};
		//shared_ptr<Point>  const & pt( int i ) const { return M_points[i]; }
		//shared_ptr<Edge> const & edg(int i) const{return M_edges[i];};
		std::array<Point,3> all_pts()const{return M_points;}
		std::array<Edge,3> all_edges()const{return M_edges;}
		Point const & pt( int i ) const { return M_points[i]; }
		Edge const & edg(int i) const{return M_edges[i];};
		Eigen::Matrix2d Jac()const;
		Eigen::Matrix2d invJac()const;
		double detJ()const;
		Point map(Point const & p)const;
		Point invmap(Point const & p)const;
		int const & neigh(int i)const{return M_neigh[i];};
		int const & neighedges(int i)const{return M_neighedges[i];};
		void setNeigh(int i, int j){M_neigh[i]=j;}
		void setNeighedges(int i, int j){M_neighedges[i]=j; M_edges[i].bcId() = Bc::Internal;}
		void printNeigh()const{std::cout << M_id << " -> " << M_neigh[0] << '\t' << M_neigh[1] << '\t' << M_neigh[2] << std::endl;};
		bool intriangle(const Point &)const;
		//std::shared_ptr<Triangle> neigh(unsigned short int i){return M_neigh[i];}


	    private:
		//std::array<shared_ptr<Point>,3>  M_points;
		//std::array<shared_ptr<Edge> ,3>   M_edges;
		void setJac();
		void setinvJac();
		void setdetJ();
		std::array<Point,3>  M_points;
		std::array<Edge,3>   M_edges;
		std::array<int,3> M_neigh;
		std::array<int, 3> M_neighedges;
		Eigen::Matrix2d M_jac;
		Eigen::Matrix2d M_invjac;
		double M_detJ;
		void M_set_edges();

	};
	std::ostream & operator<<(std::ostream &, Triangle const & );
	std::ostream & operator<<(std::ostream &, Point const &);
    }

}


#endif

/**
 * @file Geometry.hpp
 * @brief Header file for the geometrical entities
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
    /**
     * @brief Base class for geometrical entities
     */
    class Entity
    {
	public:
	    typedef unsigned int Id;
	    Entity():M_reg(NVAL),M_id(NVAL),M_bcId(Bc::Unassigned){}
	    /**
	     * @brief Check if element has assigned Id
	     *
	     * @return TRUE if Id is not assigned
	     */
	    bool unassignedId()const {return M_id==NVAL;}
	    /**
	     * @brief Check if element has assigned boundary condition
	     *
	     * @return TRUE if Bc is not assigned
	     */
	    bool unassignedBc()const {return M_bcId==Bc::Unassigned;}
	    /**
	     * @brief Check if element has assigned region
	     *
	     * @return TRUE if region is not assigned
	     */
	    bool unassignedReg()const {return M_reg==NVAL;}
	    //void setid(Id const i){M_id=i;};
	    //void setbcId(Bc const i){M_bcId=i;};
	    //void setreg(Id reg) {M_reg= reg;};
	    /**
	     * @brief Get region id
	     *
	     * @return Reference to region id
	     */
	    Id & reg(){return M_reg;}
	    /**
	     * @brief Get region id (const version)
	     *
	     * @return Constant reference to region id
	     */
	    Id const & reg()const {return M_reg;}
	    /**
	     * @brief Get element id (const version)
	     *
	     * @return Constant reference to element id
	     */
	    Id & id(){return M_id;}
	    /**
	     * @brief Get element id (const version)
	     *
	     * @return Constant reference to element id
	     */
	    Id const & id()const {return M_id;}
	    /**
	     * @brief Get boundary id (const version)
	     *
	     * @return Constant reference to boundary id
	     */
	    Bc & bcId(){return M_bcId;}
	    /**
	     * @brief Get boundary id (const version)
	     *
	     * @return Constant reference to boundary id
	     */
	    Bc const & bcId()const {return M_bcId;}
	protected:
	    /**
	     * @brief Region number
	     */
	    Id M_reg;
	    /**
	     * @brief Id of the element
	     */
	    Id M_id;
	    /**
	     * @brief Boundary Id of the element
	     */
	    Bc M_bcId;
    };

    namespace Geo
    {
	using std::shared_ptr;
	static const int dim = 2;

	/**
	 * @brief Class describing points
	 */
	class Point: public Entity{
	    public:
		/**
		 * @brief Constructor taking the two coordinates of the point. By default, a point is initialized to the origin
		 *
		 * @param x x-coordinate
		 * @param y y-coordinate
		 */
		Point(const double& x=0, const double& y=0):Entity(){M_coord[0]=x;M_coord[1]=y;};
		/**
		 * @brief Copy constructor. Everything is copied
		 *
		 * @param p A point
		 */
		Point(const Point & p){M_coord[0] = p.M_coord[0]; M_coord[1] = p.M_coord[1]; M_id = p.M_id; M_reg = p.M_reg; M_bcId = p.M_bcId;};
		/**
		 * @brief Constructor taking an Eigen fixed size vector. The components are copied, ids are not assigned
		 *
		 * @param v
		 */
		Point(const Eigen::Vector2d& v){M_coord[0] = v(0); M_coord[1] = v(1);};
		virtual ~Point(){};
		/**
		 * @brief first coordinate
		 *
		 * @return x-coordinate
		 */
		double x()const{return M_coord[0];};
		/**
		 * @brief second coordinate
		 *
		 * @return y-coordinate
		 */
		double y()const{return M_coord[1];};
		/**
		 * @brief Reference to first coordinate
		 *
		 * @return Reference to x-coord
		 */
		double& x(){return M_coord[0];};
		/**
		 * @brief Reference to second coordinate
		 *
		 * @return Reference to y-coord
		 */
		double& y(){return M_coord[1];};
		/**
		 * @brief Operator summing two points
		 *
		 * @param a first point
		 * @param b second point
		 *
		 * @return sum of points
		 */
		friend Point operator+(const Point& a, const Point& b);
		/**
		 * @brief Sum a point and en eigen vector
		 *
		 * @param a vector
		 * @param b Point
		 *
		 * @return Point
		 */
		friend Point operator+(const Eigen::Vector2d& a, const Point& b);
		/**
		 * @brief Sum a point and en eigen vector
		 *
		 * @param a Point
		 * @param b vector
		 *
		 * @return Point
		 */
		friend Point operator+(const Point& a,const Eigen::Vector2d& b);
		/**
		 * @brief Operator sutracting two points
		 *
		 * @param a first point
		 * @param b second point
		 *
		 * @return difference of points
		 */
		friend Point operator-(const Point& a, const Point& b);
		/**
		 * @brief Operator sutracting a point and a vector
		 *
		 * @param a first point
		 * @param b vector
		 *
		 * @return point: difference of point a vector
		 */
		friend Point operator-(const Eigen::Vector2d& a, const Point& b);
		friend Point operator-(const Point&,const Eigen::Vector2d&);
		friend Point operator*(const double &, const Point &);
		/**
		 * @brief assignemnt operator
		 *
		 */
		Point & operator=(const Point&);
		/**
		 * @brief Multuply a point by a scalar	
		 *
		 * @param a scalar
		 *
		 * @return a both with both coordinates multiplied
		 */
		Point operator *(const double & a)const;
		/**
		 * @brief vector-style dot product between points
		 *
		 * @param a first point
		 * @param b second point
		 *
		 * @return a scalar
		 */
		friend double dot(const Point & a, const Point &b){return a.M_coord[0]*b.M_coord[0] + a.M_coord[1]*b.M_coord[1]; };
		/**
		 * @brief Norm of the vector from the origin to the point
		 *
		 * @return the euler norm
		 */
		double norm()const{return std::sqrt(M_coord[0]*M_coord[0]+ M_coord[1]*M_coord[1]);};
		//Eigen::Vector2d data()const{return M_coord;}
		/**
		 * @brief convert the point into an Eigen::vector. Useful for matrix tranformations
		 *
		 * @return the eigen vector with the coordinates as components
		 */
		Eigen::Vector2d toEig()const{Eigen::Vector2d v; v << M_coord[0], M_coord[1]; return v;};
		//Eigen::Vector2d eig_v()const{Eigen::Vector2d v; v(0) = M_coord[0]; v(1) = M_coord[1]; return v;};
	    private:
		std::array<double,dim> M_coord;
		//Eigen::Vector2d M_coord;
	};

	/**
	 * @brief Class describing an edge 
	 */
	class Edge : public Entity
	{
	    public:
		/**
		 * @brief Default constructor: extremal points are initialized to 0
		 *
		 */
		Edge():Entity(){M_extr[0] = M_extr[1] = Point();};
		//Edge(){M_extr[0] = M_extr[1] = nullptr;};
		/**
		 * @brief Constructor
		 *
		 * @param a One extreme
		 * @param b The other extreme
		 */
		Edge(const Point& a, const Point& b):Entity(){M_extr[0] = a; M_extr[1]= b;};
		//Edge(Point& a, Point& b){M_extr[0] = shared_ptr<Point>(&a); M_extr[1]= shared_ptr<Point>(&b);};
		/**
		 * @brief copy constructor
		 *
		 * @param e
		 */
		Edge(const Edge& e){ M_extr[1] = e.M_extr[1]; M_extr[0] = e.M_extr[0]; M_bcId = e.M_bcId; M_id = e.M_id;};
		virtual ~Edge(){};
		/**
		 * @brief Length of the edge
		 *
		 */
		double length()const{return (M_extr[0]-M_extr[1]).norm();};
		/**
		 * @brief Normal to the edge
		 *
		 * @return The normalized vector
		 */
		Eigen::Vector2d normal()const;
		//double length()const{return (M_extr[0]-M_extr[1]).norm();};
		/**
		 * @brief Assignement
		 *
		 *
		 */
		Edge& operator=(const Edge&);
	    private:
		std::array<Point, 2> M_extr;
		//std::array<shared_ptr<Point>, 2> M_extr;

	};


	/**
	 * @brief Class describing a triangle
	 */
	class Triangle : public Entity
	{
	    public:
		static const int numVertices=3;
		Triangle(); //Constructs an empty triangle
		/**
		 * @brief Create a triangle from three points
		 *
		 * @param a,b,c The three points
		 */
		Triangle(  const Point& a,  const Point& b,  const Point& c); //Points are given (by reference)
		//Triangle( const Point&,  const Point&,  const Point&, const Edge&, const Edge&,const  Edge& ); //Points are given (by reference)
		/**
		 * @brief Copy constructor
		 *
		 */
		Triangle( const Triangle& ) = default;
		/**
		 * @brief Assignement 
		 *
		 */
		Triangle & operator=( const Triangle& );
		virtual ~Triangle(){};
		//shared_ptr<Point>  const & pt( int i ) const { return M_points[i]; }
		//shared_ptr<Edge> const & edg(int i) const{return M_edges[i];};
		/**
		 * @brief Get all points of the triangle
		 *
		 * @return An array of three points
		 */
		std::array<Point,3> all_pts()const{return M_points;}
		/**
		 * @brief Get all edges of a triangle
		 *
		 * @return An array of the three edge
		 */
		std::array<Edge,3> all_edges()const{return M_edges;}
		/**
		 * @brief Get a point
		 *
		 * @param i Number of the point in the triangle
		 *
		 * @return The i-th point
		 */
		Point const & pt( int i ) const { return M_points[i]; }
		/**
		 * @brief Get a edge
		 *
		 * @param i Number of the edge in the triangle
		 *
		 * @return The i-th edge
		 */
		Edge const & edg(int i) const{return M_edges[i];};
		/**
		 * @brief Get Jacobian of the transformation from the reference triangle
		 *
		 * @return The Jacobian, in matrix form
		 */
		Eigen::Matrix2d Jac()const;
		/**
		 * @brief Get the inverse Jacobian of the transformation from the reference triangle
		 *
		 * @return The inverse Jacobian, in matrix form
		 */
		Eigen::Matrix2d invJac()const;
		/**
		 * @brief Get the determinant of the Jacobian of the tranformation from the reference triangle
		 *
		 * @return the determinant of the Jacobian (i.e., Area(T)*2)
		 */
		double detJ()const;
		/**
		 * @brief Map a point from its relative position in the reference triangle to the physical point
		 *
		 * @param p The point in the reference triangle
		 *
		 * @return The physical point in the actual triangle
		 */
		Point map(Point const & p)const;
		/**
		 * @brief Map a point from the physical triangle to the reference one
		 *
		 * @param p The physical point
		 *
		 * @return The point in the reference triangle
		 */
		Point invmap(Point const & p)const;
		/**
		 * @brief Ged index of neighboring triangle on edge i
		 *
		 * @param i Edge of the triangle
		 *
		 * @return Index of the neighbor
		 */
		int const & neigh(int i)const{return M_neigh[i];};
		/**
		 * @brief Index of the edge i in the nieghboring triangle
		 *
		 * @param i Edge of the present triangle
		 *
		 * @return Index of the edge in the neighboring triangle
		 */
		int const & neighedges(int i)const{return M_neighedges[i];};
		/**
		 * @brief Set neighboring triangle
		 *
		 * @param i Edge of the current triangle
		 * @param j Index of the neighbor
		 */
		void setNeigh(int i, int j){M_neigh[i]=j;}
		/**
		 * @brief Set index of the edge of the nieghboring triangle
		 *
		 * @param i Edge of the current trinagle
		 * @param j Edge in the neighboring triangle
		 */
		void setNeighedges(int i, int j){M_neighedges[i]=j; M_edges[i].bcId() = Bc::Internal;}
		/**
		 * @brief Pirnt neighbors for the current triangle
		 */
		void printNeigh()const{std::cout << M_id << " -> " << M_neigh[0] << '\t' << M_neigh[1] << '\t' << M_neigh[2] << std::endl;};
		/**
		 * @brief Check if point p is in triangle
		 *
		 * @param p The point
		 *
		 * @return TRUE if the point is in the triangle
		 */
		bool intriangle(const Point & p)const;
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

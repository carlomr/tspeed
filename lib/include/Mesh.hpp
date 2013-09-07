#ifndef __MESH_HPP__
#define __MESH_HPP__ 1
#include<string>
#include<fstream>
#include<iostream>
#include<vector>
#include<algorithm>
#include<map>
#include<Eigen/StdVector>
#include"Geometry.hpp"

namespace Tspeed
{
    class Mesh
    {	
	public:
	    typedef unsigned int size_type;
	    explicit Mesh(const std::string);
	    Geo::Triangle const & operator[](size_t i)const{return M_tri[i];}
	    Geo::Triangle & operator[](size_t i){return M_tri[i];}
	    //std::vector<Geo::Triangle> get_tria()const{return M_tri;}
	    //int connectivity(const unsigned int ie,const unsigned int iedg)const{return M_conn[ie*3+iedg];}
	    //Mesh(const Mesh &);
	    //Mesh & operator =(Mesh const &);
	    std::vector<Geo::Triangle> const & elements()const{return M_tri;}
	    std::vector<Geo::Triangle> & elements(){return M_tri;}
	    ~Mesh(){};
	    void stats()const;
	    unsigned int ne()const{return M_tri.size();};
	    void printallNeigh()const{for(auto& ie:M_tri) ie.printNeigh();};
	    //friend std::ostream & operator<<(std::ostream &, Mesh const&);
	    std::vector<std::pair<unsigned int,unsigned int>> M_bed_map;
	    
	private:
	    std::map<std::string, Bc> M_BcMap;
	    std::map<Bc,std::string> M_BcMap_inv;
	    std::vector<Geo::Triangle> M_tri;
	    std::vector<Geo::Point> M_points;
	    std::vector<Geo::Edge> M_bed;
	    void M_neighbors();
	    std::string  M_fname;
	    //int M_verbose;
	    //std::vector<Geo::Edges> M_edg;
	    //std::vector<int> M_conn;
    };
    typedef std::shared_ptr<Mesh> Mesh_ptr;

}


#endif


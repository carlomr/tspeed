/**
 * @file Mesh.hpp
 * @brief Header file for the mesh
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
#ifndef __MESH_HPP__
#define __MESH_HPP__ 1
#include<string>
#include<fstream>
#include<iostream>
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
	    typedef std::vector<Geo::Triangle,Eigen::aligned_allocator<Eigen::Vector2d>> AlignedVecT;
	    typedef std::vector<Geo::Edge,Eigen::aligned_allocator<Eigen::Vector2d>> AlignedVecE;
	    typedef std::vector<Geo::Point,Eigen::aligned_allocator<Eigen::Vector2d>> AlignedVecP;
	    /**
	     * @brief Generate mesh from a Gmsh mesh
	     *
	     * @param fileName Name of the .msh file containing the mesh
	     */
	    explicit Mesh(const std::string fileName);
	    /**
	     * @brief Get triangle with index i (const)
	     *
	     * @param i index of the triangle
	     *
	     * @return Constant reference to triangle
	     */
	    Geo::Triangle const & operator[](size_t i)const{return M_tri[i];}
	    /**
	     * @brief Get triangle with index i ( non-const)
	     *
	     * @param i index of the triangle
	     *
	     * @return Reference to triangle
	     */
	    Geo::Triangle & operator[](size_t i){return M_tri[i];}
	    //std::vector<Geo::Triangle> get_tria()const{return M_tri;}
	    //int connectivity(const unsigned int ie,const unsigned int iedg)const{return M_conn[ie*3+iedg];}
	    //Mesh(const Mesh &);
	    //Mesh & operator =(Mesh const &);
	    /**
	     * @brief Get all triangles in the mesh (const)
	     *
	     * @return Constant reference to a vector of triangles
	     */
	    AlignedVecT const & elements()const{return M_tri;}
	    /**
	     * @brief Get all triangles in the mesh (non-const)
	     *
	     * @return Reference to a vector of triangles
	     */
	    AlignedVecT & elements(){return M_tri;}
	    ~Mesh(){};
	    /**
	     * @brief Print stats about the mesh (e.g. number of elements, anisotropy etc.)
	     */
	    void stats()const;
	    /**
	     * @brief Get number of triangles in the mesh
	     *
	     */
	    unsigned int ne()const{return M_tri.size();};
	    /**
	     * @brief Print neighbors structure
	     */
	    void printallNeigh()const{for(auto& ie:M_tri) ie.printNeigh();};
	    //friend std::ostream & operator<<(std::ostream &, Mesh const&);
	    std::vector<std::pair<unsigned int,unsigned int>> M_bed_map;
	    
	private:
	    std::map<std::string, Bc> M_BcMap;
	    std::map<Bc,std::string> M_BcMap_inv;
	    AlignedVecT M_tri;
	    AlignedVecP M_points;
	    AlignedVecE M_bed;
	    void M_neighbors();
	    std::string  M_fname;
	    //int M_verbose;
	    //std::vector<Geo::Edges> M_edg;
	    //std::vector<int> M_conn;
    };
    /**
     * @brief Shared pointer to an element of type mesh
     */
    typedef std::shared_ptr<Mesh> Mesh_ptr;

}


#endif


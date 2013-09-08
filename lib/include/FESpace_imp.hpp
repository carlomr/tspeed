/**
 * @file FESpace_imp.hpp
 * @brief Implementation of the functional space class methods
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
#ifndef __FESPACE_IMP_HPP__
#define __FESPACE_IMP_HPP__ 1
namespace Tspeed{


    template<int N, typename Q, typename S>
	FESpace<N,Q,S>::FESpace(Mesh_ptr m):M_mesh(m),M_quad(),M_sh(),M_nln((N+1)*(N+2)/2),M_nqn2d(Q::nqn2d),M_ne(M_mesh->ne())
    {
	for(int i=0; i<M_nln; ++i)
	{
	    M_B.col(i)    = M_sh.phi(i,M_quad.int_nodes().col(0), M_quad.int_nodes().col(1));
	    M_G[i].row(0) = M_sh.grad(i,M_quad.int_nodes().col(0), M_quad.int_nodes().col(1)).col(0).transpose();
	    M_G[i].row(1) = M_sh.grad(i,M_quad.int_nodes().col(0), M_quad.int_nodes().col(1)).col(1).transpose();
	    M_vert_basis(i,0) = M_sh.phi(i,0.2,0.2);
	    M_vert_basis(i,1) = M_sh.phi(i,0.8,0.2);
	    M_vert_basis(i,2) = M_sh.phi(i,0.2,0.8);
	    M_vert_basis(i,3) = M_sh.phi(i,0.5,0.5);
	    for(int iedg=0; iedg<3; ++iedg)
	    {
		M_Bedge[iedg].col(i)   = M_sh.phi (i, M_quad.edge_nodes(iedg).col(0), M_quad.edge_nodes(iedg).col(1));
		M_Gedge[i*3+iedg].row(0)   = M_sh.grad(i, M_quad.edge_nodes(iedg).col(0), M_quad.edge_nodes(iedg).col(1)).col(0);
		M_Gedge[i*3+iedg].row(1) = M_sh.grad(i, M_quad.edge_nodes(iedg).col(0), M_quad.edge_nodes(iedg).col(1)).col(1);
	    }
	}
	base_mass_and_inv();
	std::cout << "FESpace :: " << (N+1)*(N+2)/2*M_mesh->ne() << " dof per component." << std::endl;

    };
    template<int N, typename Q, typename S>
	void FESpace<N,Q,S>::base_mass_and_inv()
	{
	    if(S::is_orthonormal)
	    {
		M_base_mass = Eigen::MatrixXd::Identity(S::gdl, S::gdl);
	    }
	    else
	    {
		M_base_mass = Eigen::MatrixXd::Zero(S::gdl, S::gdl);
		for(int i=0; i<S::gdl ; ++i)
		{
		    for (int j = 0; j<S::gdl ; ++j)
		    {   
			for(int k = 0; k<Q::nqn2d ; ++k)
			    M_base_mass(i,j) += M_B(k,i)*M_B(k,j)*M_quad.iweight(k);
		    }
		}
	    }
	    M_base_invmass = M_base_mass.inverse();
	};

    template<int N, typename Q, typename S>
	void FESpace<N,Q,S>::points_out(std::string const &fname)const
	{
	    std::ofstream outf(fname);
	    Geo::Point p0(0.2,0.2),p1(0.8,0.2),p2(0.2,0.8),p3(0.5,0.5);

	    for(auto ie: M_mesh->elements())
	    {
		outf << ie.map(p0).x() << " " << ie.map(p0).y() << std::endl;
		outf << ie.map(p1).x() << " " << ie.map(p1).y() << std::endl;
		outf << ie.map(p2).x() << " " << ie.map(p2).y() << std::endl;
		outf << ie.map(p3).x() << " " << ie.map(p3).y() << std::endl;
	    }
	    outf.close();
	}
    template<int N, typename Q, typename S>
	void FESpace<N,Q,S>::field_out(std::string const &fname, Eigen::VectorXd const &uh, unsigned int step)const
	{
	    //std::ofstream outf_s1(fname+"u_x_"+std::to_string(step)+".field");
	    //std::ofstream outf_s2(fname+"u_y_"+std::to_string(step)+".field");
	    std::ofstream outf_v(fname+"u_"+std::to_string(step)+".field");
	    double valx, valy;
	    for(auto ie: M_mesh->elements())
	    {
		for(int i = 0; i<4 ; ++i)
		{
		    valx = M_vert_basis.col(i).dot(uh.segment(ie.id()*M_nln,M_nln));
		    valy = M_vert_basis.col(i).dot(uh.segment(ie.id()*M_nln+M_ne*M_nln,M_nln));
		    //outf_s1 << valx << std::endl;
		    //outf_s2 << valy << std::endl;
		    outf_v << valx << " " << valy << std::endl;
		}
	    }
	    //outf_s1.close();
	    //outf_s2.close();
	    outf_v.close();
	}

    template<int N, typename Q, typename S>
	double FESpace<N,Q,S>::L2error(std::function<std::array<double,2>(double,double)> const & uex, Eigen::VectorXd const & uh)const
	{
	    Eigen::MatrixXd local_exact = Eigen::MatrixXd::Zero(M_nqn2d,2); 
	    Eigen::MatrixXd local_aprox = Eigen::MatrixXd::Zero(M_nqn2d,2); 
	    std::array<double,2 > F;
	    Geo::Point p;
	    double error = 0;
	    double error_loc = 0;
	    for(auto ie: M_mesh->elements())	
	    {
		error_loc = 0;
		local_aprox.col(0) = M_B*uh.segment(ie.id()*M_nln, M_nln);
		local_aprox.col(1) = M_B*uh.segment(M_nln*M_ne+ie.id()*M_nln, M_nln);
		for(int i=0; i<M_nqn2d; ++i)
		{
		    p = ie.map(M_quad.inode(i));
		    F = uex(p.x(), p.y());
		    //local_exact(i,0) = F[0];
		    //local_exact(i,1) = F[1];
		    error_loc += ((local_aprox(i,0) - F[0])*(local_aprox(i,0)-F[0])+(local_aprox(i,1) - F[1])*(local_aprox(i,1)-F[1]))*M_quad.iweight(i);
		}
		//std::cout << "********" << error_loc << std::endl;
		error +=  error_loc*ie.detJ();
	    }
	    return std::sqrt(error);

	};
    template<int N, typename Q, typename S>
	Eigen::VectorXd FESpace<N,Q,S>::inverse_transform(std::function<std::array<double,2>(double,double)> const & fun)const
	{
	    Eigen::VectorXd uh(2*M_ne*M_nln);
	    unsigned int startIndex;
	    unsigned int size = M_nln;
	    Eigen::VectorXd rhs(M_nln);
	    for(auto ie: M_mesh->elements())
	    {
		startIndex = ie.id()*M_nln;
		rhs = M_loc_rhs(ie,fun);
		//std::cout << ie.id() << ": " << std::endl;
		//std::cout << ie.Jac() << std::endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
		//std::cout << rhs << std::endl;
		//std::cout << std::endl;
		//invM = 1./ie.detJ()*Eigen::MatrixXd::Identity(M_nln,M_nln);
		uh.segment(startIndex, size) = 1./ie.detJ()*M_base_invmass * rhs.segment(0, size);
		uh.segment(startIndex + M_ne*M_nln, size) = 1./ie.detJ()*M_base_invmass * rhs.segment(M_nln, size);


	    }
	    return uh;
	}
    template<int N, typename Q, typename S>
	Eigen::VectorXd FESpace<N,Q,S>::M_loc_rhs(Geo::Triangle const & ie , std::function<std::array<double,2>(double,double)> const & fun)const
	{
	    double dx;
	    Geo::Point p;
	    std::array<double,2> F;
	    Eigen::VectorXd  out(2*M_nln);
	    out = Eigen::VectorXd::Zero(2*M_nln);
	    for(unsigned int i=0; i<Q::nqn2d; ++i)
	    {
		dx = std::abs(ie.detJ())*M_quad.iweight(i);
		p = ie.map(M_quad.inode(i));
		F = fun(p.x(),p.y());
		out.segment(0,M_nln) += F[0]*dx*M_B.row(i);
		out.segment(M_nln,M_nln) += F[1]*dx*M_B.row(i);
	    }
	    //std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" <<std::endl;
	    //std::cout << *out << std::endl;
	    //std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" <<std::endl;
	    return out;
	}


   }
#endif

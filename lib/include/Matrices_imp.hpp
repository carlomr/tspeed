/**
 * @file Matrices_imp.hpp
 * @brief Implementation of the matrices for the method - templated part
 *
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
namespace Tspeed{
    template <int N, typename Q, typename T>
	void Matrices::sigmaeps_edge(Geo::Triangle const & ie, FESpace_ptr<N,Q,T> Xh,std::array<Eigen::Matrix2d,6*Q::nqn1d*T::gdl> &s,std::array<Eigen::Matrix2d,6*Q::nqn1d*T::gdl> &e, double lambda, double mu  )
	{
	    Eigen::Vector2d thisGrad;
	    //auto kien = [](unsigned int k, unsigned int i, unsigned int e, unsigned int n){return k*T::gdl*6+i*6+e*2+n;};

	    for(auto iedg: ie.all_edges())
	    {
		for(unsigned int k=0; k<Q::nqn1d; ++k) 
		{
		    for(unsigned int i=0; i<Xh->nln(); ++i)
		    {
			//std::cout << Xh->grad(k,i) << std::endl;
			//std::cout << ie.invJac() << std::endl;
			thisGrad = Xh->g_edge(k,i,iedg.id()).transpose() * ie.invJac();
			e[kien(T::gdl,k,i,iedg.id(),0)] << thisGrad(0), .5*thisGrad(1), .5*thisGrad(1), 0;
			e[kien(T::gdl,k,i,iedg.id(),1)] << 0., .5*thisGrad(0), .5*thisGrad(0), thisGrad(1);
			for(unsigned int n : {0,1})
			{
			    s[kien(T::gdl,k,i,iedg.id(),n)] = CTensorProduct(e[kien(T::gdl,k,i,iedg.id(),n)], lambda, mu);
			}
		    }    
		}
	    }

	}
template <int N, typename Q, typename T>
	void Matrices::sigmaeps(Geo::Triangle const & ie, FESpace_ptr<N,Q,T> Xh,std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> &s,std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> &e, double lambda, double mu  )
	{
	    Eigen::Vector2d thisGrad;
	    //auto kin = [](unsigned int k, unsigned int i, unsigned int n){return k*T::gdl*2+i*2+n;};
	    for(unsigned int k=0; k<Q::nqn2d; ++k) 
	    {
		for(unsigned int i=0; i<Xh->nln(); ++i)
		{
		    //std::cout << Xh->grad(k,i) << std::endl;
		    //std::cout << ie.invJac() << std::endl;
		    thisGrad = Xh->grad(k,i).transpose() * ie.invJac();
		    e[kin(T::gdl,k,i,0)] << thisGrad(0), .5*thisGrad(1), .5*thisGrad(1), 0;
		    e[kin(T::gdl,k,i,1)] << 0., .5*thisGrad(0), .5*thisGrad(0), thisGrad(1);
		    for(unsigned int n : {0,1})
		    {
			s[kin(T::gdl,k,i,n)] << (lambda + 2*mu) * e[kin(T::gdl,k,i,n)](0,0) + lambda * e[kin(T::gdl,k,i,n)](1,1) ,
			    2*mu*e[kin(T::gdl,k,i,n)](0,1) , 2*mu*e[kin(T::gdl,k,i,n)](0,1) , (lambda + 2*mu) * e[kin(T::gdl,k,i,n)](1,1) + lambda * e[kin(T::gdl,k,i,n)](0,0);
		    }


		}    
	    }

	}
    template<int N, typename Q, typename T>
	void Matrices::stress_loc(Geo::Triangle const & ie, std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> &s,std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> &e, FESpace_ptr<N,Q,T> Xh)
    {
	//unsigned int startIndex=ie.id()*T::gdl;
	double dx;
	for(unsigned int k=0; k<Q::nqn2d; ++k)
	{
	    dx = std::abs(ie.detJ())*Xh->quad().iweight(k);
	    for(unsigned int j=0; j<T::gdl; ++j)
	    {
		for(unsigned int i=0; i< T::gdl; ++i)
		{
		    for(int m:{0,1})
			for(int n:{0,1})
			{
			    A.component(m,n).block(ie.id(), ie.id())(i,j) += mat_dot(s[kin(T::gdl,k,j,n)] , e[kin(T::gdl,k,i,m)])*dx;
			}
		}
	    }
	}

    }
    template<int N, typename Q, typename T> 
	void Matrices::interelement(FESpace_ptr<N,Q,T> Xh, Geo::Triangle const & ie,std::array<Eigen::Matrix2d,3*2*Q::nqn1d*T::gdl> &s_e)
	{
	    unsigned int kk;
	    double ds;
	    Eigen::Vector2d n_e;
	    double h_scaled;
	    int neig_ie;
	    int neigedge_ie;
	    for(auto iedg: ie.all_edges())
	    {
		neig_ie = ie.neigh(iedg.id());
		if(neig_ie>=0)
		{
		    h_scaled = 1./iedg.length();
		    n_e = iedg.normal();

		    neigedge_ie = ie.neighedges(iedg.id());
		    for(unsigned int k=0; k<Q::nqn1d; ++k)
		    {
			kk=Q::nqn1d - k - 1;
			ds = iedg.length()*Xh->quad().eweight(k);
			for(unsigned int i=0; i<T::gdl; ++i)
			{
			    for(unsigned int j=0; j<T::gdl; ++j)
			    {
				for(int m:{0,1})
				{
				    for(int n:{0,1})
				    {	
					//I.component(m,n).block(ie.id(), ie.id())(i,j) += 0.1;
					//I.component(m,n).block(ie.id(), neig_ie)(i,j) -= 0.1;
					I.component(m,n).block(ie.id(), ie.id())(i,j) -= 0.5*s_e[kien(T::gdl,k,i,iedg.id(),m)].row(n).dot(n_e)*Xh->b_edge(k,iedg.id())(j)*ds;
					I.component(m,n).block(ie.id(), neig_ie)(i,j) += 0.5*s_e[kien(T::gdl,k,i,iedg.id(),m)].row(n).dot(n_e)*Xh->b_edge(kk,neigedge_ie)(j)*ds;
				    }
				}
			    }
			}
		    }
		    //I.component(1,0).block(ie.id(), neig_ie) = I.component(0,1).block(ie.id(), neig_ie);
		}
		//std::cout << h_scaled << std::endl;
	    }
	    //I.component(1,0).block(ie.id(),ie.id()) = I.component(0,1).block(ie.id(), ie.id());
	    //std::cout<< S.component(0,0).block(0, 0) << std::endl;
	    //std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;
	}
template<int N, typename Q, typename T> 
	void Matrices::stability(FESpace_ptr<N,Q,T> Xh, Geo::Triangle const & ie,Parameters const & p)
	{
	    unsigned int kk;
	    double ds;
	    Eigen::Vector2d n_e;
	    double h_scaled;
	    int neig_ie;
	    int neigedge_ie;
	    double avg_l, avg_m;
	    for(auto iedg: ie.all_edges())
	    {
		neig_ie = ie.neigh(iedg.id());
		if(neig_ie>=0)
		{
		    h_scaled = 1./iedg.length();
		    n_e = iedg.normal();

		    neigedge_ie = ie.neighedges(iedg.id());
		    avg_l = p.avg_p("lambda",ie.id(), neig_ie );
		    avg_m = p.avg_p("mu",ie.id(), neig_ie );
		    for(unsigned int k=0; k<Q::nqn1d; ++k)
		    {
			kk=Q::nqn1d - k - 1;
			ds = iedg.length()*Xh->quad().eweight(k);
			S.component(0,0).block(ie.id(), ie.id()) += ((avg_l+2*avg_m)*n_e(0)*n_e(0)+2*avg_m*n_e(1)*n_e(1))*h_scaled*Xh->b_edge(k,iedg.id())*Xh->b_edge(k,iedg.id()).transpose()*ds;
			S.component(1,1).block(ie.id(), ie.id()) += ((avg_l+2*avg_m)*n_e(1)*n_e(1)+2*avg_m*n_e(0)*n_e(0))*h_scaled*Xh->b_edge(k,iedg.id())*Xh->b_edge(k,iedg.id()).transpose()*ds;
			S.component(0,1).block(ie.id(), ie.id()) +=  (avg_l*n_e(0)*n_e(1))*h_scaled*Xh->b_edge(k,iedg.id())*Xh->b_edge(k,iedg.id()).transpose()*ds;

			S.component(0,0).block(ie.id(), neig_ie) -= ((avg_l+2*avg_m)*n_e(0)*n_e(0)+2*avg_m*n_e(1)*n_e(1))*h_scaled*Xh->b_edge(k,iedg.id())*Xh->b_edge(kk,neigedge_ie).transpose()*ds;
			S.component(1,1).block(ie.id(), neig_ie) -= ((avg_l+2*avg_m)*n_e(1)*n_e(1)+2*avg_m*n_e(0)*n_e(0))*h_scaled*Xh->b_edge(k,iedg.id())*Xh->b_edge(kk,neigedge_ie).transpose()*ds;
			S.component(0,1).block(ie.id(), neig_ie) -=  (avg_l*n_e(0)*n_e(1))*h_scaled*Xh->b_edge(k,iedg.id())*Xh->b_edge(kk,neigedge_ie).transpose()*ds;
		    }
		    S.component(1,0).block(ie.id(), neig_ie) = S.component(0,1).block(ie.id(), neig_ie);
		}
		//std::cout << h_scaled << std::endl;
	    }
	    S.component(1,0).block(ie.id(),ie.id()) = S.component(0,1).block(ie.id(), ie.id());
	    //std::cout<< S.component(0,0).block(0, 0) << std::endl;
	    //std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;
	}
    template<int N, typename Q, typename T>
	void Matrices::mass(FESpace_ptr<N,Q,T> Xh, Geo::Triangle const & ie, double rho)
    {
	M.component(1).block(ie.id(),ie.id()) = rho*ie.detJ()*Xh->base_mass();
	invM.component(1).block(ie.id(),ie.id()) = 1./rho*1./ie.detJ()*Xh->base_invmass();
	M.component(0).block(ie.id(),ie.id()) = rho*ie.detJ()*Xh->base_mass();
	invM.component(0).block(ie.id(),ie.id()) = 1./rho*1./ie.detJ()*Xh->base_invmass();
    }

        template <int N, typename Q, typename T >
	Matrices::Matrices(FESpace_ptr<N,Q,T> Xh, Parameters const & P):A(Xh->mesh(), Xh->nln()),S(Xh->mesh(), Xh->nln()),I(Xh->mesh(), Xh->nln()),M(Xh->mesh(), Xh->nln()),invM(Xh->mesh(), Xh->nln())
	{
	    std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> sigma;
	    std::array<Eigen::Matrix2d,3*2*Q::nqn1d*T::gdl> sigma_edge;
	    std::array<Eigen::Matrix2d,2*Q::nqn2d*T::gdl> eps;
	    std::array<Eigen::Matrix2d,3*2*Q::nqn1d*T::gdl> eps_edge;
	    for(auto ie: Xh->mesh()->elements())
	    {
		sigmaeps(ie, Xh, sigma, eps, P.lambda(ie.id()), P.mu(ie.id()));
		sigmaeps_edge(ie, Xh, sigma_edge, eps_edge, P.lambda(ie.id()), P.mu(ie.id()));
		stress_loc(ie, sigma, eps, Xh);
		stability(Xh, ie, P);
		interelement(Xh, ie, sigma_edge);
		mass(Xh, ie, P.rho(ie.id()));
	    }

	    //std::cout << I.component(1,0).block(1,0) << std::endl;
	    //std::cout << "++++++++++++++++++++++++++++++" << std::endl;
	    //std::cout << I.component(0,1).block(0,1) << std::endl;
	    //std::cout << "++++++++++++++++++++++++++++++" << std::endl;
	    I.symmetrize();
	    //std::cout << I.component(0,1).block(1,0) << std::endl;
	}
}

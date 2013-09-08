/**
 * @file ShapeFunctions_imp.hpp
 * @brief Implementation of the shape functions
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
namespace Tspeed
{
    Eigen::ArrayXd jacobi_polynomial(int N, int alpha, int beta, Eigen::ArrayXd const & z);
    template<int N>
	Dubiner<N>::Dubiner()
	{
	   this->M_phi.reserve((N+1)*(N+2)/2);
	   this->M_grad.reserve((N+1)*(N+2)/2);
	   for(int j = 0; j<=N; ++j)
	   { 
		for(int i = 0; i<=N-j; ++i)
		{
		    double cij=std::sqrt((2.*i +1)*2.*(i+j+1)/pow(4.,i));
		    this->M_phi.push_back([=](Arr const & csi,Arr const &eta)->Arr{return cij*(pow(2.,i))*(pow((1-eta),i))*jacobi_polynomial(i,0,0,2*csi/(1-eta)-1)*jacobi_polynomial(j,2*i+1,0,2*eta-1);});
		    if (i==0 && j==0) 
			this->M_grad.push_back([](Arr const & csi, Arr const & eta)->ArrG{return ArrG::Zero(csi.size(),2);});
		    else if (i==0 && j!=0) 
			this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ArrG grad(csi.size(),2);
						grad.col(0)=Arr::Zero(csi.size(),1);
						grad.col(1)=Arr(cij*(j+2)*jacobi_polynomial(j-1,2,1,2*eta-1));
						return grad;});
		    else if (i!=0 & j==0) 
			this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ArrG grad(csi.size(),2);
						grad.col(0)=cij*pow(2.,i)*pow((1-eta),(i-1))*(i+1)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1);
						grad.col(1)=cij*pow(2,i)*(-i*pow((1-eta),(i-1))*jacobi_polynomial(i,0,0,2*csi/(1-eta)-1) +csi*pow((1-eta),(i-2))*(i+1)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1)); 
						return grad;});
		    else 
			this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ArrG grad(csi.size(),2);
						grad.col(0)=cij*pow(2,i)*pow((1-eta),(i-1))*(i+1)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1)*jacobi_polynomial(j,2*i+1,0,2*eta-1); 
						grad.col(1)=cij*pow(2,i)*(-i*pow((1-eta),(i-1))*jacobi_polynomial(i,0,0,2*csi/(1-eta)-1)*jacobi_polynomial(j,2*i+1,0,2*eta-1) +csi*pow((1-eta),(i-2))*(i+1)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1)*jacobi_polynomial(j,2*i+1,0,2*eta-1) +pow((1-eta),(i))*(2*i+j+2)*jacobi_polynomial(i,0,0,2*csi/(1-eta)-1)*jacobi_polynomial(j-1,2*i+2,1,2*eta-1)); 
						return grad;});
		}
	    }
	};
    template<int N>
	BoundaryAdapted<N>::BoundaryAdapted()
	{
	   this->M_phi.reserve((N+1)*(N+2)/2);
	   this->M_grad.reserve((N+1)*(N+2)/2);
	   this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return 1-csi-eta;});
	   this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ ArrG grad(csi.size(),2);
		   grad.col(0)=-Arr::Ones(csi.size());
		   grad.col(1)=-Arr::Ones(csi.size());
		   return grad;});

	   this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return csi ;});
	   this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ ArrG grad(csi.size(),2);
		   grad.col(0)= Arr::Ones(csi.size());
		   grad.col(1)= Arr::Zero(csi.size());
		   return grad;});
	   this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return eta;});
	   this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ ArrG grad(csi.size(),2);
		   grad.col(0)=Arr::Zero(csi.size());
		   grad.col(1)=Arr::Ones(csi.size());
		   return grad;});



	   for(int i=1; i<N;++i)
	   {
	       this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return (1-csi-eta)*csi*(1-eta).pow(i-1) * jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1);});
	       this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ArrG grad(csi.size(), 2);
	       grad.col(0)  = (1-2*csi-eta)*(1-eta).pow(i-1)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1) + (1-csi-eta)*csi*(1-eta).pow(i-2)*jacobi_polynomial(i-2,2,2,2*csi/(1-eta)-1)*(i+2);
	       grad.col(1)  = -(1-eta).pow(i-2)*csi*(i+csi-i*csi-i*eta)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1) + csi*csi*(1-csi-eta)*(1-eta).pow(i-3)*(i+2)*jacobi_polynomial(i-2,2,2,2*csi/(1-eta)-1);
	       return grad;});
	   }
	   for(int j=1; j<N; ++j)
	   {
	       this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return  (1-csi-eta)*eta*jacobi_polynomial(j-1,1,1,2*eta-1);});
	       this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ ArrG grad(csi.size(),2);
		       grad.col(0)=-eta*jacobi_polynomial(j-1,1,1,2*eta-1);;
		       grad.col(1)=(1-csi-2*eta)*jacobi_polynomial(j-1,1,1,2*eta-1) + (1-csi-eta)*eta*jacobi_polynomial(j-2,2,2,2*eta-1)*(j+2);;
		       return grad;});
	   }
	   for(int j=1; j<N; ++j)
	   {
	       this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return  csi*eta*jacobi_polynomial(j-1,1,1,2*eta-1);});
	       this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ ArrG grad(csi.size(),2);
		       grad.col(0)=eta*jacobi_polynomial(j-1,1,1,2*eta-1);
		       grad.col(1)=csi*jacobi_polynomial(j-1,1,1,2*eta-1) + eta* csi*(j+2)*jacobi_polynomial(j-2,2,2,2*eta-1);
		       return grad;});
	   }
	   for(int j=1; j<N; ++j)
	   {
	       for(int i=1; i<N; ++i)
	       {
		   if(i+j < N)
		   {
		       this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return  (1-csi-eta)*csi*eta*(1-eta).pow(i-1) * jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1) * jacobi_polynomial(j-1,2*i+1,1,2*eta-1);});
		       this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ ArrG grad(csi.size(),2);
			       grad.col(0)=eta*(1-eta).pow(i-2)*jacobi_polynomial(j-1,2*i+1,1,2*eta-1)*( (1-2*csi-eta)*(1-eta)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1)+ (1-csi-eta)*csi*(i+2)*jacobi_polynomial(i-2,2,2,2*csi/(1-eta)-1));
			       grad.col(1)=csi*(1-eta).pow(i-3)* ( ( (1-eta).pow(2)*(1-csi-2*eta)-(1-eta)*eta*(1-csi-eta)*(i-1) ) *jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1)*jacobi_polynomial(j-1,2*i+1,1,2*eta-1) + csi*(1-csi-eta)*eta*(i+2)*jacobi_polynomial(i-2,2,2,2*csi/(1-eta)-1)*jacobi_polynomial(j-1,2*i+1,1,2*eta-1) + (1-csi-eta)*eta*(1-eta).pow(2)*(j+2*i+2)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1)*jacobi_polynomial(j-2,2*i+2,2,2*eta-1) );
			       return grad;});

		   }
	       }
	   }
	}

}

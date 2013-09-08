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



	   for(int i=1; i<N-1;++i)
	   {
	       this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return (1-csi-eta)*csi*pow((1-eta),(i-1)) * jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1);});
	       this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ArrG grad(csi.size(), 2);
	       grad.col(0)  = (1-2*csi-eta)*pow((1-eta),(i-1))*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1) + (1-csi-eta)*csi*pow((1-eta),(i-2))*jacobi_polynomial(i-2,2,2,2*csi/(1-eta)-1)*(i+2);
	       grad.col(1)  = -pow((1-eta),(i-2))*csi*(i+csi-i*csi-i*eta)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1) + csi*csi*(1-csi-eta)*pow((1-eta),(i-3))*(i+2)*jacobi_polynomial(i-2,2,2,2*csi/(1-eta)-1);
	       return grad;});
	   }
	   for(int j=1; j<N-1; ++j)
	   {
	       this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return  (1-csi-eta)*eta*jacobi_polynomial(j-1,1,1,2*eta-1);});
	       this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ ArrG grad(csi.size(),2);
		       grad.col(0)=-eta*jacobi_polynomial(j-1,1,1,2*eta-1);;
		       grad.col(1)=(1-csi-2*eta)*jacobi_polynomial(j-1,1,1,2*eta-1) + (1-csi-eta)*eta*jacobi_polynomial(j-2,2,2,2*eta-1)*(j+2);;
		       return grad;});
	   }
	   for(int j=1; j<N-1; ++j)
	   {
	       this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return  csi*eta*jacobi_polynomial(j-1,1,1,2*eta-1);});
	       this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ ArrG grad(csi.size(),2);
		       grad.col(0)=eta*jacobi_polynomial(j-1,1,1,2*eta-1);
		       grad.col(1)=csi*jacobi_polynomial(j-1,1,1,2*eta-1) + eta* csi*(j+2)*jacobi_polynomial(j-2,2,2,2*eta-1);
		       return grad;});
	   }
	   for(int j=1; j<N-1; ++j)
	   {
	       for(int i=1; i<N-1; ++i)
	       {
		   if(i+j < N)
		   {
		       this->M_phi.push_back([=](Arr const & csi, Arr const & eta)->Arr{return  (1-csi-eta)*csi*eta*pow((1-eta),(i-1)) * jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1) * jacobi_polynomial(j-1,2*i+1,1,2*eta-1);});
		       this->M_grad.push_back([=](Arr const & csi, Arr const & eta)->ArrG{ ArrG grad(csi.size(),2);
			       grad.col(0)=eta*pow((1-eta),(i-2))*jacobi_polynomial(j-1,2*i+1,1,2*eta-1)*( (1-2*csi-eta)*(1-eta)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1)+ (1-csi-eta)*csi*(i+2)*jacobi_polynomial(i-2,2,2,2*csi/(1-eta)-1));
			       grad.col(1)=csi*pow((1-eta),(i-3))* ( ( pow((1-eta),2)*(1-csi-2*eta)-(1-eta)*eta*(1-csi-eta)*(i-1) ) *jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1)*jacobi_polynomial(j-1,2*i+1,1,2*eta-1) + csi*(1-csi-eta)*eta*(i+2)*jacobi_polynomial(i-2,2,2,2*csi/(1-eta)-1)*jacobi_polynomial(j-1,2*i+1,1,2*eta-1) + (1-csi-eta)*eta*pow((1-eta),2)*(j+2*i+2)*jacobi_polynomial(i-1,1,1,2*csi/(1-eta)-1)*jacobi_polynomial(j-2,2*i+2,2,2*eta-1) );
			       return grad;});

		   }
	       }
	   }
	}

}

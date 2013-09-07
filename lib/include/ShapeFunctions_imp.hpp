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
	}

}

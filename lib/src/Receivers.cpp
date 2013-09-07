#include"Receivers.hpp"

namespace Tspeed{

    void Receivers::add(double const & x, double const &y, unsigned int const & ir, unsigned int const & step)
    {
	//std::cout << M_nel << " " << step+1 << std::endl;
	M_val.resize(M_nel, step+1);
	M_val(ir, step) = std::array<double,2>{{x,y}};

    }
    void Receivers::write(std::string const & fn)const
    {
	std::ofstream outf;
	for(unsigned int i = 0; i<M_nel; ++i )
	{
	    outf.open(fn+"_rcv_"+std::to_string(i)+".out");
	    for(unsigned int j = 0; j<M_val.cols(); ++j )
	    {
		outf << M_val(i,j)[0] << " " << M_val(i,j)[1] << std::endl;
	    }
	    outf.close();
	}
    }

}

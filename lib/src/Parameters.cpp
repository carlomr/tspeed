#include"FESpace.hpp"
namespace Tspeed{
 void Parameters::setp(std::string const & p, const int lab, const double lambda)
    {
	unsigned int count = 0;
	for(auto ie: M_mesh->elements())
	{

	    if(ie.reg() == lab)
	   {
	       (*M_map[p])(ie.id()) = lambda;
	       ++count;
	   }
	}
	std::cout << "Param :: set "<<p << " = " << lambda << " on " << count << " elements." << std::endl;
    };
 double Parameters::avg_p(std::string const & p,int i, int j)const
 {
     double loc_p = (*M_map.at(p))(i);
     double out_p = (*M_map.at(p))(j);
     if(loc_p+out_p == 0)
	 return 0;
     else
	 return 2*loc_p*out_p/(loc_p+out_p);
 }

}

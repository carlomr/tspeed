#include"Force.hpp"
namespace Tspeed{
    Eigen::VectorXd PointWiseForce::eval(const double & t)const
    {
		
	std::array<double, 2> F = M_f(t);
	Vec rhs=Vec::Zero(M_totSize);
	rhs.segment(M_startIndex1,M_size) = F[0] * M_shape[0];
	rhs.segment(M_startIndex2,M_size) = F[1] * M_shape[0];
	return rhs;
    }

}

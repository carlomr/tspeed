#include"TSPEED.hpp"
#include<iostream>
//#include<memory>

using namespace Tspeed;

int main()
{
    const double dt = 2.5e-4;
    const double tmax = 1;
    double t=0;
    double error = 0;

    //Mesh Th(std::string("square_Nnull_4tria.msh"));
    Mesh Th(std::string("./Gmsh_meshes/Lamb_new_unstruct.msh"));
    Th.stats();

    FESpace<4> Xh(Th);

    Parameters p(Th);

    //p.setp("lambda", 6, 2);
    //p.setp("mu", 6, 1);
    //p.setp("rho", 6, 1);
    p.setp("lambda", 6, 6.8270e+09);
    p.setp("mu", 6,6.8265125e+09);
    p.setp("rho", 6, 2000 );

    //Receivers r(Xh, Geo::Point(0,0));
    //Receivers r2(Xh, std::string("Prova.rcv"));
    Receivers r2(Xh, std::string("lamb.rcv"));

    
    //std::shared_ptr<Force> f( new PointWiseForce([](const double & t){return std::array<double,2>{{(1-2*M_PI*M_PI*(100)*t*t) * exp(-M_PI*M_PI*(100)*t*t),(1-2*M_PI*M_PI*(100)*t*t) * exp(-M_PI*M_PI*(100)*t*t)}};}, Geo::Point(0,0), Xh));
    std::shared_ptr<Force> f( new PointWiseForce([](const double & t){return std::array<double,2>{{0,-(1-2*M_PI*M_PI*(100)*(t-0.1)*(t-0.1)) * exp(-M_PI*M_PI*(100)*(t-0.1)*(t-0.1))}};}, Geo::Point(1500,1950), Xh));
    
    LeapFrog TA(Xh, p, r2);

    TA.set_dt(dt);
    TA.set_tmax(tmax);
    TA.add_force(f);
    TA.set_penalty(1);
    TA.set_initial_u(Xh,[](double x, double y){return std::array<double,2>{{0,0}};});
    TA.set_initial_v(Xh,[](double x, double y){return std::array<double,2>{{0,0}};});
    //TA.set_initial_u(Xh,[](double x, double y){return std::array<double,2>{{sin(M_PI*x)* sin(M_PI*y), cos(M_PI*x)*cos(M_PI*y)}};});
    //TA.set_initial_v(Xh,[](double x, double y){return std::array<double,2>{{ -sin(1e-4)*sin(M_PI*x)* sin(M_PI*y), -sin(1e-4)*cos(M_PI*x)*cos(M_PI*y)}};});

    TA.first_step();
    TA.eval_receivers();

    int step = 0;
    while(TA.is_running())
    {
	t+=dt;
	TA.step(t);
	if(++step%2 == 0)
	    TA.eval_receivers();
	//error = Xh.L2error([t](double x, double y){return std::array<double,2>{{ cos(t) *  sin(M_PI*x)* sin(M_PI*y), cos(t)*cos(M_PI*x)*cos(M_PI*y) }};}, TA.get_uh() ) ;
	//std::cout << "Time: " << t<< ", error L2 = " << error << std::endl;
    }
    TA.write_receivers("lamb_f50_N8_newF_oldP");

}


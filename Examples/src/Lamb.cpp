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
    Mesh_ptr Th(new Mesh(std::string("./Meshes/Lamb_new_unstruct.msh")));
    Th->stats();

    FESpace_ptr<5> Xh(new FESpace<5>(Th));

    Parameters p(Th);

    p.setp("lambda", 6, 6.8270e+09);
    p.setp("mu", 6,6.8265125e+09);
    p.setp("rho", 6, 2000 );

    Receivers r2(Xh, std::string("Receivers/lamb.rcv"));

    
    std::shared_ptr<Force> f( new PointWiseForce([](const double & t){return std::array<double,2>{{0,-(1-2*M_PI*M_PI*(100)*(t-0.1)*(t-0.1)) * exp(-M_PI*M_PI*(100)*(t-0.1)*(t-0.1))}};}, Geo::Point(1500,1950), Xh));
    
    LeapFrog TA(Xh, p, r2);

    TA.set_dt(dt);
    TA.set_tmax(tmax);
    TA.add_force(f);
    TA.set_penalty(1);
    TA.set_initial_u(Xh,[](double x, double y){return std::array<double,2>{{0,0}};});
    TA.set_initial_v(Xh,[](double x, double y){return std::array<double,2>{{0,0}};});

    TA.first_step();
    //TA.eval_receivers();

    int step = 0;
    while(TA.is_running())
    {
	t+=dt;
	TA.step(t);
	//if(++step%2 == 0)
	    //TA.eval_receivers();
    }
    //TA.write_receivers("lamb_f50_N8_newF_oldP");

}


#include"TSPEED.hpp"
#include<iostream>
//#include<memory>

using namespace Tspeed;

int main()
{
    const double dt = 5e-4;
    const double tmax = 1;
    double t=0;
    double error = 0;
    const int N = 3;

    //Mesh Th(std::string("square_Nnull_4tria.msh"));
    //Mesh_ptr Th(new Mesh(std::string("./Meshes/Lamb_new_unstruct.msh")));
    Mesh_ptr Th(new Mesh(std::string("./Meshes/Lamb_fullyunstruct3.msh")));
    Th->stats();

    FESpace_ptr<N, Dunavant<2*N>> Xh(new FESpace<N,Dunavant<2*N>>(Th));
    //FESpace_ptr<N> Xh(new FESpace<N>(Th));

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
    TA.eval_receivers();

    int step = 0;
    while(TA.is_running())
    {
	t+=dt;
	TA.step(t);
	if(++step%2 == 0)
	    TA.eval_receivers();
    }
    TA.write_receivers("lamb");

}


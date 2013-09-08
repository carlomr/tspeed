#include"TSPEED.hpp"
#include<iostream>
#include<memory>

using namespace Tspeed;
void wedge_init_param(double l, double m, double rho, double cf, double csurf , double & k, double & q, double & s, double &beta);

int main()
{
    const double dt = 5e-4;
    const double tmax = 1.0;
    double t=0;
    double error = 0;
    double k,q,s,beta;
    const double x0 = 1000;
    const double ymax = 2000;

    Mesh_ptr Th(new Mesh(std::string("./Meshes/wedge.msh")));
    Th->stats();

    //FESpace_ptr<5> Xh(new FESpace<5>(Th));
    FESpace_ptr<3,Dunavant<6>> Xh(new FESpace<3,Dunavant<6>>(Th));

    Parameters p(Th);

    p.setp("lambda", 1, 6.25e+09);
    p.setp("lambda", 2, 1.76e+10);

    p.setp("mu", 1, 3.125e+09);
    p.setp("mu", 2, 8.8e+09);

    p.setp("rho", 1, 2000);
    p.setp("rho", 2, 2200);

    Receivers r2(Xh, std::string("Receivers/wedge.rcv"));

    double csurf =  1865.05;
    wedge_init_param(1.76e+10, 8.8e+09,10,2200, 1865.05, k, q, s, beta);

    
    std::shared_ptr<Force> f( new PointWiseForce([](const double & t){return std::array<double,2>{{0,0}};}, Geo::Point(1500,1950), Xh));
    
    LeapFrog TA(Xh, p, r2);

    TA.set_dt(dt);
    TA.set_tmax(tmax);
    TA.add_force(f);
    TA.set_penalty(1);

    TA.set_initial_u(Xh,[=](double x, double y){return std::array<double,2>{{k*(exp(-q*(ymax-y))-2*q*s/(s*s+k*k)*exp(-s*(ymax-y)))* (1-beta*pow((0-(x-x0)/csurf),2))*exp(-beta*pow((0-(x-x0)/csurf),2)) , q*(exp(-q*(ymax-y))-2*k*k/(s*s+k*k)*exp(-s*(ymax-y))) * (1-beta*pow((0-(x-x0)/csurf),2))*exp(-beta*pow((0-(x-x0)/csurf),2))}};});
    TA.set_initial_v(Xh,[=](double x, double y){return std::array<double,2>{{k*(exp(-q*(ymax-y))-2*q*s/(s*s+k*k)*exp(-s*(ymax-y))) * ( 2*(beta*pow((dt+(x-x0)/csurf),2)-2)*beta*(dt+(x-x0)/csurf)*exp(-beta*pow((dt+(x-x0)/csurf),2))), q*(exp(-q*(ymax-y))-2*k*k/(s*s+k*k)*exp(-s*(ymax-y))) * ( 2*(beta*pow((dt+(x-x0)/csurf),2)-2)*beta*(dt+(x-x0)/csurf)*exp(-beta*pow((dt+(x-x0)/csurf),2)))}};});
		

    //Xh->points_out("wedge_points");

    int step = 0;
    TA.first_step();
    //Xh->field_out("wedge_field", TA.u(), step);
    TA.eval_receivers();

    while(TA.is_running())
    {
	++step;
	t+=dt;
	TA.step(t);
	if(step%2 == 0)
	    TA.eval_receivers();
	//if(step%20==0)
	    //Xh->field_out("wedge_field", TA.u(), step);
    }
    TA.write_receivers("wedge_f10_back");

}

void wedge_init_param(double lambda, double mu, double rho, double cf, double csurf ,  double & k, double & q, double & s, double &beta)
{
    beta      = M_PI*M_PI*(cf*3)*(cf*3);
    double c1 = sqrt((lambda+2*mu)/rho);
    double c2 = sqrt(mu/rho);


    double m = mu/(lambda+2*mu);

    k = (cf*3)/csurf;

    q = k*sqrt(1-m*(csurf/c2)*(csurf/c2));
    s = k*sqrt(1-(csurf/c2)*(csurf/c2));
}

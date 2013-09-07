#include"Mesh.hpp"
namespace{
    //! An helper function
    void skipInput(std::istream & f, std::string s){
	std::string  currLine;
	do {std::getline(f, currLine);}
	while( (currLine.find(s) == std::string::npos) && ! f.eof());
    }
}
//
namespace Tspeed{

    void Mesh::M_neighbors()
    {
	M_bed_map.reserve(M_bed.size());
	unsigned int v0,v1,v2,vn1,vn0,vn2;
	for(auto ie = M_tri.begin(); ie!=M_tri.end()-1; ++ie)
	{
	    v0 = ie->pt(0).id();
	    v1 = ie->pt(1).id();
	    v2 = ie->pt(2).id();
	    for(auto ie2 = ie+1; ie2 != M_tri.end(); ++ie2)
	    {

		vn0 = ie2->pt(0).id();
		vn1 = ie2->pt(1).id();
		vn2 = ie2->pt(2).id();
		
		if (v0==vn1 && v1==vn0)
		{
		    ie->setNeigh(2,ie2->id());
		    ie2->setNeigh(2,ie->id());
		    ie->setNeighedges(2,2);
		    ie2->setNeighedges(2,2);
		}
		else if (v0==vn0 && v1==vn2) 
		{
		    ie->setNeigh(2,ie2->id());
		    ie2->setNeigh(1,ie->id());
		    ie->setNeighedges(2,1);
		    ie2->setNeighedges(1,2);
		}
		else if (v0==vn2 && v1==vn1)
		{
		    ie->setNeigh(2,ie2->id());
		    ie2->setNeigh(0,ie->id());
		    ie->setNeighedges(2,0);
		    ie2->setNeighedges(0,2);
		}
		else if (v1==vn1 && v2==vn0)
		{
		    ie->setNeigh(0,ie2->id());
		    ie2->setNeigh(2,ie->id());
		    ie->setNeighedges(0,2);
		    ie2->setNeighedges(2,0);
		}
		else if (v1==vn0 && v2==vn2) 
		{
		    ie->setNeigh(0,ie2->id()); 
		    ie2->setNeigh(1,ie->id()); 
		    ie->setNeighedges(0,1);
		    ie2->setNeighedges(1,0);
		}
		else if (v1==vn2 && v2==vn1)
		{
		    ie->setNeigh(0,ie2->id());
		    ie2->setNeigh(0,ie->id());
		    ie->setNeighedges(0,0);
		    ie2->setNeighedges(0,0);
		}
		else if (v2==vn1 && v0==vn0)
		{
		    ie->setNeigh(1,ie2->id());
		    ie2->setNeigh(2,ie->id());
		    ie->setNeighedges(1,2);
		    ie2->setNeighedges(2,1);
		}
		else if (v2==vn0 && v0==vn2) 
		{
		    ie->setNeigh(1,ie2->id());
		    ie2->setNeigh(1,ie->id());
		    ie->setNeighedges(1,1);
		    ie2->setNeighedges(1,1);
		}
		else if (v2==vn2 && v0==vn1)
		{
		    ie->setNeigh(1,ie2->id());
		    ie2->setNeigh(0,ie->id());
		    ie->setNeighedges(1,0);
		    ie2->setNeighedges(0,1);
		}
	    }
	    if(ie->neigh(0)>0 && ie->neigh(1) && ie->neigh(2)>0)
		ie->bcId() = Bc::Internal;
	    for(int i:{0,1,2})
	    {
		if(ie->edg(i).bcId()==Bc::Unassigned)
		    M_bed_map.emplace_back(ie->id(), i);
	    }

	}
	auto iee = M_tri.end()-1;
	for(int i:{0,1,2})
	{
	    if(iee->edg(i).bcId()==Bc::Unassigned)
		M_bed_map.emplace_back(iee->id(), i);
	}
    }
    void Mesh::stats()const
    {
	double max_l=0;
	double min_l=std::numeric_limits<double>::max();
	Eigen::Vector3d l_edge, alpha;
	double h,k,rho,av_anis(0), anis, max_alpha, min_alpha(M_PI);
	for(auto ie : M_tri )
	{
	    int m = 0;
	    for (auto iedg : ie.all_edges())
	    {
		if(iedg.length() > max_l)
		    max_l = iedg.length();
		else if(iedg.length()< min_l)
		    min_l = iedg.length();
		l_edge(m++) = iedg.length();
	    }
	    h = l_edge.maxCoeff();
	    //std::cout <<h << std::endl;
	    k = l_edge.sum()/2.;
	    rho = std::sqrt( (k*Eigen::Vector3d::Ones() - l_edge).prod() / h);
	    av_anis += h/rho;
	    if(h/rho>anis)
		anis = h/rho;
	    alpha = Eigen::Vector3d::Zero();
	    for(int i=1;i < ie.all_pts().size(); ++i)
	    {
		alpha(i) = std::abs(std::acos(dot(ie.pt(i)-ie.pt(i-1), ie.pt(i)- ie.pt((i+1)%3))/(ie.edg(i-1).length()*ie.edg((i+1)%3).length())));
	    }
	    alpha(0) = M_PI-alpha.sum();
	    //std::cout << alpha << std::endl;
	    max_alpha = std::max(max_alpha, alpha.maxCoeff());
	    min_alpha = std::min(min_alpha, alpha.minCoeff());
	}
	av_anis = av_anis/M_tri.size();
	std::cout << "Mesh :: file name                     : "<<M_fname << std::endl<<
	             "Mesh :: number of elements            : "<<M_tri.size()<<std::endl<<
		     "Mesh :: minimumum edge length (h_min) : "<<min_l        <<std::endl<<
		     "Mesh :: maximum edge length (h_max)   : "<<max_l        <<std::endl<<
		     "Mesh :: max anisotropy ( (h/2r)_max ) : "<<anis/2       <<std::endl<<
		     "Mesh :: avg anisotropy ( (h/2r)_avg ) : "<<av_anis/2    <<std::endl<<
		     "Mesh :: max angle                     : "<<max_alpha/M_PI <<"pi" <<std::endl<<
		     "Mesh :: min angle                     : "<<min_alpha/M_PI <<"pi" <<std::endl;
    }
	
    Mesh::Mesh(const std::string fileName)
    {
	M_fname = fileName;
	std::ifstream f;
	f.open(fileName.c_str());
	std::string currLine;
	if ( !f.is_open() ) {
	    std::cerr<<
		"Mesh file does not exist or is corrupted"<<
		std::endl;
	    exit(1);
	}
	//Read the number of points and elements
	//if(this->M_verbose)
	    //std::cout << "[Mesh::Mesh] Reading mesh data" << std::endl;
	// Scan file lines
	// Skip until data or end of file
	skipInput(f,std::string("$Nodes"));
	// If the string is not found, abort
	if(f.eof()||f.fail()){
	    std::cerr << "FILE ERROR! Cannot find Nodes" << std::endl;
	    exit(2);
	}
	// Read number of points
	int nP;
	f >> nP;
	//std::cout<<"Number of points="<<nP<<std::endl;
	M_points.reserve(nP);
	// Fill point data structures
	double x,y;
	int id;
	for(int i = 0; i < nP; ++i) {
	    f >>id>> x>>y;
	    f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	    //std::cout<<x<<" " <<y<<" "<< id<<" "<<pl.size()<<std::endl;
	    //std::cout<<P[0]<<" " <<P[1]<<" "<< P.bcId()<<" "<<pl.size()<<std::endl;
	    M_points.emplace_back(x,y);
	    M_points[i].id() = i; //my IDs start from 0 !!!!!
	}
	skipInput(f,std::string("$Elements"));
	int nE;
	f >> nE;
	M_tri.reserve(nE);
	M_bed.reserve(nP);
	int p1, p2, p3, el_type, num_tags, phys_ent, geom_ent, dummy;
	int idT(0), idE(0);
	for(int i=0; i<nE; ++i)
	{
	    f >> dummy >> el_type>>num_tags >> phys_ent >> geom_ent;
	    //std::cout << dummy << " " << el_type << " " << num_tags << std::endl;
	    for (int j=0; j<num_tags-2; ++j)
		f >> dummy;
	    if(el_type == 2) //triangle
	    {
		f >> p1 >> p2 >> p3;
		//std::cout << i << ": "<< p1 << " " << p2 << " " <<p3 << std::endl;
		M_tri.emplace_back(M_points[p1-1], M_points[p2-1], M_points[p3-1]);
		M_tri[idT].reg() = geom_ent;
		M_tri[idT].id() = idT;
		//M_tri[idT].printNeigh();
		++idT;
	    }
	    else if(el_type == 1) //boundary element
	    {
		f >> p1 >> p2;
		M_bed.emplace_back(M_points[p1-1], M_points[p2-1]);
		M_bed[idE].reg() = geom_ent;
		M_bed[idE].id()=idE;
		++idE;
	    }
	    f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	M_tri.resize(idT);
	M_bed.resize(idE);
	//std::cout<<"Number of elements="<<M_tri.size()<<std::endl;
	//std::cout<<"Number of bound. edges="<<M_bed.size()<<std::endl;
	//for(auto j: M_points) std::cout << j << std::endl;
	//std::cout << "_______________" << std::endl;
	//for(auto i: M_tri) std::cout << i << std::endl;
	f.close();
	M_BcMap = {{"Dirichlet", Bc::Dirichlet}, {"Neumann" ,  Bc::Neumann}, {"Unassigned", Bc::Unassigned}, {"Internal", Bc::Internal}};
	M_BcMap_inv ={{Bc::Dirichlet, "Dirichlet"}, {Bc::Neumann, "Neumann"}, {Bc::Unassigned, "Unassigned"}, {Bc::Internal, "Internal"}};
	this->M_neighbors();
    }
    
}


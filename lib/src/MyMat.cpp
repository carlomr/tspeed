/**
 * @file MyMat.cpp
 * @brief Implementation of the matrix classes
 * @author Carlo Marcati
 * @date 2013-09-08
 */
/* This program is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version. 
 *  
 *  This program is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details. 
 *  
 *  You should have received a copy of the GNU General Public License 
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include"MyMat.hpp"
namespace Tspeed
{

    BaseMat::BaseMat():M_nr(0),M_nc(0),M_nln(0){};

    void BaseMat::setblock(unsigned int i, unsigned int j, Eigen::MatrixXd const & M)
    {
	unsigned int r=M_r[i];
	for(unsigned short int k = 0; k<M_r[i+1]-r; ++k)
	    if(M_c[r+k] == j)
		M_m[r+k] = M;
    }
    Eigen::MatrixXd  & BaseMat::block(unsigned int i, unsigned int j)
    {
	unsigned int r=M_r[i];
	for(unsigned short int k = 0; k<M_r[i+1]-r; ++k)
	    if(M_c[r+k] == j)
		return M_m[r+k];
	std::cerr << "Block not found: " << i << "," << j << std::endl;
	exit(1);
    }
    Eigen::MatrixXd const & BaseMat::block(unsigned int i, unsigned int j)const
    {
	unsigned int r=M_r[i];
	for(unsigned short int k = 0; k<M_r[i+1]-r; ++k)
	    if(M_c[r+k] == j)
		return M_m[r+k];
	std::cerr << "Block not found: " << i << "," << j << std::endl;
	exit(1);
    }
    MyMat::MyMat(MyMat const & m )
    {
	M_nr = m.M_nr;
	M_nc = m.M_nc;
	M_nln = m.M_nln;
	M_r = m.M_r;
	M_c = m.M_c;
	M_m.resize(m.M_m.size());
	for(unsigned int i=0; i<M_m.size(); ++i)
	    M_m[i] = m.M_m[i];
    }
/*    MyMat::operator=(MyMat const & m )*/
    //{
	//M_nr = m.M_nr;
	//M_nc = m.M_nc;
	//M_nln = m.M_nln;
	//M_r = m.M_r;
	//M_c = m.M_c;
	//M_m.resize(m.M_m.size());
	//for(int i=0; i<M_m.size(); ++i)
	    //M_m[i] = m.M_m[i];
    /*}*/
    Eigen::VectorXd  operator*(MyMatMultiDimBlockDiag<MyMatBlockDiag> const & A, Eigen::VectorXd const &v)
    {
	//Eigen::VectorXd * out = new Eigen::VectorXd(v.size());
	Eigen::VectorXd out(v.size());
	A.vecMult(v,out);
	//out->segment(0,out->size()/2) = A.component(0,0) * v.segment(0,v.size()/2) + A.component(0,1) * v.segment(v.size()/2,v.size());
	//out->segment(out->size()/2, out->size()) = A.component(1,0) * v.segment(0,v.size()/2) + A.component(1,1) * v.segment(v.size()/2,v.size());
	return out;
    }
    Eigen::VectorXd  operator*(MyMatMultiDim<MyMat> & A, Eigen::VectorXd const &v)
    {
	Eigen::VectorXd  out(v.size());
	A.vecMult(v,out);
	//out->segment(0,out->size()/2) = A.component(0,0) * v.segment(0,v.size()/2) + A.component(0,1) * v.segment(v.size()/2,v.size());
	//out->segment(out->size()/2, out->size()) = A.component(1,0) * v.segment(0,v.size()/2) + A.component(1,1) * v.segment(v.size()/2,v.size());
	return out;
    }Eigen::VectorXd  operator*(MyMatMultiDim<MyMatBlockDiag> & A, Eigen::VectorXd const &v)
    {
	Eigen::VectorXd  out(v.size());
	A.vecMult(v,out);
	//out->segment(0,out->size()/2) = A.component(0,0) * v.segment(0,v.size()/2) + A.component(0,1) * v.segment(v.size()/2,v.size());
	//out->segment(out->size()/2, out->size()) = A.component(1,0) * v.segment(0,v.size()/2) + A.component(1,1) * v.segment(v.size()/2,v.size());
	return out;
    }
    MyMat MyMat::operator*(double const &c)const
    {
	MyMat  out (*this);
	for(unsigned int i=0; i<this->M_m.size(); ++i)
	{
	    out.M_m[i] = c * this->M_m[i];
	}
	return out;
    }
    MyMat operator*(double const & c, MyMat const & M)
    {
	return M*c;
    }
    MyMatBlockDiag::MyMatBlockDiag(Mesh_ptr  Th, unsigned int nln){
	this->M_nr=nln*Th->ne();
	this->M_nc=nln*Th->ne();
	this->M_nln=nln;
	this->M_m.reserve(Th->ne());	
	this->M_c.reserve(Th->ne());
	this->M_r.reserve(Th->ne()+1);
	unsigned int count=0;
	for(auto ie: Th->elements())
	{
	    this->M_r.push_back(count);
	    this->M_m.push_back(Eigen::MatrixXd::Zero(nln,nln));
	    this->M_c.push_back(ie.id());
	    count++;
	}
	this->M_r.push_back(count);
    }

    MyMat::MyMat(Mesh_ptr  Th, unsigned int nln){
	M_nr=nln*Th->ne();
	M_nc=nln*Th->ne();
	M_nln=nln;
	M_m.reserve(Th->ne()*4);	
	M_c.reserve(Th->ne()*4);
	M_r.reserve(Th->ne()+1);
	unsigned int count=0;
	for(auto ie: Th->elements())
	{
	    M_r.push_back(count);
	    M_m.push_back(Eigen::MatrixXd::Zero(nln,nln));
	    M_c.push_back(ie.id());
	    count++;
	    for(auto iedg: ie.all_edges())
	    {
		if(ie.neigh(iedg.id())>=0)
		{
		    M_m.push_back(Eigen::MatrixXd::Zero(nln,nln));
		    M_c.push_back(ie.neigh(iedg.id()));
		    //std::cout << "Building block " << ie.id() << "," << ie.neigh(iedg.id())  << " " << iedg.id()<< std::endl;
		    count++;
		}
	    }
	}
	M_r.push_back(count);
    }
    
    void MyMat::sumtranspose(MyMat const & ot)
    {
	for(unsigned int i = 0; i<M_r.size()-1; ++i)
	{
	    M_m[M_r[i]] += ot.M_m[M_r[i]].transpose();
	    for(unsigned int j = M_r[i]+1; j<M_r[i+1]; ++j)
	    {
		//unsigned int thisRow = i;
		//unsigned int thisCol = M_c[j];
		//std::cout << i << ", " << j << ", " << M_c[j] << ", " << M_r[i] << std::endl;
		M_m[j] += ot.block(M_c[j],i).transpose();
	    }
	}
    }void MyMat::symmetrize()
    {
	for(unsigned int i = 0; i<M_r.size()-1; ++i)
	{
	    M_m[M_r[i]] += M_m[M_r[i]].transpose().eval();
	    for(unsigned int j = M_r[i]+1; j<M_r[i+1]; ++j)
	    {
		assert(i!=j);
		if(M_c[j]>i)
		{
		    M_m[j] += this->block(M_c[j],i).transpose();
		}
		else
		{
		    M_m[j] = this->block(M_c[j],i).transpose();
		}
	    }
	}
    }
    /*void MyMatBlockDiag::vecMult(Eigen::VectorXd const & x, Eigen::VectorXd & out)const*/
    //{
	////std::cout << M_nr << std::endl;
	////std::cout << x.size() << std::endl;
	//assert(M_nr==x.size());
	//unsigned int r;
	//for(unsigned int i=0; i<M_r.size()-1; ++i)
	//{
	    //r = M_r[i];
	    //out.segment(i*M_nln,M_nln) = M_m[r]*x.segment(i*M_nln,M_nln);
	//}
    /*}*/
    void BaseMat::vecMult(Eigen::VectorXd const & x, Eigen::VectorXd & out)const
    {
	//std::cout << M_nr << std::endl;
	//std::cout << x.size() << std::endl;
	assert(M_nr==x.size());
	unsigned int r;
	for(unsigned int i=0; i<M_r.size()-1; ++i)
	{
	    r = M_r[i];

	    out.segment(i*M_nln,M_nln) = M_m[r]*x.segment(i*M_nln,M_nln);
	    for(unsigned int j=1; j<M_r[i+1]-r;++j)
	    {
		//std::cout << M_nr << " " <<  r << " " << j << " " << M_m[r+j].cols() << " " << x.segment(i,M_nln).rows() << std::endl;
		out.segment(i*M_nln,M_nln) += M_m[r+j]*x.segment(M_c[j+r]*M_nln,M_nln);

	    }
	}
    }
    MyMatMultiDim<MyMat>  operator+(MyMatMultiDim<MyMat> const &a , MyMatMultiDim<MyMat> const& b)
    {
	MyMatMultiDim<MyMat> out;// = new MyMatMultiDim<MyMat>();
	for(int i:{0,1,2,3})
	    out.M_mat[i] = std::move(a.M_mat[i] + b.M_mat[i]);
	return out;
    }
    MyMatMultiDim<MyMat>  operator+(MyMatMultiDim<MyMat> const &a, MyMatMultiDim<MyMatBlockDiag> const & b)
    {
	MyMatMultiDim<MyMat> out;// = new MyMatMultiDim<MyMat>();
	for(int i:{0,1,2,3})
	    out.M_mat[i] = std::move(a.M_mat[i] + b.M_mat[i]);
	return out;
    }
    /*MyMat & sum_part(MyMat const & a, MyMat const & b)*/
    //{
	//MyMat *out = new MyMat();
	//out->M_m = a.M_m;	
	//out->M_c = a.M_c;
	//out->M_r = a.M_r;
	//out->M_nr = a.M_nr;
	//out->M_nc = a.M_nc;
	//out->M_nln = a.M_nln;
	//for(unsigned int i=0; i<b.M_m.size(); ++i)
	//{
		//out->M_m[a.M_r[i]] += b.M_m[i];
	//}
	//return *out;
    //}MyMat & sum(MyMat const & a, MyMat const & b)
    //{
	//MyMat *out = new MyMat();
	//M_m.reserve(a.M_m.size());	
	//M_c = a.M_c;
	//M_r = a.M_r;
	//M_nr = a.M_nr;
	//M_nc = a.M_nc;
	//M_nln = a.M_nln;
	//for(unsigned int i=0; i<a.M_m.size(); ++i)
	//{
	    //out->M_m.emplace_back(a.M_m[i]+b.M_m[i]); 
	//}
	//return *out;
    //}

    MyMat MyMat::operator+=(const MyMatBlockDiag & a)
    {
	for(unsigned int i=0; i<a.size(); ++i)
	{
	    M_m[M_r[i]] += a.elem(i);
	}
	return *this;
    }
    MyMat operator+(MyMatBlockDiag const &b , MyMat a)
    {
	return a+=b;
    }
    MyMat operator+(MyMat a, MyMatBlockDiag const & b)
    {
	return a+=b;
    }

    MyMat MyMat::operator+=(MyMat const & a)
    {
	//M_m.reserve(a.M_m.size());	
	assert(M_m.size() == a.M_m.size());
	for(unsigned int i=0; i<a.M_m.size(); ++i)
	{
	    M_m[i] += a.M_m[i]; 
	}
	return *this;
    }
    MyMat operator+(MyMat a, MyMat const &b)
    {
	return a+=b;
    }
    //void MyMat::sum(MyMatBlockDiag const & m, MyMat & out)const*/
    //{
	//out.M_m.reserve(this->M_m.size());	
	//out.M_c = this->M_c;
	//out.M_r = this->M_r;
	//for(unsigned int i=0; i<this->M_m.size(); ++i)
	//{
	    //if(M_r[i] == M_c[i])
		//out.M_m.push_back(this->M_m[i]+m.M_m[M_r[i]]); 
	    //else
		//out.M_m.push_back(this->M_m[i]);
	//}
    //}
    //void MyMat::sum(MyMat const & m, MyMat & out)const
    //{
	//out.M_m.reserve(this->M_m.size());	
	//out.M_c = this->M_c;
	//out.M_r = this->M_r;
	//for(unsigned int i=0; i<this->M_m.size(); ++i)
	//{
	   //out.M_m.push_back(this->M_m[i]+m.M_m[i]); 
	//}
    /*}*/
    Eigen::VectorXd operator*(MyMatBlockDiag const & A,Eigen::VectorXd const & x)
    {
	Eigen::VectorXd out(x.size());
	A.vecMult(x,out);
	return out;
    }
    Eigen::VectorXd operator*(MyMat const & A,Eigen::VectorXd const & x)
    {
	Eigen::VectorXd out(x.size());
	A.vecMult(x,out);
	return out;
    }
    
	
}

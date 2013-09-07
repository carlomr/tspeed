#ifndef __MYMAT_HPP__
#define __MYMAT_HPP__ 1
#include<Eigen/Dense>
#include<vector>
#include"Mesh.hpp"
#include<fstream>
namespace Tspeed{
    class BaseMat
    {
	public:
	    BaseMat();
	    BaseMat(Mesh_ptr , unsigned int nln) ;
	    Eigen::MatrixXd const & block(unsigned int i, unsigned int j)const;
	    Eigen::MatrixXd & block(unsigned int i, unsigned int j);
	    void setblock(unsigned int i, unsigned int j, Eigen::MatrixXd const & M);
	    BaseMat(BaseMat const&) ;
	    virtual ~BaseMat() = default;

	    unsigned int nr()const{return M_nr;};
	    std::vector<unsigned int> const & rowInd()const{return M_r;};
	    std::vector<unsigned int> const & colInd()const{return M_c;};
	    void set_rowInd(std::vector<unsigned int>const & v){M_r=v;};
	    void set_colInd(std::vector<unsigned int>const & v){M_c=v;};
	    std::vector<Eigen::MatrixXd> const & elem()const{return M_m;};
	    std::vector<Eigen::MatrixXd> & elem(){return M_m;};
	    Eigen::MatrixXd const & elem(int i)const{return M_m[i];};
	    unsigned int size()const{return M_m.size();};
	    void vecMult(Eigen::VectorXd const &, Eigen::VectorXd &)const;
	protected:
	    unsigned int M_nr;
	    unsigned int M_nc;
	    unsigned int M_nln;
	    std::vector<Eigen::MatrixXd> M_m;
	    std::vector<unsigned int> M_r;
	    std::vector<unsigned int> M_c;

    };
    class MyMatBlockDiag : public BaseMat 
    {
	public:
	    MyMatBlockDiag():BaseMat(){};
	    MyMatBlockDiag(Mesh_ptr  , unsigned int nln);
	    MyMatBlockDiag(MyMatBlockDiag&&) = default;
	    MyMatBlockDiag(MyMatBlockDiag const&) = default;
	    MyMatBlockDiag & operator=(MyMatBlockDiag&&) = default;
	    MyMatBlockDiag & operator=(MyMatBlockDiag const&) = default;
	    virtual ~MyMatBlockDiag()noexcept(true) = default;
    };

    class MyMat : public BaseMat
    {
	public:

	    MyMat():BaseMat(){};
	    MyMat(Mesh_ptr , unsigned int nln);
	    MyMat(MyMat&&) = default;
	    MyMat(MyMat const&) ;
	    MyMat&operator=(MyMat&&) = default;
	    MyMat&operator=(MyMat const&) = default;
	    virtual ~MyMat()noexcept(true) = default;

	    void symmetrize();
	    void sumtranspose(MyMat const &);
	    MyMat operator+=(MyMat const & );
	    MyMat operator+=(MyMatBlockDiag const & );
	    MyMat operator*(double const & )const;
    };
    MyMat operator*(double const & c, MyMat const & M);
    Eigen::VectorXd operator*(MyMat const &,Eigen::VectorXd const &);
    Eigen::VectorXd operator*(MyMatBlockDiag const &,Eigen::VectorXd const &);
    MyMat operator+(MyMat a, MyMat const &b);
    MyMat operator+(MyMat a, MyMatBlockDiag const &b);
    template<typename T>
	class MyMatMultiDim
    {
	public:
	    MyMatMultiDim() = default;
	    virtual ~MyMatMultiDim() = default;
	    MyMatMultiDim(Mesh_ptr , unsigned int nln);
	    T & component(int i, int j){return M_mat[i*2+j];};
	    T const & component(int i, int j)const{return M_mat[i*2+j];};
	    void symmetrize();
	    void vecMult(Eigen::VectorXd const &, Eigen::VectorXd &)const;
	    MyMatMultiDim(MyMatMultiDim&a){for (int i:{0,1,2,3}){this->M_mat[i] = a.M_mat[i];};};
	    MyMatMultiDim&operator=(MyMatMultiDim&&) = default;
	    Eigen::VectorXd operator*(Eigen::VectorXd const &v)const;
	    //template<>
	    friend MyMatMultiDim<MyMat>  operator+(MyMatMultiDim<MyMat> const &a , MyMatMultiDim<MyMat> const& b);
	    friend MyMatMultiDim<MyMat>  operator+(MyMatMultiDim<MyMat> const &a , MyMatMultiDim<MyMatBlockDiag> const& b);
	    friend MyMatMultiDim<T> operator*(double const & x, MyMatMultiDim<T> const & A){MyMatMultiDim out; for(int i:{0,1,2,3}) out.M_mat[i] = x*A.M_mat[i]; return out;};
	private:
	    std::array<T,4> M_mat;
    };
    template<typename T>
	class MyMatMultiDimBlockDiag
    {
	public:
	    MyMatMultiDimBlockDiag() = default;
	    virtual ~MyMatMultiDimBlockDiag() = default;
	    MyMatMultiDimBlockDiag(Mesh_ptr , unsigned int nln);
	    T & component(int i){return M_mat[i];};
	    T const & component(int i)const{return M_mat[i];};
	    void vecMult(Eigen::VectorXd const &, Eigen::VectorXd &)const;
	    MyMatMultiDimBlockDiag(MyMatMultiDimBlockDiag&&) = default;
	    //MyMatMultiDimBlockDiag(MyMatMultiDimBlockDiag&) = default;
	    MyMatMultiDimBlockDiag&operator=(MyMatMultiDimBlockDiag&&) = default;
	    unsigned int nr()const{return M_mat[0].nr()+M_mat[1].nr();};
	    //friend MyMatMultiDimBlockDiag<MyMat> & operator+(MyMatMultiDimBlockDiag<MyMat> const &a , MyMatMultiDimBlockDiag<MyMatBlockDiag> const& b);
	    Eigen::VectorXd operator*(Eigen::VectorXd const &v)const;
	    friend MyMatMultiDimBlockDiag const & operator*(double const & x, MyMatMultiDimBlockDiag const & A){MyMatMultiDimBlockDiag* out = new MyMatMultiDimBlockDiag(); for(int i:{0,1,2,3}) out->M_mat[i] = x*A.M_mat[i]; return *out;};
	private:
	    std::array<T,2> M_mat;
    };
    template<typename T>
	MyMatMultiDimBlockDiag<T>::MyMatMultiDimBlockDiag(Mesh_ptr  m, unsigned int nln)
    {
	for(int i:{ 0,1 })
	{
	    M_mat[i] = T(m, nln);
	}
    }

    template<typename T>
    MyMatMultiDim<T>::MyMatMultiDim(Mesh_ptr  m, unsigned int nln)
    {
	for(int i:{ 0,1,2,3 })
	{
	    M_mat[i] = T(m, nln);
	}
    }
    template <typename T>
    void MyMatMultiDimBlockDiag<T>::vecMult(Eigen::VectorXd const & x, Eigen::VectorXd & out)const
    {
	assert(x.size()==M_mat[0].nr()*2);
	out.segment(0,M_mat[0].nr()) = M_mat[0] * x.segment(0,M_mat[0].nr());
	out.segment(M_mat[0].nr(),M_mat[0].nr()) = M_mat[1]* x.segment(M_mat[0].nr(),M_mat[0].nr());
    }
template <typename T>
    void MyMatMultiDim<T>::vecMult(Eigen::VectorXd const & x, Eigen::VectorXd & out)const
    {
	assert(x.size()==M_mat[0].nr()*2);
	int hS = M_mat[0].nr();
	out.segment(0,hS) = M_mat[0] * x.segment(0,hS)+M_mat[1]* x.segment(hS,hS);
	out.segment(hS,hS) = M_mat[2] * x.segment(0,hS)+M_mat[3]* x.segment(hS,hS);

    }
    template<typename T>
    void MyMatMultiDim<T>::symmetrize()
    {
	M_mat[0].symmetrize();
	M_mat[3].symmetrize();
	T temp = M_mat[1];
	M_mat[1].sumtranspose(M_mat[2]);
	M_mat[2].sumtranspose(temp); //TODO:cosi fa schifo
    }
    template<typename T>
    Eigen::VectorXd MyMatMultiDim<T>::operator*( Eigen::VectorXd const &v)const
    {
	Eigen::VectorXd  out(v.size());
	this->vecMult(v,out);
	return out;
    }
    template<typename T>
    Eigen::VectorXd MyMatMultiDimBlockDiag<T>::operator*(Eigen::VectorXd const &v)const
    {
	Eigen::VectorXd  out(v.size());
	this->vecMult(v,out);
	return out;
    }
    
}
#endif

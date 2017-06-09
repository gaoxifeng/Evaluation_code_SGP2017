#ifndef DHM_UTIL_H
#define DHM_UTIL_H

#include "DHMTraits.h"
#include "DHMMesh.h"
#include <boost/functional/hash.hpp>
#include <boost/property_tree/ptree.hpp>

PRJ_BEGIN

template <typename T>
void putNoOverwrite(boost::property_tree::ptree& pt,const std::string& path,T val)
{
  if(pt.find(path) == pt.not_found())
    pt.put<T>(path,val);
}
//sort a short list
template <typename T>
static void sort2(T& A,T& B)
{
  if(A>B)std::swap(A,B);
}
template <typename T>
static void sort4(T& A,T& B,T& C,T& D)
{
  //round 1
  if(A>B)std::swap(A,B);
  if(B>C)std::swap(B,C);
  if(C>D)std::swap(C,D);
  //round 2
  if(A>B)std::swap(A,B);
  if(B>C)std::swap(B,C);
  //round 3
  if(A>B)std::swap(A,B);
  //ASSERT(A<=B && B<=C && C<=D)
}
template <typename T>
static void sort4Map(T& A,T& B,T& C,T& D,T& VA,T& VB,T& VC,T& VD)
{
  //round 1
  if(A>B) {
    std::swap(A,B);
    std::swap(VA,VB);
  }
  if(B>C) {
    std::swap(B,C);
    std::swap(VB,VC);
  }
  if(C>D) {
    std::swap(C,D);
    std::swap(VC,VD);
  }
  //round 2
  if(A>B) {
    std::swap(A,B);
    std::swap(VA,VB);
  }
  if(B>C) {
    std::swap(B,C);
    std::swap(VB,VC);
  }
  //round 3
  if(A>B) {
    std::swap(A,B);
    std::swap(VA,VB);
  }
  //ASSERT(A<=B && B<=C && C<=D)
}
//make simple mesh
class DHMMeshMaker : public DHMTraits
{
public:
  static void makeCube(POSES& poses,CELLIDS& cells,const BBox<scalar>& bb,const Vec3i& nrCell);
  static void makeCube(const std::string& path,const BBox<scalar>& bb,const Vec3i& nrCell);
  static void makeSphere(const std::string& path,const scalar rad,sizeType nrPad,const Vec3i& nrCell);
  static void makeDeformedCube(const std::string& path,const Vec3i& nrCell,scalar rot);
};
//mapping function and deformation gardient kernel
struct DHMMapping {
  template <typename MAT,typename VEC>
  static void MStencil(MAT& s,const VEC& coord) {
    stencil3D(s.data(),coord[0]*0.5f+0.5f,coord[1]*0.5f+0.5f,coord[2]*0.5f+0.5f);
  }
  static Vec3 M(const DHMCell& cell,const Vec3& coord) {
    return
      interp3D(cell._verts[0]->_pos,cell._verts[1]->_pos,
               cell._verts[2]->_pos,cell._verts[3]->_pos,
               cell._verts[4]->_pos,cell._verts[5]->_pos,
               cell._verts[6]->_pos,cell._verts[7]->_pos,
               coord[0]*0.5f+0.5f,coord[1]*0.5f+0.5f,coord[2]*0.5f+0.5f);
  }
  template <typename MAT,typename VEC>
  static void JStencil(MAT& s,const VEC& coord) {
    Eigen::Matrix<typename MAT::Scalar,4,1> coefs;
#define ADD_STENCIL(COL,D0,D1,D2)	\
s(COL,D0)=coefs[0]*0.5f;	\
s(COL,0)=-coefs[0]*0.5f;	\
s(COL,D0+D1)=coefs[1]*0.5f;	\
s(COL,0+D1)=-coefs[1]*0.5f;	\
s(COL,D0+D2)=coefs[2]*0.5f;	\
s(COL,0+D2)=-coefs[2]*0.5f;	\
s(COL,D0+D1+D2)=coefs[3]*0.5f;	\
s(COL,0+D1+D2)=-coefs[3]*0.5f;
    stencil2D(coefs.data(),coord[1]*0.5f+0.5f,coord[2]*0.5f+0.5f);
    ADD_STENCIL(0,1,2,4)

    stencil2D(coefs.data(),coord[0]*0.5f+0.5f,coord[2]*0.5f+0.5f);
    ADD_STENCIL(1,2,1,4)

    stencil2D(coefs.data(),coord[0]*0.5f+0.5f,coord[1]*0.5f+0.5f);
    ADD_STENCIL(2,4,1,2)
#undef ADD_STENCIL
  }
  template <typename MAT,typename MAT3>
  static void J(const MAT& s,const DHMCell& cell,MAT3& ret) {
    typedef typename MAT3::Scalar SCALAR;
    typedef Eigen::Matrix<SCALAR,3,1> VEC3;
    VEC3 p;
    ret.setZero();
    for(char v=0; v<8; v++) {
      p=cell._verts[v]->_pos.cast<SCALAR>();
      ret.col(0)+=s(0,v)*p;
      ret.col(1)+=s(1,v)*p;
      ret.col(2)+=s(2,v)*p;
    }
  }
  template <typename MAT,typename MATP,typename MAT3>
  static void JP(const MAT& s,const MATP& pos,MAT3& ret) {
    typedef typename MAT3::Scalar SCALAR;
    typedef Eigen::Matrix<SCALAR,3,1> VEC3;
    VEC3 p;
    ret.setZero();
    for(char v=0; v<8; v++) {
      p=pos.col(v).template cast<SCALAR>();
      ret.col(0)+=s(0,v)*p;
      ret.col(1)+=s(1,v)*p;
      ret.col(2)+=s(2,v)*p;
    }
  }
  template <typename JS,typename INVF,typename DFDX>
  static void calcDFDX(const JS& Js,const INVF& invF,DFDX& dFdX) {
    dFdX.setZero();
    for(char r=0; r<3; r++)
      for(char c=0; c<3; c++)
        for(char v=0; v<8; v++)
          dFdX(r+c*3,v*3+r)=
            Js.col(v).dot(invF.col(c));
  }
  static void debugMapping();
};
//matrix
template <typename T>
static void addI3x3(vector<Eigen::Triplet<T,sizeType> >& H,sizeType r,sizeType c,T coef)
{
  H.push_back(Eigen::Triplet<T,sizeType>(r+0,c+0,coef));
  H.push_back(Eigen::Triplet<T,sizeType>(r+1,c+1,coef));
  H.push_back(Eigen::Triplet<T,sizeType>(r+2,c+2,coef));
}
template <typename T,typename MT>
static void addI(vector<Eigen::Triplet<T,sizeType> >& H,sizeType r,sizeType c,const MT& coef)
{
  sizeType nr=coef.size();
  for(sizeType i=0; i<nr; i++)
    if(coef[i] != 0.0f)
      H.push_back(Eigen::Triplet<T,sizeType>(r+i,c+i,(T)coef[i]));
}
template <typename T,typename MT>
static void addBlock(vector<Eigen::Triplet<T,sizeType> >& H,sizeType r,sizeType c,const MT& coef)
{
  sizeType nrR=coef.rows();
  sizeType nrC=coef.cols();
  for(sizeType i=0; i<nrR; i++)
    for(sizeType j=0; j<nrC; j++)
      if(coef(i,j) != 0.0f)
        H.push_back(Eigen::Triplet<T,sizeType>(r+i,c+j,coef(i,j)));
}
template <typename T,typename MT>
static void addBlockPreAlloc(vector<Eigen::Triplet<T,sizeType> >& H,sizeType& off,sizeType r,sizeType c,const MT& coef)
{
  sizeType nrR=coef.rows();
  sizeType nrC=coef.cols();
  for(sizeType i=0; i<nrR; i++)
    for(sizeType j=0; j<nrC; j++)
      if(coef(i,j) != 0.0f)
        H[off++]=Eigen::Triplet<T,sizeType>(r+i,c+j,coef(i,j));
}
template <typename T,typename T2>
static void addSparseBlock(vector<Eigen::Triplet<T,sizeType> >& H,sizeType r,sizeType c,const Eigen::SparseMatrix<T2,0,sizeType>& mat)
{
  for(sizeType k=0; k<mat.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T2,0,sizeType>::InnerIterator it(mat,k); it; ++it)
      H.push_back(Eigen::Triplet<T,sizeType>(it.row()+r,it.col()+c,(T)it.value()));
}
//apply constraint
template <typename T>
static void applyConstraint(vector<Eigen::Triplet<T,sizeType> >& H,const boost::unordered_set<sizeType>& cs)
{
  for(sizeType i=0; i<(sizeType)H.size();)
    if(cs.find(H[i].row()) != cs.end() || cs.find(H[i].col()) != cs.end()) {
      H[i]=H.back();
      H.pop_back();
    } else i++;
}
template <typename VEC>
static void applyConstraint(VEC& H,const boost::unordered_set<sizeType>& cs)
{
  for(boost::unordered_set<sizeType>::const_iterator
      beg=cs.begin(),end=cs.end(); beg!=end; beg++)
    H[*beg]=0.0f;
}
//uniform validity check
struct DHMUniformValidityBound {
  static Eigen::Matrix<scalarD,27,1> bernstein(const Vec3d& pos);
  static Eigen::Matrix<scalarD,27,1> calcValidBound(const DHMCell& cell);
  static bool isValid(const DHMCell& cell);
  static void debugCorrect();
  static void debugTightness();
};

PRJ_END

#endif

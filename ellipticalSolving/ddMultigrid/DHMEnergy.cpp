#include "DHMEnergy.h"
#include "DHMUtil.h"

USE_PRJ_NAMESPACE

void DHMEnergyPool::debugEnergy()
{
  //build random cell
  DHMCell cell(0);
  for(char v=0; v<8; v++) {
    cell._verts[v].reset(new DHMVertex(v));
    cell._verts[v]->_pos.setRandom();
  }

  //debug energy
  {
    DHMARAPEnergy e(cell,1.0f);
    e.findDegree<5>();
    e.debugHessGrad();
  }
  {
    DHMConformalEnergy e(cell,1.0f);
    e.findDegree<5>();
    e.debugHessGrad();
  }
  {
    DHMElasticEnergy e(cell,1.0f,1.0f);
    e.findDegree<5>();
    e.debugHessGrad();
  }
  {
    DHMConsistencyEnergy e(cell,1.0f);
    e.findDegree<5>();
    e.debugHessGrad();
  }
  {
    DHMCell cellF(1);
    for(char v=0; v<8; v++) {
      cellF._verts[v].reset(new DHMVertex(v));
      cellF._verts[v]->_pos.setRandom();
    }
    Mat3d DFDC=Mat3d::Zero();
    DFDC(1,0)=1.0f;
    DFDC(0,1)=1.0f;
    DFDC(2,2)=1.0f;
    DHMGradConsistencyEnergy e(cellF,cell,DFDC,1.0f);
    e.findDegree<5>();
    e.debugHessGrad();
  }
}
//add energy
void DHMEnergyPool::addVDEnergy(boost::shared_ptr<DHMEnergy<24> > VEnergy)
{
  _VEnergys.push_back(VEnergy);
}
void DHMEnergyPool::clearEnergy()
{
  _VEnergys.clear();
}
//build
void DHMEnergyPool::buildHessian(TRIPS& trips)
{
  //initialize parallel adder
  sizeType nrVE=(sizeType)_VEnergys.size();
  DHMEnergy<24>::HESS hess;
  trips.assign(nrVE*576,Eigen::Triplet<scalarD,sizeType>(0,0,0.0f));

  //add hessian parallel
  OMP_PARALLEL_FOR_I(OMP_PRI(hess))
  for(sizeType i=0; i<nrVE; i++) {
    const DHMCell& cell=_VEnergys[i]->getCell();
    hess=_VEnergys[i]->hessInt();
    sizeType off=i*576;
    for(char vi=0; vi<8; vi++)
      for(char vj=0; vj<8; vj++) {
        sizeType idI=cell._verts[vi]->_index;
        sizeType idJ=cell._verts[vj]->_index;
        if(_fixed.find(idI) == _fixed.end() && _fixed.find(idJ) == _fixed.end())
          addBlockPreAlloc(trips,off,idI*3,idJ*3,hess.block<3,3>(vi*3,vj*3));
      }
  }

  //add to trips
  sizeType j=0;
  for(sizeType i=0; i<nrVE*576; i++)
    if(trips[i].value() != 0.0f)
      trips[j++]=trips[i];
  trips.resize(j);
}
void DHMEnergyPool::buildHessian(SMAT& hess)
{
  TRIPS trips;
  buildHessian(trips);
  for(boost::unordered_set<sizeType>::const_iterator
      beg=_fixed.begin(),end=_fixed.end(); beg!=end; beg++)
    addI<scalarD,Vec3d>(trips,3*(*beg),3*(*beg),Vec3d::Ones());
  hess.buildFromTripletsDepulicate(trips,0.0f);
}
void DHMEnergyPool::buildDeriv(COL& deriv)
{
  //initialize parallel adder
  sizeType nrVE=(sizeType)_VEnergys.size();
  DHMEnergy<24>::GRAD grad;

  //add deriv parallel
  OMP_PARALLEL_FOR_I(OMP_PRI(grad))
  for(sizeType i=0; i<nrVE; i++) {
    const DHMCell& cell=_VEnergys[i]->getCell();
    grad=_VEnergys[i]->gradInt();
    for(char vi=0; vi<8; vi++) {
      sizeType idI=cell._verts[vi]->_index;
      if(_fixed.find(idI) == _fixed.end())
        for(char d=0; d<3; d++) {
          scalarD& coeff=deriv[idI*3+d];
          OMP_ATOMIC_
          coeff+=grad[vi*3+d];
        }
    }
  }
}
void DHMEnergyPool::addDiagonal(SMAT& hess,scalar eps)
{
  sizeType nr=std::min(hess.rows(),hess.cols());
  for(sizeType i=0; i<nr; i++)
    hess.addToElement(i,i,eps);
}
void DHMEnergyPool::ensureDD(SMAT& hess)
{
  for(sizeType i=0; i<hess.rows(); i++) {
    scalarD DVal=0.0f;
    for(ConstSMIterator<scalarD> beg=hess.begin(i),end=hess.end(i); beg!=end; ++beg)
      DVal+=fabs(*beg);
    DVal-=hess(i,i);
    if(DVal > 0.0f)
      hess.addToElement(i,i,DVal);
  }
}
//arap energy
static Mat3d findR(const Eigen::JacobiSVD<Mat3d>& svd)
{
  for(char d=0; d<3; d++)
    ASSERT(svd.singularValues()[d] >= 0.0f)
    Mat3d U=svd.matrixU();
  Mat3d V=svd.matrixV();
  if(U.determinant()*V.determinant() < 0.0f) {
    Mat3d::Index id;
    svd.singularValues().minCoeff(&id);
    U.col(id)*=-1.0f;
  }
  Mat3d R=U*V.transpose();
  ASSERT(R.determinant() > 0.0f)
  return R;
}
DHMARAPEnergy::DHMARAPEnergy(const DHMCell& cell,scalarD weight)
  :DHMEnergy(cell),_weight(weight)
{
  setDegree(10);
  for(char v=0; v<8; v++)
    _pos0.col(v)=cell._verts[v]->_pos.cast<scalarD>();
}
DHMARAPEnergy::GRAD DHMARAPEnergy::grad(const Vec3d& pos) const
{
  Mat3d J0,invJ0;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::JP(Js,_pos0,J0);
  invJ0=J0.inverse();

  //build J and find conformal energy
  Mat3d J;
  DHMMapping::J(Js,_cell,J);
  J=J*invJ0;
  Eigen::JacobiSVD<Mat3d> svd(J,Eigen::ComputeFullU|Eigen::ComputeFullV);
  J=J-findR(svd);

  //build grad
  Eigen::Matrix<scalarD,9,24> dFdX;
  DHMMapping::calcDFDX(Js,invJ0,dFdX);
  return dFdX.transpose()*Eigen::Map<Eigen::Matrix<scalarD,9,1> >(J.data())*
         J0.determinant()*_weight;
}
DHMARAPEnergy::HESS DHMARAPEnergy::hess(const Vec3d& pos) const
{
  Mat3d J0,invJ0;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::JP(Js,_pos0,J0);
  invJ0=J0.inverse();

  //build hess
  Eigen::Matrix<scalarD,9,24> dFdX;
  DHMMapping::calcDFDX(Js,invJ0,dFdX);
  return (dFdX.transpose()*dFdX)*J0.determinant()*_weight;
}
//conformal energy
DHMConformalEnergy::DHMConformalEnergy(const DHMCell& cell,scalarD weight)
  :DHMEnergy(cell),_weight(weight)
{
  setDegree(10);
}
DHMConformalEnergy::GRAD DHMConformalEnergy::grad(const Vec3d& pos) const
{
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);

  //build J and find conformal energy
  Mat3d J;
  DHMMapping::J(Js,_cell,J);
  Eigen::JacobiSVD<Mat3d> svd(J,Eigen::ComputeFullU|Eigen::ComputeFullV);
  Mat3d R=findR(svd)*
          //svd.singularValues().maxCoeff();
          std::pow(svd.singularValues().prod(),(scalarD)(1.0f/3.0f));
  J-=R;

  //build grad
  GRAD ret;
  for(char vi=0; vi<8; vi++)
    ret.block<3,1>(vi*3,0)=
      J.col(0)*Js(0,vi)+J.col(1)*Js(1,vi)+J.col(2)*Js(2,vi);
  return ret*_weight;
}
DHMConformalEnergy::HESS DHMConformalEnergy::hess(const Vec3d& pos) const
{
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);

  //build hess
  HESS ret=HESS::Zero();
  for(char vi=0; vi<8; vi++)
    for(char vj=vi; vj<8; vj++) {
      scalarD val=Js.col(vi).dot(Js.col(vj));
      ret.block<3,3>(vi*3,vj*3).diagonal().setConstant(val);
      ret.block<3,3>(vj*3,vi*3).diagonal().setConstant(val);
    }
  return ret*_weight;
}
//elastic energy
DHMElasticEnergy::DHMElasticEnergy(const DHMCell& cell,scalarD lambda,scalarD mu)
  :DHMEnergy(cell)
{
  setDegree(10);
  for(char v=0; v<8; v++)
    _pos0.col(v)=cell._verts[v]->_pos.cast<scalarD>();

#define delta(A,B) (A == B ? 1 : 0)
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
        for(int l=0;l<3;l++)
          _C(i+j*3,k+l*3)=
          lambda*delta(i,j)*delta(k,l)+
          mu*(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k));
#undef delta
}
DHMElasticEnergy::GRAD DHMElasticEnergy::grad(const Vec3d& pos) const{
  Mat3d J0,invJ0;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::JP(Js,_pos0,J0);
  invJ0=J0.inverse();

  //build J and find conformal energy
  Mat3d J;
  DHMMapping::J(Js,_cell,J);
  J=J*invJ0;
  J=(J+J.transpose()).eval()*0.5f;

  Eigen::Matrix<scalarD,9,24> dFdX,dFdXT;
  DHMMapping::calcDFDX(Js,invJ0,dFdX);
  dFdXT=dFdX;
#define swap(A,B) {Eigen::Matrix<scalarD,1,24> tmp=A;A=B;B=tmp;}
  swap(dFdXT.row(1),dFdXT.row(3));
  swap(dFdXT.row(2),dFdXT.row(6));
  swap(dFdXT.row(5),dFdXT.row(7));
#undef swap
  dFdX=(dFdX+dFdXT).eval()*0.5f;
  return (dFdX.transpose()* _C*Eigen::Map<const Eigen::Matrix<scalarD,9,1> >(J.data()))*J0.determinant();
}
DHMElasticEnergy::HESS DHMElasticEnergy::hess(const Vec3d& pos) const{
  Mat3d J0,invJ0;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::JP(Js,_pos0,J0);
  invJ0=J0.inverse();

  //build hess
  Eigen::Matrix<scalarD,9,24> dFdX,dFdXT;
  DHMMapping::calcDFDX(Js,invJ0,dFdX);
  dFdXT=dFdX;
#define swap(A,B) {Eigen::Matrix<scalarD,1,24> tmp=A;A=B;B=tmp;}
  swap(dFdXT.row(1),dFdXT.row(3));
  swap(dFdXT.row(2),dFdXT.row(6));
  swap(dFdXT.row(5),dFdXT.row(7));
#undef swap
  dFdX=(dFdX+dFdXT).eval()*0.5f;
  return (dFdX.transpose()*_C*dFdX)*J0.determinant();
}
//consistency energy
DHMConsistencyEnergy::DHMConsistencyEnergy(const DHMCell& cell,scalarD weight)
  :DHMEnergy(cell),_weight(weight)
{
  setDegree(10);
  for(char v=0; v<8; v++)
    _pos0.col(v)=cell._verts[v]->_pos.cast<scalarD>();
}
DHMConsistencyEnergy::GRAD DHMConsistencyEnergy::grad(const Vec3d& pos) const
{
  Eigen::Matrix<scalarD,8,1> Ms;
  DHMMapping::MStencil(Ms,pos);

  //build M and consistency energy
  Vec3d x0=Vec3d::Zero();
  Vec3d x1=Vec3d::Zero();
  for(char v=0; v<8; v++) {
    x0+=_pos0.col(v)*Ms[v];
    x1+=_cell._verts[v]->_pos.cast<scalarD>()*Ms[v];
  }

  GRAD grad;
  for(char v=0; v<8; v++)
    grad.block<3,1>(v*3,0)=Ms[v]*(x1-x0);
  return grad*_weight;
}
DHMConsistencyEnergy::HESS DHMConsistencyEnergy::hess(const Vec3d& pos) const
{
  Eigen::Matrix<scalarD,8,1> Ms;
  DHMMapping::MStencil(Ms,pos);

  //build M and consistency energy
  HESS ret=HESS::Zero();
  for(char vi=0; vi<8; vi++)
    for(char vj=0; vj<8; vj++)
      ret.block<3,3>(vi*3,vj*3).diagonal().setConstant(Ms[vi]*Ms[vj]);
  return ret*_weight;
}
//gradient consistency energy
DHMGradConsistencyEnergy::DHMGradConsistencyEnergy(const DHMCell& cellC,const DHMCell& cellF,const Mat3d& DFDC,scalarD weight)
  :DHMEnergy(cellC),_cellF(cellF),_weight(weight),_DFDC(DFDC)
{
  setDegree(10);
}
DHMGradConsistencyEnergy::GRAD DHMGradConsistencyEnergy::grad(const Vec3d& pos) const
{
  //fine
  Mat3d JF,invJF;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,_DFDC*pos);
  DHMMapping::J(Js,_cellF,JF);
  JF*=2.0f;
  invJF=JF.inverse();

  //coarse
  Mat3d JC;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::J(Js,_cell,JC);
  JC=invJF.transpose()*(invJF*JC-Mat3d::Identity());

  //build grad
  GRAD ret;
  for(char vi=0; vi<8; vi++)
    ret.block<3,1>(vi*3,0)=
      JC.col(0)*Js(0,vi)+JC.col(1)*Js(1,vi)+JC.col(2)*Js(2,vi);
  return ret*_weight;
}
DHMGradConsistencyEnergy::HESS DHMGradConsistencyEnergy::hess(const Vec3d& pos) const
{
  //fine
  Mat3d JF,invJF;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,_DFDC*pos);
  DHMMapping::J(Js,_cellF,JF);
  JF*=2.0f;
  invJF=JF.inverse();

  //coarse
  DHMMapping::JStencil(Js,pos);

  //build hess
  HESS ret=HESS::Zero();
  Mat3d coef=invJF.transpose()*invJF;
  for(char vi=0; vi<8; vi++)
    for(char vj=vi; vj<8; vj++) {
      scalarD val=Js.col(vi).dot(Js.col(vj));
      ret.block<3,3>(vi*3,vj*3)=val*coef;
      ret.block<3,3>(vj*3,vi*3)=val*coef;
    }
  return ret*_weight;
}
//divergence
DHMDivEnergy::DHMDivEnergy(const DHMCell& cell,scalarD weight)
  :DHMEnergy(cell),_weight(weight)
{
  setDegree(10);
  for(char v=0; v<8; v++)
    _pos0.col(v)=cell._verts[v]->_pos.cast<scalarD>();
}
DHMDivEnergy::GRAD DHMDivEnergy::grad(const Vec3d& pos) const
{
  Mat3d J0,invJ0;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::JP(Js,_pos0,J0);
  invJ0=J0.inverse();

  //build grad
  Mat3d J=Mat3d::Identity();
  Eigen::Matrix<scalarD,9,24> dFdX;
  DHMMapping::calcDFDX(Js,invJ0,dFdX);
  return dFdX.transpose()*Eigen::Map<Eigen::Matrix<scalarD,9,1> >(J.data())*
         J0.determinant()*_weight;
}
DHMDivEnergy::HESS DHMDivEnergy::hess(const Vec3d& pos) const
{
  HESS ret;
  ret.setZero();
  return ret;
}

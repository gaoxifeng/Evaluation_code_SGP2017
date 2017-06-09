#ifndef DHM_ENERGY_H
#define DHM_ENERGY_H

#include "DHMTraits.h"
#include "DHMMesh.h"
#include "GaussLegendre.h"
#include "CommonFile/ImplicitFuncInterface.h"

PRJ_BEGIN

template <int rows>
class DHMEnergy : public DHMEnergyTraits<scalarD,rows>
{
public:
  using typename DHMEnergyTraits<scalarD,rows>::COL;
  using typename DHMEnergyTraits<scalarD,rows>::HESS;
  using typename DHMEnergyTraits<scalarD,rows>::GRAD;
  struct HessFunc {
    HessFunc(const DHMEnergy& energy):_energy(energy) {}
    HESS operator()(const Vec3d& pos) const {
      return _energy.hess(pos);
    }
    const DHMEnergy& _energy;
  };
  struct GradFunc {
    GradFunc(const DHMEnergy& energy):_energy(energy) {}
    GRAD operator()(const Vec3d& pos) const {
      return _energy.grad(pos);
    }
    const DHMEnergy& _energy;
  };
  DHMEnergy(const DHMCell& cell):_cell(cell) {}
  virtual ~DHMEnergy() {}
  virtual HESS hessInt() const {
    HessFunc f(*this);
    return GaussLegendreIntegral<HESS,scalarD>::integrate3D(f,_deg);
  }
  virtual GRAD gradInt() const {
    GradFunc f(*this);
    return GaussLegendreIntegral<GRAD,scalarD>::integrate3D(f,_deg);
  }
  virtual HESS hess(const Vec3d& pos) const=0;
  virtual GRAD grad(const Vec3d& pos) const=0;
  const DHMCell& getCell() const {
    return _cell;
  }
  void setDegree(sizeType deg) {
    _deg=deg;
  }
  template <int MAXDEG>
  sizeType findDegree() {
    Eigen::Matrix<scalar,3,8> pos;
    for(char v=0; v<8; v++) {
      Vec3& p=_cell._verts[v]->_pos;
      pos.col(v)=p;
      p.setRandom();
    }
    vector<HESS,Eigen::aligned_allocator<HESS> > h;
    vector<GRAD,Eigen::aligned_allocator<GRAD> > g;
    for(sizeType deg=0; deg<=MAXDEG; deg++) {
      setDegree(deg);
      h.push_back(hessInt());
      g.push_back(gradInt());
      if(deg > 0 && (h[deg]-h[deg-1]).norm() < 1E-9f) { // && (g[deg]-g[deg-1]).norm() < 1E-9f)
        INFOV("Degree of %s=%ld",typeid(*this).name(),deg-1)
        return deg-1;
      }
    }
    INFOV("Degree of %s>=%d",typeid(*this).name(),MAXDEG-1)
    return 0;
  }
  void debugHessGrad() {
    Vec3d pos=Vec3d::Random();
    HESS h=hess(pos);
    GRAD g=grad(pos);
    for(char v=0; v<8; v++)
      for(char d=0; d<3; d++) {
        const scalar& p=_cell._verts[v]->_pos[d];
        scalar pTmp=p;
#define DELTA 1E-7f
        const_cast<scalar&>(p)+=DELTA;
        GRAD hN=(grad(pos)-g)/DELTA;
        INFOV("Hess Err: %f %f",h.col(v*3+d).norm(),(h.col(v*3+d)-hN).norm())
        const_cast<scalar&>(p)=pTmp;
#undef DELTA
      }
  }
  void debugHessGrad(Eigen::Matrix<scalarD,rows,1>& val) {
    Vec3d pos=Vec3d::Random();
    HESS h=hess(pos);
    GRAD g=grad(pos);
    for(sizeType v=0; v<val.size(); v++) {
      const scalarD& p=val[v];
      scalarD pTmp=p;
#define DELTA 1E-7f
      const_cast<scalarD&>(p)+=DELTA;
      GRAD hN=(grad(pos)-g)/DELTA;
      INFOV("Hess Err: %f %f",h.col(v).norm(),(h.col(v)-hN).norm())
      const_cast<scalarD&>(p)=pTmp;
#undef DELTA
    }
  }
  void debugHessGrad(COL& val) {
    Vec3d pos=Vec3d::Random();
    HESS h=hess(pos);
    GRAD g=grad(pos);
    for(sizeType v=0; v<8; v++) {
      const scalarD& p=val[_cell._verts[v]->_index];
      scalarD pTmp=p;
#define DELTA 1E-7f
      const_cast<scalarD&>(p)+=DELTA;
      GRAD hN=(grad(pos)-g)/DELTA;
      INFOV("Hess Err: %f %f",h.col(v).norm(),(h.col(v)-hN).norm())
      const_cast<scalarD&>(p)=pTmp;
#undef DELTA
    }
  }
protected:
  const DHMCell& _cell;
  sizeType _deg;
};
class DHMEnergyPool : public DHMEnergyTraits<scalarD,0>
{
  typedef vector<boost::shared_ptr<DHMEnergy<24> > > VENERGYS;
public:
  static void debugEnergy();
  //add energy
  template <typename T>
  T& getVDEnergy(sizeType i) {
    return *boost::dynamic_pointer_cast<T>(_VEnergys[i]);
  }
  void addVDEnergy(boost::shared_ptr<DHMEnergy<24> > VEnergy);
  void clearEnergy();
  //build
  void buildHessian(TRIPS& trips);
  void buildHessian(SMAT& hess);
  void buildDeriv(COL& deriv);
  void addDiagonal(SMAT& hess,scalar eps);
  void ensureDD(SMAT& hess);
  //fixed
  boost::unordered_set<sizeType> _fixed;
  VENERGYS _VEnergys;
};
//arap energy
class DHMARAPEnergy : public DHMEnergy<24>
{
public:
  using typename DHMEnergy<24>::COL;
  using typename DHMEnergy<24>::HESS;
  using typename DHMEnergy<24>::GRAD;
  DHMARAPEnergy(const DHMCell& cell,scalarD weight);
  GRAD grad(const Vec3d& pos) const;
  HESS hess(const Vec3d& pos) const;
  Eigen::Matrix<scalarD,3,8> _pos0;
  scalarD _weight;
};
//conformal energy
class DHMConformalEnergy : public DHMEnergy<24>
{
public:
  using typename DHMEnergy<24>::COL;
  using typename DHMEnergy<24>::HESS;
  using typename DHMEnergy<24>::GRAD;
  DHMConformalEnergy(const DHMCell& cell,scalarD weight);
  GRAD grad(const Vec3d& pos) const;
  HESS hess(const Vec3d& pos) const;
  scalarD _weight;
};
//elastic energy
class DHMElasticEnergy : public DHMEnergy<24>
{
public:
  using typename DHMEnergy<24>::COL;
  using typename DHMEnergy<24>::HESS;
  using typename DHMEnergy<24>::GRAD;
  DHMElasticEnergy(const DHMCell& cell, scalarD lambda, scalarD mu);
  GRAD grad(const Vec3d& pos) const;
  HESS hess(const Vec3d& pos) const;
  Eigen::Matrix<scalarD,3,8> _pos0;
  Eigen::Matrix<scalarD,9,9> _C;
};
//consistency energy
class DHMConsistencyEnergy : public DHMEnergy<24>
{
public:
  using typename DHMEnergy<24>::COL;
  using typename DHMEnergy<24>::HESS;
  using typename DHMEnergy<24>::GRAD;
  DHMConsistencyEnergy(const DHMCell& cell,scalarD weight);
  GRAD grad(const Vec3d& pos) const;
  HESS hess(const Vec3d& pos) const;
  Eigen::Matrix<scalarD,3,8> _pos0;
  scalarD _weight;
};
//gradient consistency energy
class DHMGradConsistencyEnergy : public DHMEnergy<24>
{
public:
  using typename DHMEnergy<24>::COL;
  using typename DHMEnergy<24>::HESS;
  using typename DHMEnergy<24>::GRAD;
  DHMGradConsistencyEnergy(const DHMCell& cellC,const DHMCell& cellF,const Mat3d& DFDC,scalarD weight);
  GRAD grad(const Vec3d& pos) const;
  HESS hess(const Vec3d& pos) const;
  const DHMCell& _cellF;
  scalarD _weight;
  Mat3d _DFDC;
};
//divergence
class DHMDivEnergy : public DHMEnergy<24>
{
public:
  using typename DHMEnergy<24>::COL;
  using typename DHMEnergy<24>::HESS;
  using typename DHMEnergy<24>::GRAD;
  DHMDivEnergy(const DHMCell& cell,scalarD weight);
  GRAD grad(const Vec3d& pos) const;
  HESS hess(const Vec3d& pos) const;
  Eigen::Matrix<scalarD,3,8> _pos0;
  scalarD _weight;
};

PRJ_END

#endif

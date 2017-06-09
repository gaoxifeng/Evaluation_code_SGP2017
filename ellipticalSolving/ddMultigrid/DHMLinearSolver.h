#ifndef DHM_LINEAR_SOLVER_H
#define DHM_LINEAR_SOLVER_H

#include "DHMField.h"
#include "DHMCoarser.h"
#include "DHMEnergy.h"

PRJ_BEGIN

//DHMLinearSystem
struct DHMLinearSystem : public DHMTraits, public DHMEnergyTraits<scalarD> {
  typedef Eigen::Matrix<scalarD,8,8> CONS;
  struct Filter {
    bool operator()(sizeType i) const {
      return i>=(sizeType)_valid.size() || _valid[i];
    }
    vector<bool> _valid;
  };
  //basic iterative solver
  DHMLinearSystem(const DHMMesh& mesh,const SMAT* prolong);
  void buildRBTag();
  void smooth(const COL& rhs,COL& x,const Filter& f=Filter()) const;
  void dampedJacobi(const COL& rhs,COL& x,const Filter& f=Filter()) const;
  void gaussSeidel(const COL& rhs,COL& x,const Filter& f=Filter()) const;
  void RBGS(const COL& rhs,COL& x,char tag,const Filter& f=Filter()) const;
  void solve(const COL& rhs,COL& x,scalarD thres=1E-5f,sizeType nrIterInner=10) const;
  //system builder
  template <typename ENERGY>
  void buildSystem() {
    TRIPS trips;
    sizeType nrC=_mesh.nrCell();
    trips.resize(nrC*64);

    CELLID vids;
    typename ENERGY::HESS hess,cons;
    OMP_PARALLEL_FOR_I(OMP_PRI(vids,hess,cons))
    for(sizeType i=0; i<nrC; i++) {
      //calculate hessian
      const DHMCell& cell=_mesh.getCell(i);
      hess=ENERGY(cell).hessInt();
      //handle constraints
      if(findVid(vids,cons,cell))
        hess=cons.transpose()*(hess*cons).eval();
      //insert entries
      sizeType off=i*64;
      for(char r=0; r<8; r++)
        for(char c=0; c<8; c++)
          trips[off++]=Eigen::Triplet<scalarD,sizeType>(vids[r],vids[c],hess(r,c));
    }

    sizeType nrVertNC=_mesh.nrVertNC();
    _lhs.resize(nrVertNC);
    _lhs.buildFromTripletsDepulicate(trips,0.0f);
    //ASSERT(sys._lhs.isSymmetric(1E-9f));
  }
  void lumpSystem();
  static void insertVid(sizeType vid,CELLID& vids,sizeType r,CONS& cons,scalar coef);
  static bool findVid(CELLID& vids,CONS& cons,const DHMCell& cell);
  //data
  const DHMMesh& _mesh;
  const SMAT* _prolong;
  SMAT _lhs;
  //tagging
  vector<sizeType> _tag,_tagOff;
};
//Basic Solver
class DHMPressureSolver : public DHMTraits, public DHMEnergyTraits<scalarD>
{
public:
  DHMPressureSolver(const DHMMesh& mesh,DHMCVField& vel,const DHMVSField* lv=NULL);
  void precomputeDiv(FixedSparseMatrix<scalarD,Kernel<scalarD> >& rhs) const;
  const DHMLinearSystem& getSystem() const;
  virtual const DHMField<scalarD,false,1>& getPressure() const;
  virtual void precomputeLHS();
  virtual void project();
protected:
  static void precomputeRHSLevel(const DHMLinearSystem& sys,const DHMCVField& vel,const DHMField<scalarD,false,1>& pre,COL& RHS);
  virtual void solveLaplacian();
  //data
  vector<boost::shared_ptr<DHMLinearSystem> > _sys;
  DHMCVField& _vel;
  DHMField<scalarD,false,1> _pre;
  sizeType _maxIter;
  scalarD _relTol;
};
//Multigrid Solver
class DHMPressureSolverMG : public DHMPressureSolver
{
public:
  DHMPressureSolverMG(const Hierarchy& hier,DHMCVField& vel,const DHMVSField* lv=NULL);
  virtual void precomputeLHS();
protected:
  virtual void solveLaplacian();
  //data
  const Hierarchy& _hier;
};
//Dirichlet's Energy
class DHMDirichletEnergy : public DHMEnergy<8>
{
public:
  DHMDirichletEnergy(const DHMCell& cell,const DHMField<scalarD,false,1>* pre=NULL,const GRAD& comp=GRAD::Zero());
  GRAD grad(const Vec3d& pos) const;
  HESS hess(const Vec3d& pos) const;
  const DHMField<scalarD,false,1>* _pre;
  const GRAD _comp;
};

PRJ_END

#endif

#include "DHMLinearSolver.h"
#include "CommonFile/solvers/LinearSolver.h"

USE_PRJ_NAMESPACE

//DHMLinearSystem
DHMLinearSystem::DHMLinearSystem(const DHMMesh& mesh,const SMAT* prolong)
  :_mesh(mesh),_prolong(prolong) {}
void DHMLinearSystem::buildRBTag()
{
  //initialize
  _tag.clear();
  _tagOff.clear();
  vector<bool> visited(_mesh.nrVertNC(),false);

  //build tag
  _tag.push_back(0);
  _tagOff.push_back(0);
  _tagOff.push_back(1);
  visited[0]=true;
  while(true) {
    sizeType beg=_tagOff[(sizeType)_tagOff.size()-2];
    sizeType end=_tagOff[(sizeType)_tagOff.size()-1];
    for(sizeType i=beg; i<end; i++) {
      sizeType row=_tag[i];
      for(ConstSMIterator<scalarD> beg=_lhs.begin(row),end=_lhs.end(row); beg!=end; ++beg)
        if(!visited[beg.col()]) {
          _tag.push_back(beg.col());
          visited[beg.col()]=true;
        }
    }
    _tagOff.push_back((sizeType)_tag.size());
    if(_tag.size() == visited.size())
      break;
  }

  //debug write
  //DHMVSField tagField(_mesh);
  //sizeType nrTagOff=(sizeType)_tagOff.size();
  //for(char tag=0;tag<2;tag++)
  //	for(sizeType g=tag;g<nrTagOff-1;g+=2)
  //		for(sizeType row=_tagOff[g];row<_tagOff[g+1];row++)
  //			tagField.getCompRef(_tag[row])[0]=tag;
  //tagField.writeSVTK("./tag.vtk");
}
void DHMLinearSystem::smooth(const COL& rhs,COL& x,const Filter& f) const
{
  if(_tag.empty()) {
    //dampedJacobi(rhs,x,f);
    gaussSeidel(rhs,x,f);
  } else {
    RBGS(rhs,x,0,f);
    RBGS(rhs,x,1,f);
  }
}
void DHMLinearSystem::dampedJacobi(const COL& rhs,COL& x,const Filter& f) const
{
  COL xBK;
  xBK.setZero(x.size());

  const std::vector<sizeType>& rowStart=_lhs.getRowOffset();
  const SMAT::ROW& value=_lhs.getValue();

  scalarD D,w=2.0f/3.0f;
  sizeType i,j,C,nrRow=_lhs.rows();
  OMP_PARALLEL_FOR_I(OMP_PRI(i,j,C,D))
  for(i=0; i<nrRow; i++)
    if(f(i)) {
      D=0.0f;
      xBK[i]=rhs[i];
      for(j=rowStart[i]; j<rowStart[i+1]; ++j) {
        C=value[j].second;
        if(C == i)
          D=value[j].first;
        xBK[i]-=value[j].first*x[C];
      }
      if(D > 0.0f) xBK[i]*=(w/D);
      else xBK[i]=0.0f;
    }
  x+=xBK;
}
void DHMLinearSystem::gaussSeidel(const COL& rhs,COL& x,const Filter& f) const
{
  const std::vector<sizeType>& rowStart=_lhs.getRowOffset();
  const SMAT::ROW& value=_lhs.getValue();

  scalarD D;
  sizeType i,j,C,nrRow=_lhs.rows();
  for(i=0; i<nrRow; i++)
    if(f(i)) {
      D=0.0f;
      x[i]=rhs[i];
      for(j=rowStart[i]; j<rowStart[i+1]; ++j) {
        C=value[j].second;
        if(C == i) D=value[j].first;
        else x[i]-=value[j].first*x[C];
      }
      if(D > 0.0f)x[i]/=D;
      else x[i]=0.0f;
    }
}
void DHMLinearSystem::RBGS(const COL& rhs,COL& x,char tag,const Filter& f) const
{
  const std::vector<sizeType>& rowStart=_lhs.getRowOffset();
  const SMAT::ROW& value=_lhs.getValue();

  scalarD D;
  sizeType g,i,j,row,C;
  sizeType nrTagOff=(sizeType)_tagOff.size();
  OMP_PARALLEL_FOR_I(OMP_PRI(g,i,j,row,C,D))
  for(g=tag; g<nrTagOff-1; g+=2) {
    for(row=_tagOff[g]; row<_tagOff[g+1]; row++)
      if(f(i=_tag[row])) {
        D=0.0f;
        x[i]=rhs[i];
        for(j=rowStart[i]; j<rowStart[i+1]; ++j) {
          C=value[j].second;
          if(C == i) D=value[j].first;
          else x[i]-=value[j].first*x[C];
        }
        if(D > 0.0f)x[i]/=D;
        else x[i]=0.0f;
      }
  }
}
void DHMLinearSystem::solve(const COL& rhs,COL& x,scalarD thres,sizeType nrIterInner) const
{
  COL res;
  _lhs.multiplySubtract(x,res=rhs);
  scalarD tol=std::max<scalarD>(res.norm()*thres,1E-9f);
  sizeType it=0;
  for(; res.norm() > tol; it++) {
    for(sizeType v=0; v<nrIterInner; v++)
      smooth(rhs,x);
    _lhs.multiplySubtract(x,res=rhs);
  }
  INFOV("Solve in %ld iters",it*nrIterInner);
}
//system builder
void DHMLinearSystem::lumpSystem()
{
  COL diag=COL::Zero(_lhs.rows());
  for(sizeType r=0; r<_lhs.rows(); r++)
    for(ConstSMIterator<scalarD> beg=_lhs.begin(r),end=_lhs.end(r); beg!=end; ++beg)
      diag[r]+=*beg;

  TRIPS trips;
  addI(trips,0,0,diag);
  _lhs.buildFromTriplets(trips.begin(),trips.end());
}
void DHMLinearSystem::insertVid(sizeType vid,CELLID& vids,sizeType r,CONS& cons,scalar coef)
{
  //search for entry
  char sameSlot=-1,firstSlot=-1;
  for(char v=0; v<8; v++)
    if(vids[v] == vid)
      sameSlot=v;
    else if(firstSlot == -1 && vids[v] == numeric_limits<sizeType>::max())
      firstSlot=v;
  if(sameSlot == -1) {
    ASSERT(firstSlot >= 0)
    sameSlot=firstSlot;
    vids[sameSlot]=vid;
  }
  cons(r,sameSlot)=coef;
}
bool DHMLinearSystem::findVid(CELLID& vids,CONS& cons,const DHMCell& cell)
{
  boost::shared_ptr<DHMConstrainedVertex> cv;
  cons.setZero();
  vids.setConstant(numeric_limits<sizeType>::max());
  bool isCon=false;
  for(char v=0; v<8; v++)
    if(cv=boost::dynamic_pointer_cast<DHMConstrainedVertex>(cell._verts[v])) {
      isCon=true;
      for(char c=0; c<4; c++)
        if(cv->_verts[c])
          insertVid(cv->_verts[c]->_index,vids,v,cons,cv->_weights[c]);
    } else insertVid(cell._verts[v]->_index,vids,v,cons,1.0f);
  for(char v=0; v<8; v++)
    ASSERT(vids[v] != numeric_limits<sizeType>::max());
  return isCon;
}
//Basic Solver
class LaplacianNullspace : public Nullspace<Kernel<scalarD> >
{
public:
  virtual void projectOut(Vec& r) {
    Vec delta=Vec::Ones(r.size());
    delta*=r.sum()/delta.sum();
    r-=delta;
  }
};
DHMPressureSolver::DHMPressureSolver(const DHMMesh& mesh,DHMCVField& vel,const DHMVSField* lv)
  :_vel(vel),_pre(mesh),_maxIter(10000),_relTol(1E-6f)
{
  _sys.push_back(boost::shared_ptr<DHMLinearSystem>(new DHMLinearSystem(mesh,NULL)));
}
void DHMPressureSolver::precomputeDiv(FixedSparseMatrix<scalarD,Kernel<scalarD> >& rhs) const
{
  const DHMMesh& mesh=_sys[0]->_mesh;
  const SMAT& lhs=_sys[0]->_lhs;

  vector<Eigen::Triplet<scalarD,sizeType> >trips;
  for(sizeType i=0;i<mesh.nrCell();i++)
  {
    const DHMCell& cell=mesh.getCell(i);
    DHMDivEnergy div(cell,1);
    DHMDivEnergy::GRAD grad=div.gradInt();
    for(int v=0;v<8;v++)
      addBlock(trips,i,cell._verts[v]->_index*3,grad.segment<3>(v*3).transpose());
  }
  rhs.resize(mesh.nrCell(),lhs.rows()*3);
  rhs.buildFromTriplets(trips.begin(),trips.end());
}
const DHMLinearSystem& DHMPressureSolver::getSystem() const
{
  return *_sys[0];
}
const DHMField<scalarD,false,1>& DHMPressureSolver::getPressure() const
{
  return _pre;
}
void DHMPressureSolver::precomputeLHS()
{
  _sys.front()->buildSystem<DHMDirichletEnergy>();
}
void DHMPressureSolver::project()
{
  //initialize solver
  solveLaplacian();
  //remove gradient from velocity
  const DHMMesh& mesh=_pre._mesh;
  Eigen::Matrix<scalar,8,1> cellPre;
  for(sizeType i=0; i<mesh.nrCell(); i++) {
    const DHMCell& cell=mesh.getCell(i);
    for(char v=0; v<8; v++)
      cellPre[v]=_pre.getVScalar<scalar>(0,*(cell._verts[v]));
    _vel.getCompRef(i)-=cellPre;
  }
}
void DHMPressureSolver::precomputeRHSLevel(const DHMLinearSystem& sys,const DHMCVField& vel,const DHMField<scalarD,false,1>& pre,COL& RHS)
{
  sizeType nrC=sys._mesh.nrCell();
  sizeType nrVertNC=sys._mesh.nrVertNC();
  RHS.setZero(nrVertNC);

  CELLID vids;
  DHMDirichletEnergy::GRAD grad;
  DHMDirichletEnergy::HESS cons;
  OMP_PARALLEL_FOR_I(OMP_PRI(vids,grad,cons))
  for(sizeType i=0; i<nrC; i++) {
    //calculate hessian
    const DHMCell& cell=sys._mesh.getCell(i);
    grad=DHMDirichletEnergy(cell,&pre,vel.getComp(i).cast<scalarD>()).gradInt();
    //handle constraints
    if(DHMLinearSystem::findVid(vids,cons,cell))
      grad=cons.transpose()*grad;
    //insert entries
    for(char v=0; v<8; v++) {
      scalarD& coeff=RHS[vids[v]];
      OMP_ATOMIC_
      coeff-=grad[v];
    }
  }
}
void DHMPressureSolver::solveLaplacian()
{
  //setup solver
  DHMLinearSystem& sys=*(_sys.front());
  LaplacianNullspace ns;

  PCGSolver<scalarD,Kernel<scalarD>,NoPreconSolver<scalarD> > sol;
  sol.setSolverParameters(_relTol,_maxIter);
  sol.setMatrix(sys._lhs,true);
  sol.setNullspace(&ns);
  sol.setCallback(boost::shared_ptr<Callback<scalarD,Kernel<scalarD> > >(new Callback<scalarD,Kernel<scalarD> > ));

  //setup RHS
  COL rhs,x;
  _pre.clear(0.0f);
  precomputeRHSLevel(sys,_vel,_pre,rhs);
  x.resize(rhs.size());

  //solve
  sol.solve(rhs,x);
  _pre._data.block(0,0,x.size(),1)=x;
  INFOV("RHS: %f, Nr Iter: %ld",rhs.norm(),sol.getIterationsCount());
}
//Multigrid Solver
class DHMMGPreconditioner : public Solver<scalarD,Kernel<scalarD> >, public DHMEnergyTraits<scalarD>
{
public:
  DHMMGPreconditioner() {
    //multigrid setup
    _nrS=2;
    _nrSF=10;
  }
  virtual void setMatrix(const SMAT& matrix,bool syncPrecon) {
    ASSERT(false)
  }
  virtual void setMatrix(const vector<boost::shared_ptr<DHMLinearSystem> >& sys) {
    sizeType nrLv=(sizeType)sys.size();
    _rhs.resize(nrLv-1);
    _x.resize(nrLv-1);

    _sys=&sys;
    for(sizeType i=0; i<nrLv-1; i++) {
      _rhs[i].resize((*_sys)[i+1]->_mesh.nrVertNC());
      _x[i].resize((*_sys)[i+1]->_mesh.nrVertNC());
    }

    setFilter(true);
  }
  virtual SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
    result.setZero();
    VCycle(0,rhs,result);
    return SUCCESSFUL;
  }
  //multigrid solver cycle
  void VCycle(sizeType id,const COL& rhs,COL& x) {
    DHMLinearSystem& sys=*((*_sys)[id]);
    if(id < (sizeType)_sys->size()-1) {
      //pre smoothing
      for(sizeType i=0; i<_nrS; i++)
        sys.smooth(rhs,x,(id>0 && id<=_nrALv) ? _filter : DHMLinearSystem::Filter());
      //compute residual
      COL res=rhs;
      sys._lhs.multiplySubtract(x,res);
      //restrict
      DHMLinearSystem& sysN=*((*_sys)[id+1]);
      sys._prolong->multiplyTranspose(res,_rhs[id]);
      //descend
      _x[id].setZero();
      VCycle(id+1,_rhs[id],_x[id]);
      //prolong add
      sys._prolong->multiplyAdd(_x[id],x);
      //post smoothing
      for(sizeType i=0; i<_nrS; i++)
        sys.smooth(rhs,x,(id<_nrALv) ? _filter : DHMLinearSystem::Filter());
    } else {
      //PMINRESSolver<scalarD> sol;
      //sol.setSolverParameters(_relTol,_maxIter);
      //sol.setMatrix(sys._lhs,true);
      //sol.solve(sys._rhs,sys._x);

      //final smoothing
      for(sizeType i=0; i<_nrSF; i++)
        sys.smooth(rhs,x);
    }
  }
  void setFilter(bool useFilter) {
    _nrALv=0;
    if(useFilter) {
      const vector<boost::shared_ptr<DHMLinearSystem> >& sys=*_sys;
      while(sys[_nrALv]->_prolong && sys[_nrALv]->_prolong->rows() < sys[_nrALv]->_mesh.nrVert())
        _nrALv++;

      const DHMMesh& baseMesh=sys[_nrALv]->_mesh;
      sizeType nrBaseVert=baseMesh.nrVert();
      vector<bool>& valid=_filter._valid;
      valid.assign(nrBaseVert,false);

      sizeType nrSubCell=(1<<_nrALv)*(1<<_nrALv)*(1<<_nrALv);
      const DHMMesh& finestMesh=sys.front()->_mesh;
      const vector<sizeType>& adaptiveTag=finestMesh.getAdaptive();
      for(sizeType i=0; i<(sizeType)adaptiveTag.size(); i++)
        if(adaptiveTag[i] < 0)
          for(sizeType c=0,off=-1-adaptiveTag[i]; c<nrSubCell; c++,off++) {
            const DHMCell& cell=finestMesh.getCell(off);
            for(char v=0; v<8; v++)
              if(!boost::dynamic_pointer_cast<DHMConstrainedVertex>(cell._verts[v]) && cell._verts[v]->_index < nrBaseVert)
                valid[cell._verts[v]->_index]=true;
          }

      for(sizeType iter=0; iter<_nrALv*_nrS; iter++) {
        vector<bool> validLast=valid;
        for(sizeType c=0; c<baseMesh.nrCell(); c++) {
          const DHMCell& cell=baseMesh.getCell(c);
          bool isNeigh=false;
          for(char v=0; v<8; v++)
            if(validLast[cell._verts[v]->_index]) {
              isNeigh=true;
              break;
            }
          if(isNeigh)
            for(char v=0; v<8; v++)
              valid[cell._verts[v]->_index]=true;
        }
      }

      //debug write
      //DHMVSField tagField(baseMesh);
      //for(sizeType i=0;i<nrBaseVert;i++)
      //	tagField.getCompRef(i)[0]=valid[i]?1.0f:0.0f;
      //tagField.writeSVTK("./tag.vtk");
    }
  }
  //data
  const vector<boost::shared_ptr<DHMLinearSystem> >* _sys;
  vector<COL,Eigen::aligned_allocator<COL> > _rhs;
  vector<COL,Eigen::aligned_allocator<COL> > _x;
  sizeType _nrS,_nrSF,_nrALv;
  DHMLinearSystem::Filter _filter;
};
DHMPressureSolverMG::DHMPressureSolverMG(const Hierarchy& hier,DHMCVField& vel,const DHMVSField* lv)
  :DHMPressureSolver(*(hier._mesh.front()),vel),_hier(hier)
{
  sizeType nrMesh=(sizeType)_hier._mesh.size();
  _sys.clear();
  for(sizeType i=0; i<nrMesh; i++) {
    const DHMMesh& mesh=*(_hier._mesh[i]);
    const SMAT* prolong=(i<nrMesh-1) ? &(_hier._prolong[i]) : NULL;
    _sys.push_back(boost::shared_ptr<DHMLinearSystem>(new DHMLinearSystem(mesh,prolong)));
  }
}
void DHMPressureSolverMG::precomputeLHS()
{
  for(sizeType i=0; i<(sizeType)_sys.size(); i++) {
    _sys[i]->buildSystem<DHMDirichletEnergy>();
    _sys[i]->buildRBTag();
  }
}
void DHMPressureSolverMG::solveLaplacian()
{
  DHMLinearSystem& sys=*(_sys.front());
  LaplacianNullspace ns;

  PCGSolver<scalarD,Kernel<scalarD>,DHMMGPreconditioner> sol;
  sol.setSolverParameters(_relTol,_maxIter);
  sol.setMatrix(sys._lhs,false);
  dynamic_cast<DHMMGPreconditioner*>(sol.getPre())->setMatrix(_sys);
  sol.setUseIPCG(true);
  sol.setNullspace(&ns);
  sol.setCallback(boost::shared_ptr<Callback<scalarD,Kernel<scalarD> > >(new Callback<scalarD,Kernel<scalarD> > ));

  //setup RHS
  COL rhs,x;
  _pre.clear(0.0f);
  precomputeRHSLevel(sys,_vel,_pre,rhs);
  x.resize(rhs.size());

  //solve
  sol.solve(rhs,x);
  _pre._data.block(0,0,x.size(),1)=x;
  INFOV("RHS: %f, Nr Iter: %ld",rhs.norm(),sol.getIterationsCount());
}
//Dirichlet's Energy, gradient and hessian
DHMDirichletEnergy::DHMDirichletEnergy(const DHMCell& cell,const DHMField<scalarD,false,1>* pre,const GRAD& comp)
  :DHMEnergy(cell),_pre(pre),_comp(comp)
{
  setDegree(10);
}
DHMDirichletEnergy::GRAD DHMDirichletEnergy::grad(const Vec3d& pos) const
{
  Mat3d JF;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::J(Js,_cell,JF);
  Js=JF.inverse().eval().transpose()*Js;

  GRAD compPre;
  for(char v=0; v<8; v++)
    compPre[v]=_pre->getVScalar<scalarD>(0,*(_cell._verts[v]));
  return (Js.transpose()*(Js*(compPre-_comp)))*JF.determinant();
}
DHMDirichletEnergy::HESS DHMDirichletEnergy::hess(const Vec3d& pos) const
{
  Mat3d JF;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::J(Js,_cell,JF);
  Js=JF.inverse().eval().transpose()*Js;
  return (Js.transpose()*Js)*JF.determinant();
}

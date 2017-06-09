#include "DHMSmoother.h"
#include "DHMEnergy.h"
#include "DHMUtil.h"

USE_PRJ_NAMESPACE

DHMSmoother::DHMSmoother(DHMMesh& mesh)
  :_mesh(mesh)
{
  _smesh=mesh.getSurface(&_surf);
  setParameter();
}
void DHMSmoother::setParameter()
{
  putNoOverwrite<sizeType>(*(_mesh._tree),"nrLaplacian",20);
  putNoOverwrite<scalar>(*(_mesh._tree),"dtLaplacian",0.1f);

  putNoOverwrite<sizeType>(*(_mesh._tree),"maxIterSmooth",10000);
  putNoOverwrite<scalar>(*(_mesh._tree),"relThresSmooth",1E-3f);
  putNoOverwrite<scalar>(*(_mesh._tree),"conformalWeightSmooth",0.0f);
  putNoOverwrite<scalar>(*(_mesh._tree),"consistencyWeightSmooth",1.0f);
}
bool DHMSmoother::smooth()
{
  //initialize
  DHMEnergyPool E;
  for(sizeType i=0; i<_mesh.nrCell(); i++) {
    if(_mesh._tree->get<bool>("ARAPInteriorSmooth",true))
      E.addVDEnergy(boost::shared_ptr<DHMARAPEnergy>(new DHMARAPEnergy(_mesh.getCell(i),1.0f)));
    else E.addVDEnergy(boost::shared_ptr<DHMConformalEnergy>(new DHMConformalEnergy(_mesh.getCell(i),1.0f)));
  }
  for(boost::unordered_map<sizeType,sizeType>::const_iterator
      beg=_surf.begin(),end=_surf.end(); beg!=end; beg++)
    E._fixed.insert(beg->first);

  //surface smooth
  sizeType nrLap=_mesh._tree->get<sizeType>("nrLaplacian");
  for(sizeType i=0; i<nrLap; i++) {
    smoothSurface();
    smoothInterior(E);
    if(_mesh._tree->get<bool>("writeLaplacian",true)) {
      ostringstream oss;
      oss << "./lap" << i << ".vtk";
      _mesh.writeVTK(oss.str());
    }
  }
  return makeValid();
}
//make flip free
bool DHMSmoother::updateWeight(DHMEnergyPool& EC) const
{
  sizeType nrInvalid=0;
  for(sizeType i=0; i<_mesh.nrCell(); i++)
    if(!DHMUniformValidityBound::isValid(_mesh.getCell(i))) {
      scalarD& w=EC.getVDEnergy<DHMConformalEnergy>(i)._weight;
      w=std::max<scalarD>(w,0.01f)*2.0f;
      nrInvalid++;
    }
  INFOV("Found %d invalid elements!",nrInvalid)
  return nrInvalid>0;
}
bool DHMSmoother::makeValid()
{
  //initialize
  sizeType maxIter=_mesh._tree->get<sizeType>("maxIterSmooth");
  scalar thres=_mesh._tree->get<scalar>("relThresSmooth");
  scalar conformalWeight=_mesh._tree->get<scalar>("conformalWeightSmooth");
  scalar consistencyWeight=_mesh._tree->get<scalar>("consistencyWeightSmooth");

  //build conformal/consistency energy
  DHMEnergyPool E;
  for(sizeType i=0; i<_mesh.nrCell(); i++)
    E.addVDEnergy(boost::shared_ptr<DHMConformalEnergy>(new DHMConformalEnergy(_mesh.getCell(i),conformalWeight)));
  for(sizeType i=0; i<_mesh.nrCell(); i++)
    E.addVDEnergy(boost::shared_ptr<DHMConsistencyEnergy>(new DHMConsistencyEnergy(_mesh.getCell(i),consistencyWeight)));

  //build hessian
  SMAT hess;
  hess.resize(_mesh.nrVert()*3,_mesh.nrVert()*3);
  E.buildHessian(hess);
  _sol.setSolverParameters(1E-6f,10000);
  _sol.setMatrix(hess,true);

  //main loop
  COL deriv,delta;
  scalarD deltaMax=0.0f;
  INFO("Make Valid Main Loop")
  for(sizeType i=0; i<maxIter; i++) {
    //calculate hessian and derivative
    deriv.setZero(_mesh.nrVert()*3);
    E.buildDeriv(deriv);

    //solve system for delta
    COL pos;
    _mesh.getPos(pos);
    delta.setZero(_mesh.nrVert()*3);
    _sol.solve(deriv,delta);
    pos-=delta;
    _mesh.setPos(pos);

    //check for termination condition
    if(deltaMax == 0.0f && delta.cwiseAbs().maxCoeff() > 0.0f)
      deltaMax=delta.cwiseAbs().maxCoeff();
    INFOV("Delta %f/%f solved in %d iters",delta.cwiseAbs().maxCoeff(),deltaMax,_sol.getIterationsCount())
    if(deltaMax == 0.0f || delta.cwiseAbs().maxCoeff() < thres*deltaMax) {
      if(!updateWeight(E))
        return true;
      else {
        E.buildHessian(hess);
        _sol.setMatrix(hess,true);
      }
    }
  }
  return false;
}
//smooth surface/interior
void normalizeLaplacian(FixedSparseMatrix<scalarD,Kernel<scalarD> >& m)
{
  for(sizeType r=0; r<m.rows(); r++) {
    scalarD diag=fabs(m(r,r));
    for(SMIterator<scalarD> beg=m.begin(r),end=m.end(r); beg!=end; ++beg)
      *beg/=diag;
  }
}
scalar weightLaplacian(const ObjMesh& mesh,int v0,int v1,const vector<int>& tss)
{
  scalar w=0.0f;
  ASSERT_MSG(tss.size() == 2,"We only handle closed manifold mesh!")
  const Vec3& p0=mesh.getV()[v0];
  const Vec3& p1=mesh.getV()[v1];

  const Vec3i& tri0=mesh.getI()[tss[0]];
  for(char v=0; v<3; v++)
    if(tri0[v] != v0 && tri0[v] != v1) {
      const Vec3& p=mesh.getV()[tri0[v]];
      w+=1.0f/std::tan(getAngle3D<scalar>(p0-p,p1-p));
      break;
    }

  const Vec3i& tri1=mesh.getI()[tss[1]];
  for(char v=0; v<3; v++)
    if(tri1[v] != v0 && tri1[v] != v1) {
      const Vec3& p=mesh.getV()[tri1[v]];
      w+=1.0f/std::tan(getAngle3D<scalar>(p0-p,p1-p));
      break;
    }
  return w;
}
void DHMSmoother::buildK(SMAT& L) const
{
  sizeType nrV=(sizeType)_smesh.getV().size();
  sizeType nrT=(sizeType)_smesh.getI().size();

  //denominator
  vector<scalar> area(nrV,0.0f);
  for(int i=0; i<nrT; i++) {
    const Vec3i& tri=_smesh.getI()[i];
    scalar areaTri=_smesh.getArea(i);
    area[tri[0]]+=areaTri;
    area[tri[1]]+=areaTri;
    area[tri[2]]+=areaTri;
  }

  //numerator
  TRIPS trips;
  ObjMesh::EdgeMap eMap;
  _smesh.buildEdge(eMap);
  Vec3d I=Vec3d::Ones();
  for(map<std::pair<int,int>,ObjMesh::Edge,ObjMesh::EdgeMap::LSS>::const_iterator
      beg=eMap._ess.begin(),end=eMap._ess.end(); beg!=end; beg++) {
    int t0=beg->first.first;
    int t1=beg->first.second;
    scalar w=weightLaplacian(_smesh,t0,t1,beg->second._tris);

    addI(trips,t0*3,t1*3, I*w/(4.0f*area[t0]));
    addI(trips,t0*3,t0*3,-I*w/(4.0f*area[t0]));
    addI(trips,t1*3,t0*3, I*w/(4.0f*area[t1]));
    addI(trips,t1*3,t1*3,-I*w/(4.0f*area[t1]));
  }
  L.resize(nrV*3);
  L.buildFromTripletsDepulicate(trips,0.0f);
  if(_mesh._tree->get<bool>("normalizeLaplacian",true))
    normalizeLaplacian(L);
}
void DHMSmoother::smoothSurface()
{
  //initialize mean curvature flow
  scalar dt=_mesh._tree->get<scalar>("dtLaplacian");
  if(dt == 0.0f)
    return;
  scalar beta0=_smesh.getVolume();
  sizeType nrV=(sizeType)_smesh.getV().size();

  //get pos0
  COL pos0,pos1;
  pos0.resize(nrV*3);
  pos1.resize(nrV*3);
  for(sizeType r=0; r<nrV; r++)
    pos0.block<3,1>(r*3,0)=_smesh.getV()[r].cast<scalarD>();

  //build LHS
  SMAT L;
  buildK(L);
  L.mul(-dt);
  for(sizeType r=0; r<nrV*3; r++)
    L.addToElement(r,r,1.0f);

  //solve for new pos
  _sol.setSolverParameters(1E-6f,10000);
  _sol.setMatrix(L,true);
  _sol.solve(pos0,pos1);
  for(sizeType r=0; r<nrV; r++)
    _smesh.getV()[r]=pos1.block<3,1>(r*3,0).cast<scalar>();
  //preserve volume
  scalar coef=std::pow(beta0/_smesh.getVolume(),(scalar)(1.0f/3.0f));
  for(sizeType r=0; r<nrV; r++)
    _smesh.getV()[r]*=coef;

  //write back
  INFOV("Laplacian smooth solved in %d iters",_sol.getIterationsCount())
  for(boost::unordered_map<sizeType,sizeType>::const_iterator
      beg=_surf.begin(),end=_surf.end(); beg!=end; beg++)
    _mesh.getVert(beg->first)._pos=_smesh.getV()[beg->second];
}
void DHMSmoother::smoothInterior(DHMEnergyPool& E)
{
  //initialize mean curvature flow
  sizeType maxIter=_mesh._tree->get<sizeType>("maxIterSmooth");
  scalar thres=_mesh._tree->get<scalar>("relThresSmooth");

  //build hessian
  SMAT hess;
  hess.resize(_mesh.nrVert()*3,_mesh.nrVert()*3);
  E.buildHessian(hess);
  _sol.setSolverParameters(1E-6f,10000);
  _sol.setMatrix(hess,true);

  //main loop
  COL deriv,delta;
  scalarD deltaMax=0.0f;
  INFO("Smooth Interior Main Loop")
  for(sizeType i=0; i<maxIter; i++) {
    //calculate hessian and derivative
    deriv.setZero(_mesh.nrVert()*3);
    E.buildDeriv(deriv);

    //solve system for delta
    COL pos;
    _mesh.getPos(pos);
    delta.setZero(_mesh.nrVert()*3);
    _sol.solve(deriv,delta);
    pos-=delta;
    _mesh.setPos(pos);

    //check for termination condition
    if(deltaMax == 0.0f && delta.cwiseAbs().maxCoeff() > 0.0f)
      deltaMax=delta.cwiseAbs().maxCoeff();
    INFOV("Delta %f/%f solved in %d iters",delta.cwiseAbs().maxCoeff(),deltaMax,_sol.getIterationsCount())
    if(deltaMax == 0.0f || delta.cwiseAbs().maxCoeff() < thres*deltaMax)
      break;
  }
}
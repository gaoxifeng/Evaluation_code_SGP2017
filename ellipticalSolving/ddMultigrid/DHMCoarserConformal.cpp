#include "DHMCoarserConformal.h"
#include "DHMEnergy.h"
#include "DHMUtil.h"
#include "CommonFile/solvers/PMinres.h"

USE_PRJ_NAMESPACE

DHMCoarserGeomConformal::DHMCoarserGeomConformal(DHMMesh& meshF,DHMMesh& meshC,const SMAT& P)
  :DHMCoarserGeom(meshF,meshC,P)
{
  setParameter();
}
void DHMCoarserGeomConformal::setParameter()
{
  putNoOverwrite<sizeType>(*(_meshF._tree),"maxIterConformal",1000);
  putNoOverwrite<scalar>(*(_meshF._tree),"relThresConformal",1E-3f);
}
bool DHMCoarserGeomConformal::updateWeight(DHMEnergyPool& EC) const
{
  sizeType nrInvalid=0;
  for(sizeType i=0; i<_meshC.nrCell(); i++)
    if(!DHMUniformValidityBound::isValid(_meshC.getCell(i))) {
      scalarD& w=EC.getVDEnergy<DHMConformalEnergy>(i)._weight;
      w=std::max<scalarD>(w,0.01f)*2.0f;
      nrInvalid++;
    }
  INFOV("Found %ld invalid elements!",nrInvalid)
  return nrInvalid>0;
}
FORCE_INLINE Vec3d localTID(const DHMCell& cell,sizeType vid)
{
  for(char v=0; v<8; v++)
    if(cell._verts[v]->_index == vid)
      return Vec3d(v&1?1:-1,v&2?1:-1,v&4?1:-1);
  ASSERT(false);
  return Vec3d::Constant(numeric_limits<scalarD>::max());
}
Mat3d DHMCoarserGeomConformal::findTransfer(const DHMCell& cellF,const DHMCell& cellC) const
{
  char d=0;
  Vec2i pair[4];
  for(char v=0; v<8; v++) {
    sizeType idF=cellF._verts[v]->_index;
    if(_P.numElement(idF) == 1) {
      sizeType idC=_P.getValue()[_P.getRowOffset()[idF]].second;
      pair[d++]=Vec2i(idF,idC);
      break;
    }
  }
  ASSERT(d == 1)

  for(char v=0; v<8; v++) {
    sizeType idF=cellF._verts[v]->_index;
    if(_P.numElement(idF) == 2) {
      sizeType idC0=_P.getValue()[_P.getRowOffset()[idF]+0].second;
      sizeType idC1=_P.getValue()[_P.getRowOffset()[idF]+1].second;
      ASSERT(idC0 == pair[0][1] || idC1 == pair[0][1])
      pair[d++]=Vec2i(idF,idC0 == pair[0][1]?idC1:idC0);
    }
  }
  ASSERT(d == 4)

  Mat3d TF,TC;
  for(char c=0; c<3; c++) {
    TF.col(c)=(localTID(cellF,pair[c+1][0])-localTID(cellF,pair[0][0]));
    TC.col(c)=(localTID(cellC,pair[c+1][1])-localTID(cellC,pair[0][1]));
  }
  TF=TF*TC.inverse();
  //ASSERT((TF*TF.transpose()-Mat3d::Identity()).norm() < 1E-6f)
  return TF;
}
void DHMCoarserGeomConformal::coarsen()
{
  //initialize
  sizeType maxIter=_meshF._tree->get<sizeType>("maxIterConformal");
  scalarD thres=_meshF._tree->get<scalarD>("relThresConformal");

  //build conformal energy with default weight 0.01f
  DHMEnergyPool EC;
  for(sizeType i=0; i<_meshC.nrCell(); i++)
    EC.addVDEnergy(boost::shared_ptr<DHMConformalEnergy>(new DHMConformalEnergy(_meshC.getCell(i),0.0f)));
  //build consistency energy
  boost::unordered_map<sizeType,sizeType> pMap;
  buildParentMap(pMap);
  for(sizeType i=0; i<_meshF.nrCell(); i++) {
    if(_meshF.getPadding() == 1 && _meshF.getCell(i)._layer < 0)
      continue;
    if(pMap.find(i) == pMap.end())
      continue;
    const DHMCell& cellC=_meshC.getCell(pMap.find(i)->second);
    const DHMCell& cellF=_meshF.getCell(i);
    Mat3d DFDC=findTransfer(cellF,cellC);
    EC.addVDEnergy(boost::shared_ptr<DHMGradConsistencyEnergy>(new DHMGradConsistencyEnergy(cellC,cellF,DFDC,1.0f)));
  }

  //main loop
  SMAT hess;
  COL deriv,delta;
  scalarD deltaMax=0.0f;
  PMINRESSolver<scalarD,Kernel<scalarD>,NoPreconSolver<scalarD> > sol;
  sol.setSolverParameters(1E-6f,10000);
  INFO("Conformal Coarsening Loop")
  for(sizeType it=0; it<maxIter; it++) {
    //calculate hessian and derivative
    hess.resize(_meshC.nrVert()*3,_meshC.nrVert()*3);
    deriv.setZero(_meshC.nrVert()*3);
    EC.buildHessian(hess);
    EC.buildDeriv(deriv);

    //solve system for delta
    COL pos;
    _meshC.getPos(pos);
    delta.setZero(_meshC.nrVert()*3);
    sol.setMatrix(hess,true);
    sol.solve(deriv,delta);
    pos-=delta;
    _meshC.setPos(pos);

    //check for termination condition
    if(deltaMax == 0.0f && delta.cwiseAbs().maxCoeff() > 0.0f)
      deltaMax=delta.cwiseAbs().maxCoeff();
    INFOV("Delta %f/%f solved in %ld iters",delta.cwiseAbs().maxCoeff(),deltaMax,sol.getIterationsCount())
    if((deltaMax == 0.0f || delta.cwiseAbs().maxCoeff() < thres*deltaMax) && !updateWeight(EC))
      break;

    //ostringstream oss;
    //oss << "./iter" << it << ".vtk";
    //_meshC.writeVTK(oss.str());
  }
}

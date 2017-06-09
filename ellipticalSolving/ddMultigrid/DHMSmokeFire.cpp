#include "DHMSmokeFire.h"

USE_PRJ_NAMESPACE

//Solve Homogeneous-Neumann-Laplacian Equation
DHMLaplacianSolver::DHMLaplacianSolver(const DHMMesh& mesh,DHMVSField& t):_t(t)
{
  _sys.reset(new DHMLinearSystem(mesh,NULL));
  _sys->buildSystem<DHMScalarConsistencyEnergy>();
  _sys->buildRBTag();
  _kdt=0.0f;
  _M=_sys->_lhs;
}
void DHMLaplacianSolver::advance(scalar dt,scalar kappa)
{
  if(_kdt != dt*kappa) {
    _kdt=dt*kappa;
    _sys->buildSystem<DHMDirichletEnergy>();
    _sys->_lhs.mul(_kdt);
    _sys->_lhs.add(_M);
  }

  if(_kdt == 0.0f)
    return;

  //compute RHS
  Cold RHS(_M.rows());
  Cold X=_t._data.block(0,0,_M.rows(),1).cast<scalarD>();
  _M.multiply(X,RHS);
  _sys->solve(RHS,X);
  _t._data.block(0,0,_M.rows(),1)=X.cast<scalar>();
}
//Solve Fire Equation
const scalar DHMSmokeFire::ABSOLUTE=273.16f;
DHMSmokeFire::DHMSmokeFire(Hierarchy& hier):_hier(hier)
{
  const DHMMesh& fmesh=*(hier._mesh.front());
  _cv.reset(new DHMCVField(fmesh));
  _f.reset(new DHMVSField(fmesh));
  _s.reset(new DHMVSField(fmesh));
  _t.reset(new DHMVSField(fmesh));
  _time=0.0f;

  sizeType lv=0;
  while(lv < (sizeType)hier._prolong.size() && hier._prolong[lv].rows() != hier._mesh[lv]->nrVert())
    lv++;
  _adv.reset(new DHMAdvection(*(hier._mesh[lv])));
  _adv->setAdaptiveMesh(&fmesh,lv);

  //_sol.reset(new DHMPressureSolver(*(hier._mesh.front()),*_cv));
  _sol.reset(new DHMPressureSolverMG(hier,*_cv));
  _sol->precomputeLHS();

  _scatter.reset(new DHMLaplacianSolver(*(hier._mesh.front()),*_t));
  setDefaultParam();
}
void DHMSmokeFire::setDefaultParam(bool reset)
{
  //constitutive model
  putNoOverwrite(_tree,"r",46.0f);
  putNoOverwrite(_tree,"c",3000.0f);
  putNoOverwrite(_tree,"ks",1.0f);
  putNoOverwrite(_tree,"kt",(1700.0f+ABSOLUTE));
  putNoOverwrite(_tree,"kappa",0.01f);
  //environment
  putNoOverwrite(_tree,"Tmax",(1700.0f+ABSOLUTE));
  putNoOverwrite(_tree,"T0",(20.0f+ABSOLUTE));
  putNoOverwrite(_tree,"buoy",0.15f);
  putNoOverwrite(_tree,"damping",0.0f);
  putNoOverwrite(_tree,"gx",0.0f);
  putNoOverwrite(_tree,"gy",0.0f);
  putNoOverwrite(_tree,"gz",-9.81f);
  putNoOverwrite<sizeType>(_tree,"stepper",EULER);
  putNoOverwrite(_tree,"useMacCormack",true);
  putNoOverwrite(_tree,"useAdaptiveBasis",true);

  if(reset) {
    _f->clear(0.0f);
    _s->clear(0.0f);
    _t->clear(_tree.get<scalar>("T0"));
    _cv->clear(0.0f);
  }
}
void DHMSmokeFire::advance(scalar dt)
{
  scalar r=_tree.get<scalar>("r");
  scalar c=_tree.get<scalar>("c");
  scalar ks=_tree.get<scalar>("ks");
  scalar kt=_tree.get<scalar>("kt");
  scalar kappa=_tree.get<scalar>("kappa");
  scalar Tmax=_tree.get<scalar>("Tmax");
  scalar T0=_tree.get<scalar>("T0");
  scalar buoy=_tree.get<scalar>("buoy");
  scalar damping=_tree.get<scalar>("damping");
  Vec3 gravity(_tree.get<scalar>("gx"),_tree.get<scalar>("gy"),_tree.get<scalar>("gz"));
  DHMStepper stepper=(DHMStepper)_tree.get<sizeType>("stepper");
  bool useMacCormack=_tree.get<bool>("useMacCormack");
  bool useAdaptiveBasis=_tree.get<bool>("useAdaptiveBasis");

  //external force:
  //damp, buoyancy, density
  TBEGT
  sizeType nrC=_hier._mesh.front()->nrCell();
  sizeType nrP=_adv->getNrSamplePerCell();
  const vector<DHMAdvection::AdvPoint>& pts=_adv->getSamplePoints();
  DHMAdvection::COL samples(nrC*nrP*3);
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrC; i++) {
    DHMCVField::COMPREF comp=_cv->getCompRef(i);
    comp*=exp(-damping*dt);
    for(sizeType p=0,off=i*nrP; p<nrP; p++,off++) {
      const DHMAdvection::AdvPoint& pt=pts[off];
      samples.block<3,1>(off*3,0)=(buoy*(_t->getScalar<scalar>(0,pt._cid,pt._crd)-T0)*dt)*-gravity;
    }
  }
  DHMCVField force(*(_hier._mesh.front()));
  DHMFieldConverter::convertAdaptive(force,samples,useAdaptiveBasis);
  _cv->_data+=force._data;
  TENDT("EF",_tree)

  //combustion model:
  //fuel damping
  //fuel to soot
  //fuel to temperature
  //temperature damping
  TBEGT
  scalar coef=pow(T0-Tmax,scalar(4.0f));
  sizeType nrVertNC=_hier._mesh.front()->nrVertNC();
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrVertNC; i++) {
    scalar& sVal=_s->_data[i];
    scalar& tVal=_t->_data[i];
    scalar& fVal=_f->_data[i];

    scalar expVal=exp(-r*dt);
    sVal+=ks*fVal*(1.0f-expVal);
    scalar TCubic=3.0f*c*pow(std::max<scalar>(tVal-T0,0.0f),3);
    scalar deltaTCubic=coef*TCubic/(dt*TCubic+coef);
    tVal=pow(deltaTCubic/(3.0f*c),scalar(1.0f/3.0f))+T0;
    tVal+=kt*fVal*(1.0f-expVal);
    fVal*=expVal;
  }
  TENDT("CM",_tree)

  //scatter temperature
  TBEGT
  _scatter->advance(dt,kappa);
  TENDT("ST",_tree)

  //enforce divergence free
  TBEGT
  _sol->project();
  TENDT("PRJ",_tree)

  //advect
  TBEGT
  _adv->setUseAdaptiveBasis(useAdaptiveBasis);
  _adv->_source.resize(3);
  DHMSource::parseSource(_tree,_time,_adv->_source[0],"sourceF");
  DHMSource::parseSource(_tree,_time,_adv->_source[1],"sourceS");
  DHMSource::parseSource(_tree,_time,_adv->_source[2],"sourceT");
  DHMAdvection::CDHMVSFieldPTR vsPtr[3]= {_f.get(),_s.get(),_t.get()};
  DHMAdvection::DHMVSFieldPTR vsOutPtr[3]= {_f.get(),_s.get(),_t.get()};
  if(useMacCormack)
    _adv->advectMacCormack(*_cv,_cv.get(),vsPtr,vsOutPtr,3,dt,stepper);
  else
    _adv->advectSemiLagrangian(*_cv,_cv.get(),vsPtr,vsOutPtr,3,dt,stepper);
  TENDT("ADV",_tree)

  //advance time
  _time+=dt;
  if(_tree.get<bool>("profile",true))
    profileWrite();
}
void DHMSmokeFire::writeFrame(ostream& os) const
{
  writeBinaryData(_cv->_data,os);
  writeBinaryData(_f->_data,os);
  writeBinaryData(_s->_data,os);
  writeBinaryData(_t->_data,os);
  writeBinaryData(_time,os);
}
void DHMSmokeFire::readFrame(istream& is)
{
  readBinaryData(_cv->_data,is);
  readBinaryData(_f->_data,is);
  readBinaryData(_s->_data,is);
  readBinaryData(_t->_data,is);
  readBinaryData(_time,is);
}
void DHMSmokeFire::profileWrite() const
{
  INFO("Frm")
  INFOV("\tEF: %f",_tree.get<scalarD>("EF"))
  INFOV("\tCM: %f",_tree.get<scalarD>("CM"))
  INFOV("\tST: %f",_tree.get<scalarD>("ST"))
  INFOV("\tPRJ: %f",_tree.get<scalarD>("PRJ"))
  INFOV("\tADV: %f",_tree.get<scalarD>("ADV"))
}
const DHMVSField& DHMSmokeFire::getFuel() const
{
  return *_f;
}
const DHMVSField& DHMSmokeFire::getSoot() const
{
  return *_s;
}
const DHMVSField& DHMSmokeFire::getTemp() const
{
  return *_t;
}
const DHMCVField& DHMSmokeFire::getVel() const
{
  return *_cv;
}
DHMVSField& DHMSmokeFire::getFuel()
{
  return *_f;
}
DHMVSField& DHMSmokeFire::getSoot()
{
  return *_s;
}
DHMVSField& DHMSmokeFire::getTemp()
{
  return *_t;
}
DHMCVField& DHMSmokeFire::getVel()
{
  return *_cv;
}
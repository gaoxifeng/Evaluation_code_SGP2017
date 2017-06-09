#include "DHMWater.h"

USE_PRJ_NAMESPACE

//Solve Water Equation
DHMWater::DHMWater(Hierarchy& hier):_hier(hier)
{
  const DHMMesh& fmesh=*(hier._mesh.front());
  _pset.reset(new DHMParticleSet);
  _cv.reset(new DHMCVField(fmesh));
  _lv.reset(new DHMVSField(fmesh));
  _time=0.0f;

  sizeType lv=0;
  while(lv < (sizeType)hier._prolong.size() && hier._prolong[lv].rows() != hier._mesh[lv]->nrVert())
    lv++;
  _adv.reset(new DHMAdvection(*(hier._mesh[lv])));
  _adv->setAdaptiveMesh(&fmesh,lv);
  _pset->_cRad=_adv->getCRad();

  _sol.reset(new DHMPressureSolver(*(hier._mesh.front()),*_cv,_lv.get()));
  //_sol.reset(new DHMPressureSolverMG(hier,*_cv,_lv.get()));
  _sol->precomputeLHS();
  setDefaultParam();
}
void DHMWater::setDefaultParam(bool reset)
{
  putNoOverwrite(_tree,"gx",0.0f);
  putNoOverwrite(_tree,"gy",0.0f);
  putNoOverwrite(_tree,"gz",-9.81f);
  putNoOverwrite(_tree,"FLIPWeight",0.9f);
  putNoOverwrite<sizeType>(_tree,"nrPPerDim",4);
  putNoOverwrite(_tree,"useAdaptiveBasis",true);

  if(reset) {
    _pset->clear();
    _pset->_cRad=_adv->getCRad()/(scalar)_tree.get<sizeType>("nrPPerDim")/2.0f;
    _lv->clear(0.0f);
    _cv->clear(0.0f);
  }
}
void DHMWater::advance(scalar dt)
{
  scalar FLIPWeight=_tree.get<scalar>("FLIPWeight");
  Vec3 gravity(_tree.get<scalar>("gx"),_tree.get<scalar>("gy"),_tree.get<scalar>("gz"));
  bool useAdaptiveBasis=_tree.get<bool>("useAdaptiveBasis");

  //external force:
  TBEGT
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<_pset->size(); i++) {
    ParticleN<scalar>& pt=(*_pset)[i];
    pt._vel+=dt*gravity;
  }
  TENDT("EF",_tree)

  //save old velocity field
  DHMCVField cvOld(*(_hier._mesh.front()));
  _pset->transferVel(*_cv);
  cvOld._data=_cv->_data;

  //enforce divergence free
  TBEGT
  _sol->project();
  TENDT("PRJ",_tree)

  //advect
  TBEGT
  _adv->setUseAdaptiveBasis(useAdaptiveBasis);
  _adv->_source.resize(1);
  DHMSource::parseSource(_tree,_time,_adv->_source[0],"sourceW");
  _adv->advectParticleSet(*_cv,*_lv,*_pset,dt,&cvOld,FLIPWeight);
  TENDT("ADV",_tree)

  //correct
  TBEGT
  TENDT("COR",_tree)

  //adaptive
  TBEGT
  TENDT("ADA",_tree)
}
//IO
DHMParticleSet& DHMWater::getPSet()
{
  return *_pset;
}
const DHMParticleSet& DHMWater::getPSet() const
{
  return *_pset;
}
void DHMWater::writeFrame(ostream& os) const
{
  _pset->write(os);
  writeBinaryData(_cv->_data,os);
  writeBinaryData(_lv->_data,os);
  writeBinaryData(_time,os);
}
void DHMWater::readFrame(istream& is)
{

  _pset->read(is);
  readBinaryData(_cv->_data,is);
  readBinaryData(_lv->_data,is);
  readBinaryData(_time,is);
}
void DHMWater::profileWrite() const
{
  INFO("Frm")
  INFOV("\tEF: %f",_tree.get<scalarD>("EF"))
  INFOV("\tPRJ: %f",_tree.get<scalarD>("PRJ"))
  INFOV("\tADV: %f",_tree.get<scalarD>("ADV"))
  INFOV("\tCOR: %f",_tree.get<scalarD>("COR"))
  INFOV("\tADA: %f",_tree.get<scalarD>("ADA"))
}
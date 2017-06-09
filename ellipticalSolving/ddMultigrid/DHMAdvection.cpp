#include "DHMAdvection.h"
#include "DHMUtil.h"
#include "DHMEnergy.h"
#include "DHMDiscreteGrid.h"
#include <boost/filesystem/operations.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/shared_array.hpp>

USE_PRJ_NAMESPACE

//Source Setup
//parser
void DHMSource::readStamp(const boost::property_tree::ptree& tree,vector<scalar>& stamps)
{
  istringstream iss(tree.get<std::string>("stamps","0"));
  sizeType nr;
  iss >> nr;

  stamps.resize(nr);
  for(sizeType i=0; i<nr; i++)
    iss >> stamps[i];
}
bool DHMSource::isValidSource(const boost::property_tree::ptree& tree,scalar time)
{
  vector<scalar> stamps;
  readStamp(tree,stamps);

  scalar dur=stamps.back()-stamps.front();
  if(time < stamps.front())
    return false;

  scalar frac=fmod(time-stamps.front(),dur)+stamps.front();
  if(!tree.get<bool>("loop"))
    frac=time;
  for(sizeType i=0; i<(sizeType)stamps.size(); i+=2)
    if(frac >= stamps[i] && frac <= stamps[i+1])
      return true;
  return false;
}
void DHMSource::addSourcePoint(const std::vector<DHMSource>& s,DHMParticleSet& pset,scalar cRad)
{
  BBox<scalar> bb;
  ParticleN<scalar> pt;
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<(sizeType)s.size(); i++) {
    pt._vel=s[i]._vel;
    if(s[i]._type == BOX) {
      bb=s[i]._bb;
    } else if(s[i]._type == SPHERE) {
      bb._minC=s[i]._sphere._ctr-Vec3::Constant(s[i]._sphere._rad);
      bb._maxC=s[i]._sphere._ctr+Vec3::Constant(s[i]._sphere._rad);
    }
    for(scalar x=bb._minC[0]; x<=bb._maxC[0]; x+=cRad*2.0f)
      for(scalar y=bb._minC[1]; y<=bb._maxC[1]; y+=cRad*2.0f)
        for(scalar z=bb._minC[2]; z<=bb._maxC[2]; z+=cRad*2.0f)
          if(!pset.hasNeigh(pt._pos=Vec3(x,y,z),cRad*0.95f))
            if(s[i]._type == BOX) {
              OMP_CRITICAL_
              pset.addParticle(pt);
            } else if(s[i]._type == SPHERE && s[i]._sphere.intersect(pt._pos,pt._pos)) {
              OMP_CRITICAL_
              pset.addParticle(pt);
            }
  }
}
bool DHMSource::checkSourcePoint(const Vec3& a,const Vec3& b,const std::vector<DHMSource>& s,scalarD& val)
{
  sizeType nrS=(sizeType)s.size();
  for(sizeType i=0; i<nrS; i++)
    if(s[i]._type == SPHERE && s[i]._sphere.intersect(a,b)) {
      val=s[i]._value;
      return true;
    } else if(s[i]._type == BOX && s[i]._bb.intersect(a,b)) {
      val=s[i]._value;
      return true;
    }
  return false;
}
void DHMSource::parseSingleSource(boost::property_tree::ptree& tree,scalar time,vector<DHMSource>& source,const std::string& name)
{
  if(tree.get<std::string>("name") != name || !isValidSource(tree,time))
    return;
  DHMSource s;
  s._vel=Vec3(tree.get<scalar>("velX",0.0f),tree.get<scalar>("velY",0.0f),tree.get<scalar>("velZ",0.0f));
  s._value=tree.get<scalar>("val");
  s._type=(TYPE)(tree.get<sizeType>("type"));
  Vec3 ctr(tree.get<scalar>("ctrX"),tree.get<scalar>("ctrY"),tree.get<scalar>("ctrZ"));
  if(s._type == SPHERE) {
    s._sphere=Sphere<scalar>(ctr,tree.get<scalar>("rad"));
  } else if(s._type == BOX) {
    Vec3 ext(tree.get<scalar>("extX"),tree.get<scalar>("extY"),tree.get<scalar>("extZ"));
    s._bb=BBox<scalar>(ctr-ext,ctr+ext);
  }
  source.push_back(s);
}
void DHMSource::parseSource(boost::property_tree::ptree& tree,scalar time,vector<DHMSource>& source,const std::string& name)
{
  source.clear();
  for(sizeType i=0;; i++) {
    ostringstream oss;
    oss << "source" << i;
    if(tree.find(oss.str()) != tree.not_found())
      parseSingleSource(tree.get_child(oss.str()),time,source,name);
    else break;
  }
}
void DHMSource::writeSource(const boost::property_tree::ptree& tree)
{
  for(sizeType i=0;; i++) {
    ostringstream oss;
    oss << "source" << i;
    if(tree.find(oss.str()) != tree.not_found()) {
      INFO(oss.str().c_str())
      boost::property_tree::write_xml(cout,tree.get_child(oss.str()),
                                      boost::property_tree::xml_writer_make_settings(' ',2));
    } else break;
  }
}
//creator
sizeType DHMSource::findNewId(boost::property_tree::ptree& tree)
{
  sizeType i=0;
  while(true) {
    ostringstream oss;
    oss << "source" << i;
    if(tree.find(oss.str()) == tree.not_found())
      break;
    i++;
  }
  return i;
}
sizeType DHMSource::addSphereSource(boost::property_tree::ptree& tree,const Vec3& ctr,const scalar rad,scalar val,const std::string& name)
{
  sizeType id=findNewId(tree);
  ostringstream oss;
  oss << "source" << id;
  tree.put<sizeType>(oss.str()+".type",SPHERE);
  tree.put(oss.str()+".ctrX",ctr[0]);
  tree.put(oss.str()+".ctrY",ctr[1]);
  tree.put(oss.str()+".ctrZ",ctr[2]);
  tree.put(oss.str()+".rad",rad);
  tree.put(oss.str()+".val",val);
  tree.put(oss.str()+".name",name);
  addSourceDuration(tree,id,0,1);
  setSourceLoop(tree,id,true);
  return id;
}
sizeType DHMSource::addBoxSource(boost::property_tree::ptree& tree,const BBox<scalar>& bb,scalar val,const std::string& name)
{
  sizeType id=findNewId(tree);
  ostringstream oss;
  oss << "source" << id;
  Vec3 ctr=(bb._maxC+bb._minC)/2;
  Vec3 ext=(bb._maxC-bb._minC)/2;
  tree.put<sizeType>(oss.str()+".type",BOX);
  tree.put(oss.str()+".ctrX",ctr[0]);
  tree.put(oss.str()+".ctrY",ctr[1]);
  tree.put(oss.str()+".ctrZ",ctr[2]);
  tree.put(oss.str()+".extX",ext[0]);
  tree.put(oss.str()+".extY",ext[1]);
  tree.put(oss.str()+".extZ",ext[2]);
  tree.put(oss.str()+".val",val);
  tree.put(oss.str()+".name",name);
  addSourceDuration(tree,id,0,1);
  setSourceLoop(tree,id,true);
  return id;
}
void DHMSource::setSourceParticleVel(boost::property_tree::ptree& tree,sizeType id,const Vec3& vel)
{
  ostringstream oss;
  oss << "source" << id;
  tree.put(oss.str()+".velX",vel[0]);
  tree.put(oss.str()+".velY",vel[1]);
  tree.put(oss.str()+".velZ",vel[2]);
}
//modifier
void DHMSource::addSourceDuration(boost::property_tree::ptree& tree,sizeType i,scalar t0,scalar t1)
{
  ostringstream ossChild;
  ossChild << "source" << i;
  boost::property_tree::ptree& child=tree.get_child(ossChild.str());

  vector<scalar> stamps;
  readStamp(child,stamps);

  vector<scalar> stampsOut;
  for(sizeType i=0; i<(sizeType)stamps.size(); i+=2)
    if(stamps[i+1] < t0 || stamps[i] > t1) {
      stampsOut.push_back(stamps[i]);
      stampsOut.push_back(stamps[i+1]);
    } else {
      t0=min(t0,stamps[i]);
      t1=max(t1,stamps[i+1]);
    }
  stampsOut.insert(std::upper_bound(stampsOut.begin(),stampsOut.end(),t0),t0);
  stampsOut.insert(std::upper_bound(stampsOut.begin(),stampsOut.end(),t1),t1);

  ostringstream oss;
  oss << stampsOut.size() << " ";
  for(sizeType i=0; i<(sizeType)stampsOut.size(); i++)
    oss << stampsOut[i] << " ";
  child.put("stamps",oss.str());
}
void DHMSource::setSourceLoop(boost::property_tree::ptree& tree,sizeType i,bool loop)
{
  ostringstream oss;
  oss << "source" << i << ".loop";
  tree.put(oss.str(),loop);
}
//Semi-Lagrangian back tracer
template <typename CV>
class Advecter
{
public:
  typedef DHMAdvection::AdvPoint PT;
  Advecter(const DHMAdvection& adv,const CV& cv,scalar dt,DHMStepper stepper)
    :_adv(adv),_cv(cv),_dt(dt),_stepper(stepper) {}
  bool advectEuler(const PT& in,PT& out) const {
    out=in;
    return _adv.trace(out._pos=in._pos-_cv.getVel(0,out._cid,out._crd)*_dt,out._cid,out._crd);
  }
  bool advectRK2(const PT& in,PT& out) const {
    out=in;
    bool bd=_adv.trace(out._pos=in._pos-_cv.getVel(0,out._cid,out._crd)*_dt*0.5f,out._cid,out._crd);
    bd=_adv.trace(out._pos=in._pos-_cv.getVel(0,out._cid,out._crd)*_dt,out._cid,out._crd)||bd;
    return bd;
  }
  bool advectRK4(const PT& in,PT& out) const {
    out=in;
    Vec3 k1=-_cv.getVel(0,out._cid,out._crd);
    bool bd=_adv.trace(out._pos=in._pos+k1*_dt*0.5f,out._cid,out._crd);
    Vec3 k2=-_cv.getVel(0,out._cid,out._crd);
    bd=_adv.trace(out._pos=in._pos+k2*_dt*0.5f,out._cid,out._crd)||bd;
    Vec3 k3=-_cv.getVel(0,out._cid,out._crd);
    bd=_adv.trace(out._pos=in._pos+k3*_dt,out._cid,out._crd)||bd;
    Vec3 k4=-_cv.getVel(0,out._cid,out._crd);
    bd=_adv.trace(out._pos=in._pos+(k1+2*k2+2*k3+k4)*_dt/6.0f,out._cid,out._crd)||bd;
    return bd;
  }
  bool advect(const PT& in,PT& out) const {
    switch(_stepper) {
    case EULER:
      return advectEuler(in,out);
    case RK2:
      return advectRK2(in,out);
    case RK4:
      return advectRK4(in,out);
    }
    return false;
  }
protected:
  //data
  const DHMAdvection& _adv;
  const CV& _cv;
  DHMStepper _stepper;
  scalar _dt;
};
//close point solver
struct DistCallback {
public:
  DistCallback(const Vec3& pos,const DHMMesh& mesh,sizeType& cid,Vec3& crd)
    :_pos(pos),_mesh(mesh),_cid(cid),_crd(crd) {}
  void updateDist(sizeType cid,scalar& sqrDist) {
    Vec3d crd=Vec3d::Zero();
    scalar newSqrDist=(scalar)solveSQP(_mesh.getCell(cid),crd);
    if(newSqrDist < sqrDist) {
      sqrDist=newSqrDist;
      _cid=cid;
      _crd=crd.cast<scalar>();
    }
  }
  scalarD solveSQP(const DHMCell& cell,Vec3d& crd) {
    _vpos[0]=cell._verts[0]->_pos.cast<scalarD>();
    _vpos[1]=cell._verts[1]->_pos.cast<scalarD>();
    _vpos[2]=cell._verts[2]->_pos.cast<scalarD>();
    _vpos[3]=cell._verts[3]->_pos.cast<scalarD>();
    _vpos[4]=cell._verts[4]->_pos.cast<scalarD>();
    _vpos[5]=cell._verts[5]->_pos.cast<scalarD>();
    _vpos[6]=cell._verts[6]->_pos.cast<scalarD>();
    _vpos[7]=cell._verts[7]->_pos.cast<scalarD>();

    _Jx[0]=(_vpos[1]-_vpos[0])/2;
    _Jx[1]=(_vpos[3]-_vpos[2])/2;
    _Jx[2]=(_vpos[5]-_vpos[4])/2;
    _Jx[3]=(_vpos[7]-_vpos[6])/2;

    _Jy[0]=(_vpos[2]-_vpos[0])/2;
    _Jy[1]=(_vpos[3]-_vpos[1])/2;
    _Jy[2]=(_vpos[6]-_vpos[4])/2;
    _Jy[3]=(_vpos[7]-_vpos[5])/2;

    _Jz[0]=(_vpos[4]-_vpos[0])/2;
    _Jz[1]=(_vpos[5]-_vpos[1])/2;
    _Jz[2]=(_vpos[6]-_vpos[2])/2;
    _Jz[3]=(_vpos[7]-_vpos[3])/2;

    //main loop of SQP algorithm
    Mat3d A,A2;
    Vec3d b,b2,delta,dCrd,lambda;
    scalarD alpha,thres=1E-5f;
    bool isActive[3]= {true,true,true};
    for(sizeType it=0; it<1000; it++) {
      //newton step
      lambda=(crd+Vec3d::Ones())/2;
      delta=_pos.cast<scalarD>()-
            interp3D(_vpos[0],_vpos[1],_vpos[2],_vpos[3],
                     _vpos[4],_vpos[5],_vpos[6],_vpos[7],
                     lambda[0],lambda[1],lambda[2]);
      A.col(0)=interp2D(_Jx[0],_Jx[1],_Jx[2],_Jx[3],lambda[1],lambda[2]);
      A.col(1)=interp2D(_Jy[0],_Jy[1],_Jy[2],_Jy[3],lambda[0],lambda[2]);
      A.col(2)=interp2D(_Jz[0],_Jz[1],_Jz[2],_Jz[3],lambda[0],lambda[1]);
      b2=b=A.transpose()*delta;
      A2=A=A.transpose()*A;

      //solve
      for(char i=0; i<3; i++)
        if(!isActive[i]) {
          A2.row(i).setZero();
          A2.col(i).setZero();
          A2(i,i)=1.0f;
          b2[i]=0.0f;
        }

      dCrd=A2.llt().solve(b2);
      lambda=b-A*dCrd;
      bool changed=false;
      for(char i=0; i<3; i++)
        if(!isActive[i] && crd[i]*lambda[i] < 0.0f)
          changed=isActive[i]=true;
      if(!changed && dCrd.norm() < thres)
        return delta.squaredNorm();

      //simple line search
      char minId=lineSearch(crd,dCrd,alpha);
      if(minId >= 0)
        isActive[minId]=false;
      for(char i=0; i<3; i++)
        crd[i]=std::max<scalarD>(std::min<scalarD>(crd[i]+dCrd[i]*alpha,1.0f),-1.0f);
    }
    return numeric_limits<scalar>::max();
  }
  static char lineSearch(const Vec3d& crd,const Vec3d& f,scalarD& alpha) {
    char minId=-1;
    scalarD alphaD;
    alpha=1.0f;
    for(char i=0; i<3; i++)
      if(f[i] < 0.0f) {
        alphaD=(-1.0f-crd[i])/f[i];
        if(alphaD < alpha) {
          alpha=alphaD;
          minId=i;
        }
      } else if(f[i] > 0.0f) {
        alphaD=(1.0f-crd[i])/f[i];
        if(alphaD < alpha) {
          alpha=alphaD;
          minId=i;
        }
      }
    return minId;
  }
  //data
  const Vec3& _pos;
  const DHMMesh& _mesh;
  Vec3d _vpos[8],_Jx[4],_Jy[4],_Jz[4];
  sizeType& _cid;
  Vec3& _crd;
};
//advector
DHMAdvection::DHMAdvection(const DHMMesh& mesh):_mesh(mesh)
{
  //build spatial hierarchy: leaf
  _cells.reset(new DHMDiscreteGrid(mesh));
  setAdaptiveMesh(&mesh,0);
  _useAdaptiveBasis=true;
}
void DHMAdvection::setAdaptiveMesh(const DHMMesh* meshA,sizeType lv)
{
  if(!meshA) {
    setAdaptiveMesh(&_mesh,0);
    return;
  }

  _sys.reset(new DHMLinearSystem(*meshA,NULL));
  _sys->buildSystem<DHMScalarConsistencyEnergy>();
  _sys->lumpSystem();
  _dim=sizeType(1)<<lv;

  //compute pts
  sizeType nrC=meshA->nrCell();
  sizeType nrP=GaussLegendreIntegral<scalarD,scalarD>::nrP3D(2);
  _pts.resize(nrP*nrC);
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& cell=meshA->getCell(i);
    for(sizeType p=0; p<nrP; p++) {
      AdvPoint& ptsI=_pts[i*nrP+p];
      ptsI._cid=i;
      ptsI._crd=GaussLegendreIntegral<scalarD,scalarD>::point3D(2,p).cast<scalar>();
      ptsI._pos=DHMMapping::M(cell,ptsI._crd);
    }
  }
}
void DHMAdvection::setUseAdaptiveBasis(bool use)
{
  _useAdaptiveBasis=use;
}
void DHMAdvection::advectParticleSet(const DHMCVField& cv,const DHMVSField& lv,DHMParticleSet& pset,scalar dt,const DHMCVField* cvOld,scalar FLIPWeight) const
{
  Vec3 crd,cvv;
  sizeType cid,nrP=pset.size();
  OMP_PARALLEL_FOR_I(OMP_PRI(crd,cvv,cid))
  for(sizeType i=0; i<nrP; i++) {
    //advection step
    ParticleN<scalar>& pt=pset[i];
    trace(pt._pos,cid,crd);
    if(lv.getScalar<scalar>(0,cid,crd) > 0.0f) {
      pt._pos+=pt._vel*dt;	//revert to EULER
    } else {
      //particle is in fluid body
      //if cvOld/FLIPWeight provided, we adjust velocity
      cvv=cv.getVel(0,cid,crd);
      if(cvOld)
        pt._vel=cvv*(1.0f-FLIPWeight)+(pt._vel+cvv-cvOld->getVel(0,cid,crd))*FLIPWeight;

      //advect along cv
      trace(pt._pos-cvv*dt*0.5f,cid,crd);
      if(lv.getScalar<scalar>(0,cid,crd) > 0.0f)
        pt._pos+=pt._vel*dt;	//revert to EULER
      else pt._pos+=cv.getVel(0,cid,crd)*dt;
    }

    //position correction
    trace(pt._pos,cid,crd);
    pt._pos=DHMMapping::M(_sys->_mesh.getCell(cid),crd);
  }
  //add source
  scalar rd=pset._cRad/(scalar)_dim;
  pset.buildSpatialHash(rd);
  DHMSource::addSourcePoint(_source[0],pset,rd);
}
void DHMAdvection::advectMacCormack(const DHMCVField& cv,DHMCVFieldPTR cvOut,CDHMVSFieldPTR* vs,DHMVSFieldPTR* vsOut,sizeType nr,scalar dt,DHMStepper stepper) const
{
  //Step One:
  Matd vsMINMAX;
  advectSemiLagrangian(cv,cvOut,vs,vsOut,nr,dt,stepper,&vsMINMAX);

  //Step Two:
  vector<boost::shared_ptr<DHMVSField> > vsOutTmp(nr);
  boost::shared_array<DHMVSFieldPTR> vsOutTmpPtr(new DHMVSFieldPTR[nr]);
  for(sizeType i=0; i<nr; i++) {
    vsOutTmp[i].reset(new DHMVSField(_sys->_mesh));
    vsOutTmpPtr[i]=vsOutTmp[i].get();
  }
  advectSemiLagrangian(cv,NULL,const_cast<CDHMVSFieldPTR*>(vsOut),vsOutTmpPtr.get(),nr,-dt,stepper);

  //Step Three: correct and limiting
  sizeType rows=_sys->_lhs.rows();
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<rows; i++) {
    for(sizeType c=0; c<nr; c++)
      //not on boundary and not close to source
      if(!isInf(vsMINMAX(i,c))) {
        scalar& val=vsOut[c]->getCompRef(i)[0];
        scalar e=(vsOutTmpPtr[c]->getComp(i)[0]-vs[c]->getComp(i)[0])*0.5f;
        val-=e;
        if(val < vsMINMAX(i,c) || val > vsMINMAX(i,c+nr))
          val+=e;
      }
  }
}
void DHMAdvection::advectSemiLagrangian(const DHMCVField& cv,DHMCVFieldPTR cvOut,CDHMVSFieldPTR* vs,DHMVSFieldPTR* vsOut,sizeType nr,scalar dt,DHMStepper stepper,Matd* vsMINMAX) const
{
  sizeType nrC=_sys->_mesh.nrCell();
  sizeType nrVertNC=_sys->_mesh.nrVertNC();
  sizeType nrP=GaussLegendreIntegral<scalarD,scalarD>::nrP3D(2);

  vector<Vec3d,Eigen::aligned_allocator<Vec3d> > pts(nrP);
  vector<scalarD> weights(nrP);
  for(sizeType i=0; i<nrP; i++) {
    pts[i]=GaussLegendreIntegral<scalarD,scalarD>::point3D(2,i);
    weights[i]=GaussLegendreIntegral<scalarD,scalarD>::weight3D(2,i);
  }

  COL cvs(nrC*nrP*3);
  Matd vsRHS;
  vsRHS.setZero(nrVertNC,nr);
  if(vsMINMAX) {
    vsMINMAX->resize(nrVertNC,nr*2);
    vsMINMAX->block(0,0,nrVertNC,nr).setConstant(numeric_limits<scalarD>::max());
    vsMINMAX->block(0,nr,nrVertNC,nr).setConstant(-numeric_limits<scalarD>::max());
  }

  //find RHS
  AdvPoint pt;
  CELLID vids;
  DHMScalarConsistencyEnergy::HESS cons;
  DHMScalarConsistencyEnergy::GRAD stencil;

  Eigen::Matrix<scalarD,8,-1> grad(8,nr);
  Rowd MIND(nr),MAXD(nr);

  ASSERT((sizeType)_source.size() >= nr)
  Advecter<DHMCVField> adv(*this,cv,dt,stepper);
  OMP_PARALLEL_FOR_I(OMP_PRI(pt,vids,cons,stencil) OMP_FPRI(grad,MIND,MAXD))
  for(sizeType i=0; i<nrC; i++) {
    grad.setZero();
    const DHMCell& cell=_sys->_mesh.getCell(i);
    DHMScalarConsistencyEnergy e(cell,NULL,NULL);
    MIND.setConstant(numeric_limits<scalarD>::max());
    MAXD.setConstant(-numeric_limits<scalarD>::max());
    for(sizeType p=0,off=i*nrP; p<nrP; p++,off++) {
      bool bd=adv.advect(_pts[off],pt);
      if(cvOut)
        cvs.block<3,1>(off*3,0)=cv.getVel(0,pt._cid,pt._crd);
      stencil=e.gradStencil(pts[p])*weights[p];
      for(sizeType c=0; c<nr; c++) {
        scalarD val=vs[c]->getScalar<scalarD>(0,pt._cid,pt._crd,&(MIND[c]),&(MAXD[c]));
        bool src=DHMSource::checkSourcePoint(_pts[off]._pos,pt._pos,_source[c],val);
        MIND[c]=std::min(MIND[c],(bd||src)?-numeric_limits<scalarD>::infinity():val);
        MAXD[c]=std::max(MAXD[c],(bd||src)?numeric_limits<scalarD>::infinity():val);
        grad.col(c)-=stencil*val;
      }
    }
    if(DHMLinearSystem::findVid(vids,cons,cell))
      grad=cons.transpose()*grad;
    OMP_CRITICAL_
    for(char v=0; v<8; v++) {
      vsRHS.row(vids[v])-=grad.row(v);
      if(vsMINMAX)
        for(sizeType c=0; c<nr; c++) {
          (*vsMINMAX)(vids[v],c)=min((*vsMINMAX)(vids[v],c),MIND[c]);
          (*vsMINMAX)(vids[v],nr+c)=max((*vsMINMAX)(vids[v],nr+c),MAXD[c]);
        }
    }
  }

  if(cvOut)
    DHMFieldConverter::convertAdaptive(*cvOut,cvs,_useAdaptiveBasis);
  for(sizeType r=0; r<_sys->_lhs.rows(); r++) {
    scalarD diag=_sys->_lhs(r,r);
    for(sizeType c=0; c<nr; c++)
      vsOut[c]->getCompRef(r)[0]=(scalar)(vsRHS(r,c)/diag);
  }
}
bool DHMAdvection::trace(const Vec3& pos,sizeType& cid,Vec3& crd) const
{
  scalar sqrDist,thres=_cells->eps0()*1E-3f;
  DistCallback cb(pos,_mesh,cid,crd);
  _cells->pointQuery(pos,thres*thres,cb,sqrDist);
  if(!_sys->_mesh.getAdaptive().empty())
    cid=_sys->_mesh.getAdaptive()[cid];
  if(cid < 0) {
    //this means the cell is subdivided
    //we need to find the actual subcell
    crd=(crd+Vec3::Ones())*(scalar)_dim/2.0f;
    Vec3i cOff=compMin(floor(crd),Vec3i::Constant(_dim-1));
    cid=(-1-cid)+Vec3i(_dim*_dim,_dim,1).dot(cOff);
    crd=(crd-cOff.cast<scalar>())*2.0f-Vec3::Ones();
  }
  return sqrDist > thres*thres;
}
const vector<DHMAdvection::AdvPoint>& DHMAdvection::getSamplePoints() const
{
  return _pts;
}
sizeType DHMAdvection::getNrSamplePerCell() const
{
  return GaussLegendreIntegral<scalarD,scalarD>::nrP3D(2);
}
scalar DHMAdvection::getCRad() const
{
  scalar avgEdge=0.0f;
  sizeType nrC=_mesh.nrCell();
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& c=_mesh.getCell(i);
    for(char e=0; e<12; e++) {
      Vec2i vid=c.getEdge(e);
      avgEdge+=(_mesh.getVert(vid[0])._pos-_mesh.getVert(vid[1])._pos).norm();
    }
  }
  avgEdge/=(scalar)nrC*12.0f;
  return avgEdge;
}
void DHMAdvection::debug() const
{
  boost::filesystem::create_directory("./DebugAdv");

  //build bounding box
  BBox<scalar> bb;
  for(sizeType i=0; i<_mesh.nrVert(); i++)
    bb.setUnion(_mesh.getVert(i)._pos);
  bb.enlargedEps(1.0f);

  //build line
  Vec3 crd;
  sizeType cid;
  Vec3 ctr=(bb._minC+bb._maxC)/2.0f;
  scalar rad=bb.getExtent().norm()/2.0f;
  for(sizeType i=0; i<1000; i++) {
    Vec3 pos=ctr+Vec3::Random()*rad;
    bool interior=trace(pos,cid,crd);
    Vec3 posF=DHMMapping::M(_sys->_mesh.getCell(cid),crd);

    //write
    POSES lines;
    lines.push_back(pos);
    lines.push_back(posF);

    vector<scalar> css;
    css.push_back(0.0f);
    css.push_back(1.0f);

    ostringstream oss;
    oss << "./DebugAdv/frm" << i << ".vtk";

    VTKWriter<scalar> os("Advection",oss.str(),true);
    os.appendPoints(lines.begin(),lines.end());
    os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                   VTKWriter<scalar>::IteratorIndex<Vec3i>(1,2,0),
                   VTKWriter<scalar>::LINE);
    os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
                   VTKWriter<scalar>::IteratorIndex<Vec3i>(2,0,0),
                   VTKWriter<scalar>::POINT);
    os.appendCustomPointData("Color",css.begin(),css.end());
  }

  //debug inside/outside tag
  POSES lines;
  vector<scalar> css;
  for(scalar x=bb._minC[0]; x<bb._maxC[0]; x+=_cells->eps0())
    for(scalar y=bb._minC[1]; y<bb._maxC[1]; y+=_cells->eps0())
      for(scalar z=bb._minC[2]; z<bb._maxC[2]; z+=_cells->eps0()) {
        lines.push_back(Vec3(x,y,z));
        css.push_back(trace(Vec3(x,y,z),cid,crd) ? 1.0f : 0.0f);
      }
  VTKWriter<scalar> os("Advection","./DebugAdv/InOut.vtk",true);
  os.appendPoints(lines.begin(),lines.end());
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)lines.size(),0,0),
                 VTKWriter<scalar>::POINT);
  os.appendCustomPointData("Color",css.begin(),css.end());
}
void DHMAdvection::debug(const std::string& path,const Matd& vsMINMAX) const
{
  sizeType rows=vsMINMAX.rows();
  sizeType nr=vsMINMAX.cols()/2;
  DHMVSField diff(_sys->_mesh);
  DHMVSField bdTag(_sys->_mesh);
  for(sizeType i=0; i<nr; i++) {
    for(sizeType r=0; r<rows; r++) {
      diff.getCompRef(r)[0]=(scalar)(vsMINMAX(r,i+nr)-vsMINMAX(r,i));
      if(isInf(diff.getCompRef(r)[0])) {
        diff.getCompRef(r)[0]=0.0f;
        bdTag.getCompRef(r)[0]=1.0f;
      } else bdTag.getCompRef(r)[0]=0.0f;
    }

    ostringstream oss;
    oss << path << "diff" << (char)('A'+i) << ".vtk";
    diff.writeSVTK(oss.str());

    ostringstream oss2;
    oss2 << path << "bdTag" << (char)('A'+i) << ".vtk";
    bdTag.writeSVTK(oss2.str());
  }
}
//DHMScalarConsistencyEnergy
DHMScalarConsistencyEnergy::DHMScalarConsistencyEnergy(const DHMCell& cell,const ImplicitFunc<scalar>* adv,const DHMVSField* vsOut)
  :DHMEnergy<8>(cell),_adv(adv),_vsOut(vsOut)
{
  setDegree(10);
}
DHMScalarConsistencyEnergy::GRAD DHMScalarConsistencyEnergy::grad(const Vec3d& pos) const
{
  Eigen::Matrix<scalarD,8,1> Ms;
  DHMMapping::MStencil(Ms,pos);
  scalarD s=0.0f;
  Vec3 posE=Vec3::Zero();
  for(char v=0; v<8; v++) {
    const DHMVertex& vert=*(_cell._verts[v]);
    s+=_vsOut->getVScalar<scalarD>(0,vert)*Ms[v];
    posE+=vert._pos*(scalar)Ms[v];
  }
  return gradStencil(pos)*(s-(*_adv)(posE));
}
DHMScalarConsistencyEnergy::GRAD DHMScalarConsistencyEnergy::gradStencil(const Vec3d& pos) const
{
  Mat3d JF;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::J(Js,_cell,JF);

  Eigen::Matrix<scalarD,8,1> Ms;
  DHMMapping::MStencil(Ms,pos);
  return Ms.transpose()*JF.determinant();
}
DHMScalarConsistencyEnergy::HESS DHMScalarConsistencyEnergy::hess(const Vec3d& pos) const
{
  Mat3d JF;
  Eigen::Matrix<scalarD,3,8> Js;
  DHMMapping::JStencil(Js,pos);
  DHMMapping::J(Js,_cell,JF);

  Eigen::Matrix<scalarD,8,1> Ms;
  DHMMapping::MStencil(Ms,pos);
  return Ms*Ms.transpose()*JF.determinant();
}
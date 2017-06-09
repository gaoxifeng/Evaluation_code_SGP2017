#include "DHMField.h"
#include "DHMEnergy.h"
#include "DHMLinearSolver.h"

USE_PRJ_NAMESPACE

//velocity consistency energy
class DHMVelConsistencyEnergy : public DHMEnergy<8>
{
public:
  using typename DHMEnergy<8>::COL;
  using typename DHMEnergy<8>::HESS;
  using typename DHMEnergy<8>::GRAD;
  DHMVelConsistencyEnergy(const DHMCell& cell,const VelFunc<scalar>* vel,const GRAD& comp)
    :DHMEnergy(cell),_vel(vel),_comp(comp) {
    setDegree(10);
  }
  GRAD grad(const Vec3d& pos) const {
    Mat3d JF;
    Eigen::Matrix<scalarD,3,8> Js;
    DHMMapping::JStencil(Js,pos);
    DHMMapping::J(Js,_cell,JF);
    Js=JF.inverse().eval().transpose()*Js;

    Vec3 x0=DHMMapping::M(_cell,pos.cast<scalar>());
    Vec3d vel=(*_vel)(x0).cast<scalarD>();
    return (Js.transpose()*(Js*_comp-vel))*JF.determinant();
  }
  template <typename SAMPLES>
  GRAD gradIntFast(const vector<Vec3d,Eigen::aligned_allocator<Vec3d> >& pts,const vector<scalarD>& weights,sizeType off,const SAMPLES& samples) const {
    Mat3d JF;
    Eigen::Matrix<scalarD,3,8> Js;
    sizeType nrP=(sizeType)pts.size();
    GRAD ret=GRAD::Zero();
    for(sizeType i=0; i<nrP; i++,off+=3) {
      DHMMapping::JStencil(Js,pts[i]);
      DHMMapping::J(Js,_cell,JF);
      Js=JF.inverse().eval().transpose()*Js;
      ret-=Js.transpose()*samples.template segment<3>(off).template cast<scalarD>()*(JF.determinant()*weights[i]);
    }
    return ret;
  }
  HESS hess(const Vec3d& pos) const {
    Mat3d JF;
    Eigen::Matrix<scalarD,3,8> Js;
    DHMMapping::JStencil(Js,pos);
    DHMMapping::J(Js,_cell,JF);
    Js=JF.inverse().eval().transpose()*Js;
    return (Js.transpose()*Js)*JF.determinant();
  }
  const VelFunc<scalar>* _vel;
  const GRAD& _comp;
};
bool findVid(DHMTraits::CELLID& vids,DHMLinearSystem::CONS& cons,const DHMCell& cell,const boost::unordered_map<sizeType,sizeType>& cvmap,const sizeType* rowmap)
{
  boost::shared_ptr<DHMConstrainedVertex> cv;
  cons.setZero();
  vids.setConstant(numeric_limits<sizeType>::max());
  bool isCon=false;
  for(char v=0; v<8; v++)
    if(cv=boost::dynamic_pointer_cast<DHMConstrainedVertex>(cell._verts[v])) {
      isCon=true;
      for(char c=0; c<4; c++)
        if(cv->_verts[c]) {
          sizeType cid=cvmap.find(cv->_verts[c]->_index)->second;
          DHMLinearSystem::insertVid(cid,vids,v,cons,cv->_weights[c]);
        }
    } else DHMLinearSystem::insertVid(rowmap[v],vids,v,cons,1.0f);
  for(char v=0; v<8; v++)
    ASSERT(vids[v] != numeric_limits<sizeType>::max());
  return isCon;
}
void DHMFieldConverter::convert(DHMVVField& field,const DHMCVField& func)
{
  sizeType nrC=func._mesh.nrCell();
  sizeType nrV=field.nrComp();

  field.clear(0.0f);
  COL weight=COL::Zero(nrV);
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& cell=field._mesh.getCell(i);
    boost::shared_ptr<DHMConstrainedVertex> cv;
    for(char v=0; v<8; v++) {
      Vec3 vel=func.getVel(0,i,Vec3(v&1?-1.0f:1.0f,v&2?-1.0f:1.0f,v&4?-1.0f:1.0f));
      if(cv=boost::dynamic_pointer_cast<DHMConstrainedVertex>(cell._verts[v])) {
        for(char cvi=0; cvi<4; cvi++)
          OMP_CRITICAL_
          if(cv->_verts[cvi]) {
            sizeType vid=cv->_verts[cvi]->_index;
            field.getCompRef(vid)+=vel*cv->_weights[cvi];
            weight[vid]+=cv->_weights[cvi];
          }
      } else {
        OMP_CRITICAL_ {
          sizeType vid=cell._verts[v]->_index;
          field.getCompRef(vid)+=vel;
          weight[vid]+=1.0f;
        }
      }
    }
  }

  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrV; i++)
    field.getCompRef(i)/=weight[i];
}
void DHMFieldConverter::convert(DHMVVField& field,const VelFunc<scalar>& func)
{
  const DHMMesh& mesh=field._mesh;
  sizeType nrVertNC=mesh.nrVertNC();
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrVertNC; i++)
    field.getCompRef(i)=func(mesh.getVert(i)._pos).cast<scalar>();
}
void DHMFieldConverter::convert(DHMVSField& field,const ImplicitFunc<scalar>& func)
{
  const DHMMesh& mesh=field._mesh;
  sizeType nrVertNC=mesh.nrVertNC();
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrVertNC; i++)
    field.getCompRef(i)[0]=func(mesh.getVert(i)._pos);
}
void DHMFieldConverter::convert(DHMCSField& field,const ImplicitFunc<scalar>& func)
{
  const DHMMesh& mesh=field._mesh;
  sizeType nrC=mesh.nrCell();
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrC; i++) {
    Vec3 pos=Vec3::Zero();
    for(char v=0; v<8; v++)
      pos+=mesh.getCell(i)._verts[v]->_pos;
    field.getCompRef(i)[0]=func(pos*0.125f);
  }
}
void DHMFieldConverter::convertAdaptive(DHMCVField& field,const COL& samples,bool useAdaptiveBasis)
{
  typedef multi_unordered_map<sizeType,sizeType,true> GMAP;
  typedef GMAP::const_iterator CGITER;
  GMAP groups;
  const DHMMesh& mesh=field._mesh;
  sizeType nrC=mesh.nrCell();
  field.clear(0.0f);

  //precompute gauss points
  sizeType nrP=GaussLegendreIntegral<scalarD,scalarD>::nrP3D(2);
  vector<Vec3d,Eigen::aligned_allocator<Vec3d> > pts(nrP);
  vector<scalarD> weights(nrP);
  for(sizeType i=0; i<nrP; i++) {
    pts[i]=GaussLegendreIntegral<scalarD,scalarD>::point3D(2,i);
    weights[i]=GaussLegendreIntegral<scalarD,scalarD>::weight3D(2,i);
  }

  //solve vectors for unconstrained cells
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& cell=mesh.getCell(i);
    bool cons=false;
    if(useAdaptiveBasis)
      for(char v=0; v<8; v++)
        if(boost::dynamic_pointer_cast<DHMConstrainedVertex>(cell._verts[v])) {
          //make use of the fact that index indicate BaseMesh ID
          OMP_CRITICAL_
          groups.insert(cell._index,i);
          cons=true;
          break;
        }
    if(!cons) {
      DHMVelConsistencyEnergy e(cell,NULL,DHMVelConsistencyEnergy::GRAD::Zero());
      field.getCompRef(i)=e.hessInt().ldlt().solve(-e.gradIntFast(pts,weights,i*nrP*3,samples)).eval().cast<scalar>();
    }
  }

  //solve groups
  Matd LHS;
  Cold RHS,PRE;
  GMAP::SET gcids;
  boost::unordered_map<sizeType,sizeType> cvmap;
  vector<sizeType> rowmap;

  CELLIDS vids;
  vector<DHMLinearSystem::CONS> cons;
  DHMVelConsistencyEnergy::HESS hess;
  DHMVelConsistencyEnergy::GRAD grad;

  vector<sizeType> keys;
  for(CGITER beg=groups.begin(),end=groups.end(); beg!=end; beg++)
    keys.push_back(beg->first);
  sizeType nrKey=keys.size();
  OMP_PARALLEL_FOR_I(OMP_PRI(LHS,RHS,PRE,gcids,cvmap,rowmap, vids,cons,hess,grad))
  for(sizeType k=0; k<nrKey; k++) {
    groups.get(keys[k],gcids);
    sizeType nrGC=(sizeType)gcids.size(),nrRow=0;

    //assemble indices
    cvmap.clear();
    rowmap.assign(nrGC*8,-1);
    for(sizeType i=0; i<nrGC; i++) {
      const DHMCell& cell=mesh.getCell(gcids[i]);
      for(char v=0; v<8; v++)
        if(!boost::dynamic_pointer_cast<DHMConstrainedVertex>(cell._verts[v]))
          cvmap[cell._verts[v]->_index]=rowmap[i*8+v]=nrRow++;
    }

    //build hess/grad
    RHS.setZero(nrRow);
    LHS.setZero(nrRow,nrRow);
    vids.resize(nrGC);
    cons.resize(nrGC);
    for(sizeType i=0; i<nrGC; i++) {
      const DHMCell& cell=mesh.getCell(gcids[i]);
      DHMVelConsistencyEnergy e(cell,NULL,DHMVelConsistencyEnergy::GRAD::Zero());
      findVid(vids[i],cons[i],cell,cvmap,&(rowmap[i*8]));
      grad=cons[i].transpose()*e.gradIntFast(pts,weights,gcids[i]*nrP*3,samples);
      hess=cons[i].transpose()*(e.hessInt()*cons[i]).eval();
      for(char r=0; r<8; r++) {
        RHS[vids[i][r]]-=grad[r];
        for(char c=0; c<8; c++)
          LHS(vids[i][r],vids[i][c])+=hess(r,c);
      }
    }

    //solve
    PRE=LHS.ldlt().solve(RHS);
    for(sizeType i=0; i<nrGC; i++) {
      for(char v=0; v<8; v++)
        grad[v]=PRE[vids[i][v]];
      field.getCompRef(gcids[i])=(cons[i]*grad).cast<scalar>();
    }
  }
}
void DHMFieldConverter::convertAdaptive(DHMCVField& field,const VelFunc<scalar>& func,bool useAdaptiveBasis)
{
  const DHMMesh& mesh=field._mesh;
  sizeType nrC=mesh.nrCell();
  sizeType nrP=GaussLegendreIntegral<scalarD,scalarD>::nrP3D(2);

  COL samples(nrC*nrP*3);
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& cell=mesh.getCell(i);
    for(sizeType p=0; p<nrP; p++) {
      Vec3 pt=GaussLegendreIntegral<scalarD,scalarD>::point3D(2,p).cast<scalar>();
      samples.block<3,1>((i*nrP+p)*3,0)=func(DHMMapping::M(cell,pt));
    }
  }
  convertAdaptive(field,samples,useAdaptiveBasis);
}
void DHMFieldConverter::convertAdaptiveDirect(DHMCVField& field,const COL& samples)
{
  typedef multi_unordered_map<sizeType,sizeType,true> GMAP;
  typedef GMAP::const_iterator CGITER;
  GMAP groups;
  const DHMMesh& mesh=field._mesh;
  sizeType nrC=mesh.nrCell();
  field.clear(0.0f);

  //precompute gauss points
  sizeType nrP=GaussLegendreIntegral<scalarD,scalarD>::nrP3D(2);
  vector<Vec3d,Eigen::aligned_allocator<Vec3d> > pts(nrP);
  vector<scalarD> weights(nrP);
  for(sizeType i=0; i<nrP; i++) {
    pts[i]=GaussLegendreIntegral<scalarD,scalarD>::point3D(2,i);
    weights[i]=GaussLegendreIntegral<scalarD,scalarD>::weight3D(2,i);
  }

  //solve vectors for unconstrained cells
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& cell=mesh.getCell(i);
    bool cons=false;
    for(char v=0; v<8; v++)
      if(boost::dynamic_pointer_cast<DHMConstrainedVertex>(cell._verts[v])) {
        //make use of the fact that index indicate BaseMesh ID
        OMP_CRITICAL_
        groups.insert(cell._index,i);
        cons=true;
        break;
      }
    if(!cons) {
      DHMVelConsistencyEnergy e(cell,NULL,DHMVelConsistencyEnergy::GRAD::Zero());
      field.getCompRef(i)=e.hessInt().ldlt().solve(-e.gradIntFast(pts,weights,i*nrP*3,samples)).eval().cast<scalar>();
    }
  }

  //solve groups
  Matd LHS;
  Cold RHS,PRE;
  GMAP::SET gcids;
  boost::unordered_map<sizeType,sizeType> cvmap;

  vector<sizeType> keys;
  for(CGITER beg=groups.begin(),end=groups.end(); beg!=end; beg++)
    keys.push_back(beg->first);
  sizeType nrKey=keys.size();
  OMP_PARALLEL_FOR_I(OMP_PRI(LHS,RHS,PRE,gcids,cvmap))
  for(sizeType k=0; k<nrKey; k++) {
    groups.get(keys[k],gcids);
    sizeType nrGC=(sizeType)gcids.size(),nrRow=nrGC*8;

    //assemble indices
    cvmap.clear();
    for(sizeType i=0; i<nrGC; i++) {
      const DHMCell& cell=mesh.getCell(gcids[i]);
      for(char v=0; v<8; v++)
        if(!boost::dynamic_pointer_cast<DHMConstrainedVertex>(cell._verts[v]))
          cvmap[cell._verts[v]->_index]=i*8+v;
        else nrRow++;
    }

    //build hess/grad
    RHS.setZero(nrRow);
    LHS.setZero(nrRow,nrRow);
    for(sizeType i=0; i<nrGC; i++) {
      DHMVelConsistencyEnergy e(mesh.getCell(gcids[i]),NULL,DHMVelConsistencyEnergy::GRAD::Zero());
      RHS.block<8,1>(i*8,0)=-e.gradIntFast(pts,weights,gcids[i]*nrP*3,samples);
      LHS.block<8,8>(i*8,i*8)=e.hessInt();
    }

    //build constraints
    sizeType off=nrGC*8;
    for(sizeType i=0; i<nrGC; i++) {
      const DHMCell& cell=mesh.getCell(gcids[i]);
      boost::shared_ptr<DHMConstrainedVertex> cv;
      for(char v=0; v<8; v++)
        if(cv=boost::dynamic_pointer_cast<DHMConstrainedVertex>(cell._verts[v])) {
          LHS(off,i*8+v)=LHS(i*8+v,off)=-1.0f;
          for(char cvid=0; cvid<4; cvid++)
            if(cv->_verts[cvid]) {
              sizeType c=cvmap.find(cv->_verts[cvid]->_index)->second;
              LHS(off,c)=LHS(c,off)=cv->_weights[cvid];
            }
          off++;
        }
    }
    ASSERT(off == nrRow)

    //solve
    PRE=LHS.householderQr().solve(RHS);
    for(sizeType i=0; i<nrGC; i++)
      field.getCompRef(gcids[i])=PRE.block<8,1>(i*8,0).cast<scalar>();
  }
}
void DHMFieldConverter::convertAdaptiveDirect(DHMCVField& field,const VelFunc<scalar>& func)
{
  const DHMMesh& mesh=field._mesh;
  sizeType nrC=mesh.nrCell();
  sizeType nrP=GaussLegendreIntegral<scalarD,scalarD>::nrP3D(2);

  COL samples(nrC*nrP*3);
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& cell=mesh.getCell(i);
    for(sizeType p=0; p<nrP; p++) {
      Vec3 pt=GaussLegendreIntegral<scalarD,scalarD>::point3D(2,p).cast<scalar>();
      samples.block<3,1>((i*nrP+p)*3,0)=func(DHMMapping::M(cell,pt));
    }
  }
  convertAdaptiveDirect(field,samples);
}

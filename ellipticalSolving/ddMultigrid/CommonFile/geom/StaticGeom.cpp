#include "../CollisionDetection.h"
#include "../MakeMesh.h"
#include "StaticGeom.h"
#include "StaticGeomCell.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

//Geom
StaticGeomCell::StaticGeomCell(sizeType type):Serializable(type) {}
StaticGeomCell::StaticGeomCell(const Mat4& T,sizeType dim,sizeType type)
  :Serializable(type),_T(T),_invT(T.inverse()),_dim(dim) {}
void StaticGeomCell::getMesh(ObjMeshTpl<scalar>& mesh,bool ref) const
{
  mesh.getV()=_vss;
  mesh.getI()=_iss;
  mesh.setDim((int)_dim);
  if(!ref) {
    mesh.getT()=_T.block<3,3>(0,0);
    mesh.getPos()=_T.block<3,1>(0,3);
    mesh.applyTrans(Vec3::Zero());
  }
}
BBox<scalar> StaticGeomCell::getBB() const
{
  Vec3 pt;
  BBox<scalar> tmp=getBBInner(),ret;
  for(sizeType x=0; x<2; x++)
    for(sizeType y=0; y<2; y++)
      for(sizeType z=0; z<2; z++) {
        pt[0]=(x==0) ? tmp._minC[0] : tmp._maxC[0];
        pt[1]=(y==0) ? tmp._minC[1] : tmp._maxC[1];
        pt[2]=(z==0) ? tmp._minC[2] : tmp._maxC[2];
        ret.setUnion(transformHomo<scalar>(_T,pt));
      }
  return ret;
}
bool StaticGeomCell::dist(const Vec3& pt,Vec3& n) const
{
  Vec3 pt0=transformHomo<scalar>(_invT,pt);
  if(distInner(pt0,n)) {
    n=_T.block<3,3>(0,0)*n;
    return true;
  }
  return false;
}
bool StaticGeomCell::closest(const Vec3& pt,Vec3& n) const
{
  Vec3 pt0=transformHomo<scalar>(_invT,pt);
  bool inside=closestInner(pt0,n);
  n=_T.block<3,3>(0,0)*n;
  return inside;
}
scalar StaticGeomCell::rayQuery(Vec3 x0,Vec3 dir) const
{
  x0=transformHomo<scalar>(_invT,x0);
  dir=(_invT.block<3,3>(0,0)*dir).eval();
  return rayQueryInner(x0,dir);
}
bool StaticGeomCell::read(std::istream& is)
{
  readVector(_vss,is);
  readVector(_iss,is);
  readVector(_bvh,is);

  readBinaryData(_T,is);
  readBinaryData(_invT,is);
  readBinaryData(_dim,is);
  readBinaryData(_index,is);
  return is.good();
}
bool StaticGeomCell::write(std::ostream& os) const
{
  writeVector(_vss,os);
  writeVector(_iss,os);
  writeVector(_bvh,os);

  writeBinaryData(_T,os);
  writeBinaryData(_invT,os);
  writeBinaryData(_dim,os);
  writeBinaryData(_index,os);
  return os.good();
}
const Mat4& StaticGeomCell::getInvT() const
{
  return _invT;
}
const Mat4& StaticGeomCell::getT() const
{
  return _T;
}
void StaticGeomCell::setT(const Mat4& T)
{
  _T=T;
  _invT=T.inverse();
}
sizeType StaticGeomCell::dim() const
{
  return _dim;
}
bool StaticGeomCell::closestInner(const Vec3& pt,Vec3& n) const
{
  ASSERT_MSG(false,"Not supported: function closestInner!")
}
scalar StaticGeomCell::rayQueryInner(const Vec3& x0,const Vec3& dir) const
{
  ASSERT_MSG(false,"Not supported: function rayQueryInner!")
  return 1;
}
void StaticGeomCell::build()
{
  ObjMesh mesh;
  getMeshInner(mesh);
  _vss=mesh.getV();
  _iss=mesh.getI();
  _bvh.resize(_iss.size());
  for(sizeType i=0; i<(sizeType)_bvh.size(); i++) {
    Node<sizeType>& n=_bvh[i];
    n._nrCell=1;
    n._cell=i;
    n._bb.reset();
    for(sizeType v=0; v<_dim; v++)
      n._bb.setUnion(_vss[_iss[i][v]]);
    n._bb.enlargedEps(0.01f);
  }
  buildBVH<sizeType>(_bvh,_dim,-1);
}

//Static Geometry
PRJ_BEGIN
struct GeomCallback {
  GeomCallback(const Vec3& pt,Vec3& n):_pt(pt),_dist(ScalarUtil<scalar>::scalar_max),_n(n) {}
  bool validNode(const Node<boost::shared_ptr<StaticGeomCell> >& node) {
    return node._bb.contain(_pt);
  }
  void updateDist(const Node<boost::shared_ptr<StaticGeomCell> >& node) {
    Vec3 n;
    scalar dist;
    if(node._cell->dist(_pt,n) && (dist=n.norm()) < _dist) {
      _cell=node._cell;
      _dist=dist;
      _n=n;
    }
  }
  const Vec3& _pt;
  boost::shared_ptr<StaticGeomCell> _cell;
  scalar _dist;
  Vec3& _n;
};
struct LineCallback {
  LineCallback(const Vec3& x0,Vec3& dir,sizeType dim):_x0(x0),_dir(dir),_dim(dim) {}
  bool validNode(const Node<boost::shared_ptr<StaticGeomCell> >& node) {
    return true;//node._bb.intersect(_x0,_x0+_dir,_dim);
  }
  void updateDist(const Node<boost::shared_ptr<StaticGeomCell> >& node) {
    scalar s=node._cell->rayQuery(_x0,_dir);
    if(s < 1) {
      _dir*=s;
      _cell=node._cell;
    }
  }
  Vec3 _x0;
  Vec3& _dir;
  sizeType _dim;
  boost::shared_ptr<StaticGeomCell> _cell;
};
struct VertexCallback {
  VertexCallback(StaticGeomCallback& cb):_cb(cb) {}
  bool validNode(const Node<sizeType>& node) {
    if(_gA->dim() == 2) {
      OBB2D obb(_BTA.block<3,3>(0,0),_BTA.block<3,1>(0,3),node._bb);
      return obb.intersect(_bbA);
    } else {
      OBB3D obb(_BTA.block<3,3>(0,0),_BTA.block<3,1>(0,3),node._bb);
      return obb.intersect(_bbA);
    }
    return true;
  }
  void updateDist(const Node<sizeType>& node) {
    Vec3 n,vRef,vG;
    sizeType dim=_gA->dim();
    const Vec3i& I=_gB->_iss[node._cell];
    for(sizeType v=0; v<dim; v++)
      if(!_mask[I[v]]) {
        _mask[I[v]]=true;
        vG=transformHomo<scalar>(_gB->getT(),_gB->_vss[I[v]]);
        if(_gA->dist(vG,n))
          _cb.onCollideVertex(vG,n,_gB);
      }
  }
  void onCell(const Node<boost::shared_ptr<StaticGeomCell> >& nA,
              const Node<boost::shared_ptr<StaticGeomCell> >& nB) {
    _gA=nA._cell;
    _gB=nB._cell;
    _mask.assign(_gB->_vss.size(),false);
    _bbA=_gA->getBBInner();
    _BTA=_gA->getInvT()*_gB->getT();

    BVHQuery<sizeType,BBox<scalar> > queryNarrow(nB._cell->_bvh,nB._cell->dim(),-1);
    queryNarrow.pointQuery(*this);
  }
  StaticGeomCallback& _cb;
  boost::shared_ptr<StaticGeomCell> _gA,_gB;
  vector<bool> _mask;
  BBox<scalar> _bbA;
  Mat4 _BTA;
};
PRJ_END
DebugStaticGeomCallback::DebugStaticGeomCallback(const std::string& path):_os("DebugCallback",path,true) {}
void DebugStaticGeomCallback::onCollideVertex(const Vec3& x,const Vec3& n,boost::shared_ptr<StaticGeomCell> c)
{
  vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
  vss.push_back(x);
  vss.push_back(x+n);

  _os.setRelativeIndex();
  _os.appendPoints(vss.begin(),vss.end());
  _os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                  VTKWriter<scalar>::IteratorIndex<Vec3i>(1,2,0),
                  VTKWriter<scalar>::LINE,true);
}
StaticGeom::StaticGeom(sizeType dim):HasMagic(0xaabbccdd),_dim(dim)
{
  _bvh.reset(new vector<Node<boost::shared_ptr<StaticGeomCell> > >);
}
const vector<Node<boost::shared_ptr<StaticGeomCell> > >& StaticGeom::getBVH() const
{
  return *_bvh;
}
vector<Node<boost::shared_ptr<StaticGeomCell> > >& StaticGeom::getBVH()
{
  return *_bvh;
}
void StaticGeom::clear()
{
  _css.clear();
  _bvh->clear();
}
void StaticGeom::update()
{
  if(!_bvh)
    return;

  vector<Node<boost::shared_ptr<StaticGeomCell>,BBox<scalar> > >& bvh=*_bvh;
  for(sizeType i=0; i<(sizeType)bvh.size(); i++)
    if(bvh[i]._cell)
      bvh[i]._bb=bvh[i]._cell->getBB();

  BVHQuery<boost::shared_ptr<StaticGeomCell>,BBox<scalar> > handler(*_bvh,_dim,boost::shared_ptr<StaticGeomCell>());
  handler.updateBVH();
}
void StaticGeom::assemble()
{
  _bvh->clear();
  for(sizeType i=0; i<(sizeType)_css.size(); i++) {
    _css[i]->_index=i;
    Node<boost::shared_ptr<StaticGeomCell> > n;
    n._l=n._r=n._parent=-1;
    n._nrCell=1;
    n._cell=_css[i];
    n._bb=n._cell->getBB();
    _bvh->push_back(n);
  }
  buildBVH(*_bvh,_dim,boost::shared_ptr<StaticGeomCell>());
}
bool StaticGeom::dist(const Vec3& pt,Vec3& n,boost::shared_ptr<StaticGeomCell>& cell) const
{
  BVHQuery<boost::shared_ptr<StaticGeomCell> > query(*_bvh,_dim,boost::shared_ptr<StaticGeomCell>());
  GeomCallback g(pt,n);
  query.pointQuery(g);
  cell=g._cell;
  return g._dist < ScalarUtil<scalar>::scalar_max;
}
bool StaticGeom::rayQuery(const Vec3& pt0,Vec3& dir,boost::shared_ptr<StaticGeomCell>& cell,Vec3& r) const
{
  BVHQuery<boost::shared_ptr<StaticGeomCell> > query(*_bvh,_dim,boost::shared_ptr<StaticGeomCell>());
  LineCallback g(pt0,dir,_dim);
  query.pointQuery(g);
  if(g._cell) {
    cell=g._cell;
    r=transformHomo<scalar>(cell->getInvT(),g._x0+g._dir);
    return true;
  } else return false;
}
void StaticGeom::collideVertex(const StaticGeom& other,StaticGeomCallback& cb) const
{
  BVHQuery<boost::shared_ptr<StaticGeomCell>,BBox<scalar> > query(*_bvh,_dim,boost::shared_ptr<StaticGeomCell>());
  BVHQuery<boost::shared_ptr<StaticGeomCell>,BBox<scalar> > queryOther(*(other._bvh),_dim,boost::shared_ptr<StaticGeomCell>());
  VertexCallback vcb(cb);
  query.interBodyQuery(queryOther,vcb);
}
void StaticGeom::addGeomCell(boost::shared_ptr<StaticGeomCell> c)
{
  _css.push_back(c);
}
void StaticGeom::addGeomBox(const Mat4& trans,const BBox<scalar>& bb,scalar depth)
{
  Mat4 T=Mat4::Identity();
  T.block<3,1>(0,3)=(bb._maxC+bb._minC)/2.0f;
  addGeomBox(trans*T,bb.getExtent()*0.5f,depth);
}
void StaticGeom::addGeomBox(const Mat4& trans,const Vec3& ext,scalar depth)
{
  depth=std::max<scalar>(depth,0.0f);
  _css.push_back(boost::shared_ptr<StaticGeomCell>
                 (new BoxGeomCell(trans,_dim,ext,depth,depth > 0.0f)));
}
void StaticGeom::addGeomBox(const OBBTpl<scalar,2>& obb,scalar depth)
{
  ASSERT(_dim == 2)
  Mat4 m=Mat4::Identity();
  m.block<2,2>(0,0)=obb._rot;
  m.block<2,1>(0,3)=obb._trans;
  addGeomBox(m,Vec3(obb._ext[0],obb._ext[1],0.0f),depth);
}
void StaticGeom::addGeomBox(const OBBTpl<scalar,3>& obb,scalar depth)
{
  ASSERT(_dim == 3)
  Mat4 m=Mat4::Identity();
  m.block<3,3>(0,0)=obb._rot;
  m.block<3,1>(0,3)=obb._trans;
  addGeomBox(m,obb._ext,depth);
}
void StaticGeom::addGeomCylinder(const Mat4& trans,scalar rad,scalar y,bool capsule)
{
  if(capsule)
    _css.push_back(boost::shared_ptr<StaticGeomCell>(new CylinderGeomCell(trans,_dim,rad,y)));
  else _css.push_back(boost::shared_ptr<StaticGeomCell>(new CapsuleGeomCell(trans,_dim,rad,y)));
}
void StaticGeom::addGeomPlane(const Mat4& trans,const Vec4& plane,scalar ext)
{
  scalar alpha=-plane[3]/plane.block<3,1>(0,0).squaredNorm();
  Vec3 p0=plane.block<3,1>(0,0)*alpha;
  Quatf q;
  q.setFromTwoVectors(Vec3::Unit(1).cast<scalarF>(),plane.block<3,1>(0,0).normalized().cast<scalarF>());

  Mat4 T=Mat4::Identity();
  T.block<3,1>(0,3)=p0-plane.block<3,1>(0,0).normalized()*ext;
  T.block<3,3>(0,0)=q.cast<scalar>().toRotationMatrix();
  if(_dim == 3)
    _css.push_back(boost::shared_ptr<StaticGeomCell>(new CylinderGeomCell(trans*T,_dim,ext,ext)));
  else _css.push_back(boost::shared_ptr<StaticGeomCell>(new BoxGeomCell(trans*T,_dim,Vec3::Constant(ext))));
}
void StaticGeom::addGeomPlane(const Mat4& trans,const PlaneTpl<scalar>& plane,scalar ext)
{
  Vec4 p;
  p.block<3,1>(0,0)=plane._n;
  p[3]=-plane._n.dot(plane._x0);
  addGeomPlane(trans,p,ext);
}
void StaticGeom::addGeomSphere(const Vec3& ctr,scalar rad,scalar depth)
{
  Mat4 T=Mat4::Identity();
  T.block<3,1>(0,3)=ctr;
  depth=std::max<scalar>(depth,0.0f);
  _css.push_back(boost::shared_ptr<StaticGeomCell>(new SphereGeomCell(T,_dim,rad,depth,depth > 0.0f)));
}
void StaticGeom::addGeomMesh(const Mat4& trans,const ObjMeshTpl<scalar>& mesh,scalar depth,bool insideOut)
{
  _css.push_back(boost::shared_ptr<StaticGeomCell>(new ObjMeshGeomCell(trans,mesh,depth,insideOut)));
}
void StaticGeom::addGeomMesh(const Mat4& trans,const std::string& path,scalar depth,bool insideOut)
{
  ObjMeshTpl<scalar> mesh;
  boost::filesystem::ifstream is(path);
  mesh.read(is,false,false);
  mesh.smooth();
  mesh.makeUniform();
  addGeomMesh(trans,mesh,depth,insideOut);
}
void StaticGeom::addGeomSolidBox(const Mat4& trans,const Vec3& ext,scalar thick)
{
  ObjMesh mesh;
  if(_dim == 2)MakeMesh::makeBox2D(mesh,ext,thick);
  else MakeMesh::makeBox3D(mesh,ext,thick);
  addGeomMesh(trans,mesh);
}
void StaticGeom::addGeomSolidSphere(const Mat4& trans,scalar rad,scalar thick)
{
  ObjMesh mesh;
  if(_dim == 2)MakeMesh::makeSphere2D(mesh,rad,32,thick);
  else MakeMesh::makeSphere3D(mesh,rad,32,thick);
  addGeomMesh(trans,mesh);
}
//IO
void StaticGeom::writeVTK(const std::string& path) const
{
  ObjMeshTpl<scalar> mesh;
  VTKWriter<scalar> os("Geom",path,true);
  boost::filesystem::path components=boost::filesystem::path(path).parent_path()/"geomComponents/";
  boost::filesystem::create_directory(components);
  for(sizeType i=0; i<(sizeType)_css.size(); i++) {
    _css[i]->getMesh(mesh);
    ostringstream oss;
    oss << components.string() << "/comp" << i << ".obj";
    mesh.write(oss.str());
    boost::filesystem::ofstream povSS(boost::filesystem::path(oss.str()).replace_extension(".pov"));
    mesh.writePov(povSS,false);
    boost::filesystem::ofstream spovSS(boost::filesystem::path(oss.str()).replace_extension(".spov"));
    mesh.writePov(spovSS,true);
    mesh.writeVTK(os,false,false);
  }
}
void StaticGeom::writeBVH() const
{
  writeBVHByLevel<boost::shared_ptr<StaticGeomCell> >(*_bvh,boost::shared_ptr<StaticGeomCell>());
}
bool StaticGeom::write(std::ostream& os) const
{
  if(!HasMagic::writeMagic(os))
    return false;
  IOData dat;
  registerType(dat);
  writeBinaryData(_dim,os);
  writeVector(_css,os,&dat);
  writeVector(*_bvh,os,&dat);
  return os.good();
}
bool StaticGeom::read(std::istream& is)
{
  _css.clear();
  _bvh->clear();
  if(!HasMagic::readMagic(is))
    return false;

  IOData dat;
  registerType(dat);
  readBinaryData(_dim,is);
  readVector(_css,is,&dat);
  readVector(*_bvh,is,&dat);
  return is.good();
}
bool StaticGeom::write(const boost::shared_ptr<StaticGeomCell>& cell,std::ostream& os)
{
  IOData dat;
  registerType(dat);
  writeBinaryData(cell,os,&dat);
  return os.good();
}
bool StaticGeom::read(boost::shared_ptr<StaticGeomCell>& cell,std::istream& is)
{
  IOData dat;
  registerType(dat);
  readBinaryData(cell,is,&dat);
  return is.good();
}
void StaticGeom::registerType(IOData& dat)
{
  dat.registerType<BoxGeomCell>();
  dat.registerType<SphereGeomCell>();
  dat.registerType<CylinderGeomCell>();
  dat.registerType<CapsuleGeomCell>();
  dat.registerType<ObjMeshGeomCell>();
}
void StaticGeom::debugRayQuery(const std::string& path,boost::shared_ptr<StaticGeomCell> cell,sizeType nr)
{
  StaticGeom geom(cell->dim());
  for(sizeType i=0; i<nr; i++) {
    boost::shared_ptr<StaticGeomCell> newCell=
      boost::dynamic_pointer_cast<StaticGeomCell>(cell->copy());

    Mat4 T=Mat4::Identity();
    T.block<3,3>(0,0)=cell->dim() == 3 ?
                      makeRotation<scalar>(Vec3::Random()) :
                      makeRotation<scalar>(Vec3::Unit(2)*(scalar)rand()/(scalar)RAND_MAX);
    T.block(0,3,cell->dim(),1).setRandom();
    newCell->setT(T);
    geom.addGeomCell(newCell);
  }
  geom.assemble();
  {
    ostringstream oss;
    oss << "./data/geom" << path << ".vtk";
    geom.writeVTK(oss.str());
  }

  ostringstream oss;
  oss << "./data/ray" << path << ".vtk";
  VTKWriter<scalar> os("rayQuery",oss.str(),true);
  vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
  vector<scalar> css;

  Vec3 x,dir,r;
  BBox<scalar> bb=geom._bvh->back()._bb;
  scalar norm=bb.getExtent().norm();
  for(sizeType j=0; j<100; j++) {
    x.setZero();
    x.block(0,0,cell->dim(),1).setRandom();
    x=x.normalized()*norm*2;
    dir=Vec3::Zero()-x;

    boost::shared_ptr<StaticGeomCell> ICell;
    geom.rayQuery(x,dir,ICell,r);

    vss.push_back(x);
    vss.push_back(x+dir);

    css.push_back(ICell ? 1.0f : 0.0f);
    css.push_back(ICell ? 1.0f : 0.0f);
  }
  os.appendPoints(vss.begin(),vss.end());
  os.appendCustomPointData("intersect",css.begin(),css.end());
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>(vss.size()/2,2,0),
                 VTKWriter<scalar>::LINE);
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>(vss.size(),0,0),
                 VTKWriter<scalar>::POINT);
}
void StaticGeom::debugVertexQuery(const std::string& path,boost::shared_ptr<StaticGeomCell> cell,sizeType nr)
{
  StaticGeom geomA(cell->dim());
  StaticGeom geomB(cell->dim());
  for(sizeType i=0; i<nr; i++) {
    boost::shared_ptr<StaticGeomCell> newCell=
      boost::dynamic_pointer_cast<StaticGeomCell>(cell->copy());

    Mat4 T=Mat4::Identity();
    T.block<3,3>(0,0)=cell->dim() == 3 ?
                      makeRotation<scalar>(Vec3::Random()) :
                      makeRotation<scalar>(Vec3::Unit(2)*(scalar)rand()/(scalar)RAND_MAX);
    T.block(0,3,cell->dim(),1).setRandom();
    newCell->setT(T);
    if(i > nr/2)
      geomA.addGeomCell(newCell);
    else geomB.addGeomCell(newCell);
  }
  geomA.assemble();
  geomB.assemble();
  {
    ostringstream oss;
    oss << "./data/geomA" << path << ".vtk";
    geomA.writeVTK(oss.str());
  }
  {
    ostringstream oss;
    oss << "./data/geomB" << path << ".vtk";
    geomB.writeVTK(oss.str());
  }
  ostringstream oss;
  oss << "./data/vert" << path << ".vtk";
  DebugStaticGeomCallback dsgcb(oss.str());
  geomA.collideVertex(geomB,dsgcb);
}
void StaticGeom::debugRayQuery()
{
  boost::shared_ptr<BoxGeomCell> box2(new BoxGeomCell(Mat4::Identity(),2,Vec3(0.1f,0.1f,0.0f)));
  debugRayQuery("box2",box2);
  debugVertexQuery("box2",box2);
  boost::shared_ptr<BoxGeomCell> box3(new BoxGeomCell(Mat4::Identity(),3,Vec3(0.1f,0.1f,0.1f)));
  debugRayQuery("box3",box3);
  debugVertexQuery("box3",box3);
  boost::shared_ptr<SphereGeomCell> sphere2(new SphereGeomCell(Mat4::Identity(),2,0.1f));
  debugRayQuery("sphere2",sphere2);
  debugVertexQuery("sphere2",sphere2);
  boost::shared_ptr<SphereGeomCell> sphere3(new SphereGeomCell(Mat4::Identity(),3,0.1f));
  debugRayQuery("sphere3",sphere3);
  debugVertexQuery("sphere3",sphere3);
  boost::shared_ptr<CylinderGeomCell> cylinder3(new CylinderGeomCell(Mat4::Identity(),3,0.1f,0.2f));
  debugRayQuery("cylinder3",cylinder3);
  debugVertexQuery("cylinder3",cylinder3);
  boost::shared_ptr<CapsuleGeomCell> capsule2(new CapsuleGeomCell(Mat4::Identity(),2,0.1f,0.2f));
  debugRayQuery("capsule2",capsule2);
  debugVertexQuery("capsule2",capsule2);
  boost::shared_ptr<CapsuleGeomCell> capsule3(new CapsuleGeomCell(Mat4::Identity(),3,0.1f,0.2f));
  debugRayQuery("capsule3",capsule3);
  debugVertexQuery("capsule3",capsule3);
  ObjMesh mesh;
  boost::filesystem::ifstream mis("./data/bunny.obj");
  mesh.read(mis,false,false);
  mesh.getScale()=5.0f;
  mesh.applyTrans(Vec3::Zero());
  boost::shared_ptr<ObjMeshGeomCell> mesh3(new ObjMeshGeomCell(Mat4::Identity(),mesh,1000.0f,false));
  debugRayQuery("mesh",mesh3,4);
  debugVertexQuery("mesh",mesh3,16);
}

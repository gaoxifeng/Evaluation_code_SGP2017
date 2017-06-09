#ifndef STATIC_GEOM_H
#define STATIC_GEOM_H

#include "../MathBasic.h"
#include "BVHBuilder.h"

PRJ_BEGIN

template <typename T,typename BBOX> struct Node;
template <typename T> class ObjMeshTpl;
template <typename T,typename BBOX> class BVHQuery;
template <typename T,int dim> class OBBTpl;
template <typename T> class PlaneTpl;

struct StaticGeomCell : public Serializable {
  friend struct VertexCallback;
public:
  StaticGeomCell(sizeType type);
  StaticGeomCell(const Mat4& T,sizeType dim,sizeType type);
  virtual ~StaticGeomCell() {}
  template <typename T> void setUserData(T* data) {
    _userData=(void*)(data);
  }
  template <typename T> T* getUserData() {
    return (T*)(_userData);
  }
  void getMesh(ObjMeshTpl<scalar>& mesh,bool ref=false) const;
  virtual BBox<scalar> getBB() const;
  virtual bool dist(const Vec3& pt,Vec3& n) const;
  virtual bool closest(const Vec3& pt,Vec3& n) const;
  virtual scalar rayQuery(Vec3 x0,Vec3 dir) const;
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  const Mat4& getInvT() const;
  const Mat4& getT() const;
  void setT(const Mat4& T);
  sizeType dim() const;
  sizeType _index;
protected:
  virtual void getMeshInner(ObjMeshTpl<scalar>& mesh) const=0;
  virtual BBox<scalar> getBBInner() const=0;
  virtual bool distInner(const Vec3& pt,Vec3& n) const=0;
  virtual bool closestInner(const Vec3& pt,Vec3& n) const;
  virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  void build();
  //vertex
  std::vector<Vec3,Eigen::aligned_allocator<Vec3> > _vss;
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _iss;
  std::vector<Node<sizeType> > _bvh;
  //data
  Mat4 _T,_invT;
  sizeType _dim;
  void* _userData;
};
class StaticGeomCallback
{
public:
  virtual void onCollideVertex(const Vec3& x,const Vec3& n,boost::shared_ptr<StaticGeomCell> c) =0;
};
class DebugStaticGeomCallback : public StaticGeomCallback
{
public:
  DebugStaticGeomCallback(const std::string& path);
  virtual void onCollideVertex(const Vec3& x,const Vec3& n,boost::shared_ptr<StaticGeomCell> c);
  VTKWriter<scalar> _os;
};
class StaticGeom : public HasMagic
{
public:
  //method
  StaticGeom(sizeType dim);
  const vector<Node<boost::shared_ptr<StaticGeomCell>,BBox<scalar> > >& getBVH() const;
  vector<Node<boost::shared_ptr<StaticGeomCell>,BBox<scalar> > >& getBVH();
  sizeType nrG() const {
    return (sizeType)_css.size();
  }
  const StaticGeomCell& getG(sizeType i) const {
    return *(_css[i]);
  }
  boost::shared_ptr<StaticGeomCell> getGPtr(sizeType i) const {
    return _css[i];
  }
  StaticGeomCell& getG(sizeType i) {
    return *(_css[i]);
  }
  void clear();
  void update();
  void assemble();
  bool dist(const Vec3& pt,Vec3& n,boost::shared_ptr<StaticGeomCell>& cell) const;
  bool rayQuery(const Vec3& pt0,Vec3& dir,boost::shared_ptr<StaticGeomCell>& cell,Vec3& r) const;
  void collideVertex(const StaticGeom& other,StaticGeomCallback& cb) const;
  //geometry
  void addGeomCell(boost::shared_ptr<StaticGeomCell> c);
  void addGeomBox(const Mat4& trans,const BBox<scalar>& bb,scalar depth=0.0f);
  void addGeomBox(const Mat4& trans,const Vec3& ext,scalar depth=0.0f);
  void addGeomBox(const OBBTpl<scalar,2>& obb,scalar depth=0.0f);
  void addGeomBox(const OBBTpl<scalar,3>& obb,scalar depth=0.0f);
  void addGeomCylinder(const Mat4& trans,scalar rad,scalar y,bool capsule=false);
  void addGeomPlane(const Mat4& trans,const Vec4& plane,scalar ext=100.0f);
  void addGeomPlane(const Mat4& trans,const PlaneTpl<scalar>& plane,scalar ext=100.0f);
  void addGeomSphere(const Vec3& ctr,scalar rad,scalar depth=0.0f);
  void addGeomMesh(const Mat4& trans,const ObjMeshTpl<scalar>& mesh,scalar depth=0.0f,bool insideOut=false);
  void addGeomMesh(const Mat4& trans,const std::string& path,scalar depth=0.0f,bool insideOut=false);
  void addGeomSolidBox(const Mat4& trans,const Vec3& ext,scalar thick);
  void addGeomSolidSphere(const Mat4& trans,scalar rad,scalar thick);
  //IO
  void writeVTK(const std::string& path) const;
  void writeObj(const std::string& path) const;
  void writeBVH() const;
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  //tools
  static bool write(const boost::shared_ptr<StaticGeomCell>& cell,std::ostream& os);
  static bool read(boost::shared_ptr<StaticGeomCell>& cell,std::istream& is);
  static void registerType(IOData& dat);
  static void debugRayQuery(const std::string& path,boost::shared_ptr<StaticGeomCell> cell,sizeType nr=100);
  static void debugVertexQuery(const std::string& path,boost::shared_ptr<StaticGeomCell> cell,sizeType nr=100);
  static void debugRayQuery();
private:
  //mesh data
  vector<boost::shared_ptr<StaticGeomCell> > _css;
  boost::shared_ptr<vector<Node<boost::shared_ptr<StaticGeomCell>,BBox<scalar> > > > _bvh;
  sizeType _dim;
};

PRJ_END

#endif

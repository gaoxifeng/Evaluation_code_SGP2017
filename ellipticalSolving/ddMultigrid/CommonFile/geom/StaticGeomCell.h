#ifndef STATIC_GEOM_CELL_H
#define STATIC_GEOM_CELL_H

#include "StaticGeom.h"
#include "../GridBasic.h"

PRJ_BEGIN

template <typename T,typename TI,typename TG>struct Grid;
typedef Grid<scalar,scalar,vector<scalar,Eigen::aligned_allocator<scalar> > > ScalarField;

struct BoxGeomCell : public StaticGeomCell {
  BoxGeomCell();
  BoxGeomCell(const Mat4& T,sizeType dim,const Vec3& ext,scalar depth=0.0f,bool insideOut=false);
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  boost::shared_ptr<Serializable> copy() const;
protected:
  virtual void getMeshInner(ObjMeshTpl<scalar>& mesh) const;
  virtual BBox<scalar> getBBInner() const;
  virtual bool distInner(const Vec3& pt,Vec3& n) const;
  virtual bool closestInner(const Vec3& pt,Vec3& n) const;
  virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  Vec3 _ext;
  scalar _depth;
  bool _insideOut;
};
struct SphereGeomCell : public StaticGeomCell {
  friend struct CapsuleGeomCell;
  SphereGeomCell();
  SphereGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar depth=0.0f,bool insideOut=false);
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  boost::shared_ptr<Serializable> copy() const;
protected:
  virtual void getMeshInner(ObjMeshTpl<scalar>& mesh) const;
  virtual BBox<scalar> getBBInner() const;
  virtual bool distInner(const Vec3& pt,Vec3& n) const;
  virtual bool closestInner(const Vec3& pt,Vec3& n) const;
  virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  scalar _rad;
  scalar _depth;
  bool _insideOut;
};
struct CylinderGeomCell : public StaticGeomCell {
  CylinderGeomCell();
  CylinderGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar y);
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  boost::shared_ptr<Serializable> copy() const;
protected:
  CylinderGeomCell(const Mat4& T,sizeType dim,sizeType type);
  virtual void getMeshInner(ObjMeshTpl<scalar>& mesh) const;
  virtual BBox<scalar> getBBInner() const;
  virtual bool distInner(const Vec3& pt,Vec3& n) const;
  virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  scalar _rad,_y;
};
struct CapsuleGeomCell : public CylinderGeomCell {
  CapsuleGeomCell();
  CapsuleGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar y);
  boost::shared_ptr<Serializable> copy() const;
protected:
  virtual void getMeshInner(ObjMeshTpl<scalar>& mesh) const;
  virtual BBox<scalar> getBBInner() const;
  virtual bool distInner(const Vec3& pt,Vec3& n) const;
  virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
};
struct ObjMeshGeomCell : public StaticGeomCell {
  ObjMeshGeomCell();
  ObjMeshGeomCell(const Mat4& trans,const ObjMeshTpl<scalar>& mesh,scalar depth,bool insideOut);
  bool read(std::istream& is);
  bool write(std::ostream& os) const;
  boost::shared_ptr<Serializable> copy() const;
  scalar depth() const;
  void updateDist(const Node<sizeType>& node,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist,scalar* minDist) const;
  void buildLevelSet(scalar cellSz, scalar off);
  scalar distLevelSet(const Vec3& pos) const;
protected:
  virtual void getMeshInner(ObjMeshTpl<scalar>& mesh) const;
  virtual BBox<scalar> getBBInner() const;
  virtual bool distInner(const Vec3& pt,Vec3& n) const;
  virtual bool closestInner(const Vec3& pt,Vec3& n) const;
  virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  void calcMinDist2D(const Vec3i& I,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist,scalar* minDist) const;
  void calcMinDist3D(const Vec3i& I,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist,scalar* minDist) const;
  ScalarField _grid;
  char _insideOut;
  scalar _depth;
};

PRJ_END

#endif

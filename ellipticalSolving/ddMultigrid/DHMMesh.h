#ifndef DHM_MESH_H
#define DHM_MESH_H

#include "DHMTraits.h"
#include "CommonFile/ObjMesh.h"
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>

PRJ_BEGIN

//elements
struct DHMVertex : public Serializable {
  DHMVertex();
  DHMVertex(sizeType index);
  virtual bool read(istream& is,IOData* dat);
  virtual bool write(ostream& os,IOData* dat) const;
  virtual boost::shared_ptr<Serializable> copy() const;
  Vec3 _pos;
  sizeType _index;
  bool _isSurface;
protected:
  DHMVertex(sizeType type,sizeType index);
};
struct DHMCell : public DHMTraits, public Serializable {
  DHMCell();
  DHMCell(sizeType index);
  void writeVTK(const std::string& path) const;
  virtual bool read(istream& is,IOData* dat);
  virtual bool write(ostream& os,IOData* dat) const;
  virtual boost::shared_ptr<Serializable> copy() const;
  static Vec2i getEdgeId(sizeType i);
  static Vec4i getFaceId(sizeType i);
  Vec2i getEdge(sizeType i) const;
  Vec4i getFace(sizeType i) const;
  Vec4i getFace(const TCOORDS& tcoords,sizeType tid,char d) const;
  boost::shared_ptr<DHMVertex> _verts[8];
  sizeType _index,_layer;
};
struct DHMFace : public Serializable {
  DHMFace();
  Vec4i getId() const;
  virtual bool read(istream& is,IOData* dat);
  virtual bool write(ostream& os,IOData* dat) const;
  virtual boost::shared_ptr<Serializable> copy() const;
  boost::shared_ptr<DHMVertex> _verts[4];
  boost::shared_ptr<DHMCell> _cells[2];
};
//some very simple but useful topological operations
class DHMMesh : public DHMTraits, public DHMEnergyTraits<scalarD>, public Serializable
{
  friend class DHMAdaptiveMesh;
public:
  //getter
  DHMMesh();
  virtual ~DHMMesh() {}
  sizeType nrVert() const;
  sizeType nrVertNC() const;
  sizeType nrCell() const;
  sizeType nrFace() const;
  DHMVertex& getVert(sizeType i);
  const DHMVertex& getVert(sizeType i) const;
  DHMCell& getCell(sizeType i);
  const DHMCell& getCell(sizeType i) const;
  DHMFace& getFace(sizeType i);
  const DHMFace& getFace(sizeType i) const;
  const vector<sizeType>& getAdaptive() const;
  void getPos(COL& pos) const;
  void setPos(const COL& pos);
  //IO
  void assembleVIndex() const;
  void reset(const VERTPTRS &verts,const CELLPTRS &cells,const CELLIDS *labels,sizeType nrPadding=numeric_limits<sizeType>::max());
  bool readVTK(const std::string& path);
  void writeVTK(const std::string& path,bool TCoords=false,const vector<scalar>* data=NULL,const vector<scalar>* vdata=NULL) const;
  void writeTVTK(const std::string& path,scalar coef,const vector<scalar>* data=NULL) const;
  virtual bool read(istream& is);
  virtual bool write(ostream& os) const;
  void writePov(const std::string& path) const;
  void writeAdaptiveEdge(const BBox<scalar>& bb,sizeType lv,const std::string& path) const;
  ObjMesh getSurface(boost::unordered_map<sizeType,sizeType>* vertMap=NULL,boost::unordered_set<Vec2i,Hash>* edgeMap=NULL) const;
  //add padding by subdivision
  sizeType getPadding() const;
  void addPadding(sizeType nrL,bool postResample=true);
  void subdivide();
  //add layers to ensure multgrid coarsening safety
  void patchPOT(bool postResample=true);
  void buildTCoord(TCOORDS& tcoords) const;
  //params
  boost::shared_ptr<boost::property_tree::ptree> _tree;
protected:
  //assembler
  void ensureRHS();
  void assembleIndex();
  void buildFace(const CELLIDS *labels=NULL);
  void buildLayer();
  sizeType buildPadding(POSSET* svert=NULL) const;
  //subdivider
  void subdivide(const DHMCell& cell,const Vec4i& face,sizeType nrSlice,PMAP& padMap,CELLPTRS& cells);
  void addLayer(TCOORDS& tcoords,sizeType tid,char d,sizeType nrL);
  void deleteRedundant();
  //resample vertices according to arc-length
  void resample(const TCOORDS& tcoords,char d,sizeType f=0,sizeType t=0);
  Vec3 resample(const vector<sizeType>& que,scalar len) const;
  //data
  VERTPTRS _verts;
  CELLPTRS _cells;
  FACEPTRS _faces;
  vector<sizeType> _adaptive;
};

PRJ_END

#endif

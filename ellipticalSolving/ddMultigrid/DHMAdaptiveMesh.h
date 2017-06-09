#ifndef DHM_ADAPTIVE_MESH_H
#define DHM_ADAPTIVE_MESH_H

#include "DHMMesh.h"

PRJ_BEGIN

struct Hierarchy;
struct DHMConstrainedVertex : public DHMVertex {
  DHMConstrainedVertex();
  DHMConstrainedVertex(sizeType index);
  virtual bool read(istream& is,IOData* dat);
  virtual bool write(ostream& os,IOData* dat) const;
  virtual boost::shared_ptr<Serializable> copy() const;
  boost::shared_ptr<DHMVertex> _verts[4];
  Vec4 _weights;
};
class DHMAdaptiveMesh : public DHMTraits, public DHMEnergyTraits<scalarD>
{
public:
  DHMAdaptiveMesh(Hierarchy& hier);
  const DHMMesh& getBaseMesh() const;
  void writeVTKLevel(const std::string& path,sizeType lv);
  void clearAdaptiveMesh();
  void createAdaptiveMesh(const BBox<scalar>& bb,sizeType lv);
  //we have several useful conventions for this call
  //convention: DHMVertex are put before any DHMConstrainedVertex
  //convention: BaseMesh DHMVertex and DHMCell are reduced
  //convention: subdivided DHMCell::_index are used to indicate BaseMesh cell index
  //convention: subdivided DHMCell::_layer/_index are useless in their conventional sense
  //convention: subdivided DHMVertex::_index are in order but DHMConstrainedVertex::_index are not
  void createAdaptiveMesh(const vector<sizeType>& cells,sizeType lv);
  void parityCheck() const;
protected:
  //data
  Hierarchy& _hier;
};

PRJ_END

#endif

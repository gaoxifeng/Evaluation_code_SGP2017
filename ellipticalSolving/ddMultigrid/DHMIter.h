#ifndef DHM_ITER_H
#define DHM_ITER_H

#include "DHMMesh.h"

PRJ_BEGIN

//IO Iterator
template <typename RET>
struct ElemIter {
  typedef RET value_type;
  ElemIter(const DHMMesh& mesh,sizeType id):_mesh(mesh),_id(id) {}
  void operator++() {
    _id++;
  }
  bool operator!=(const ElemIter& other) const {
    return _id != other._id;
  }
  virtual value_type operator*() const=0;
  const DHMMesh& _mesh;
  sizeType _id;
};
struct VertPosIter : public ElemIter<Vec3> {
  VertPosIter(const DHMMesh& mesh,sizeType id):ElemIter(mesh,id) {}
  value_type operator*() const {
    return _mesh.getVert(_id)._pos;
  }
};
struct VertTypeIter : public ElemIter<scalar> {
  VertTypeIter(const DHMMesh& mesh,sizeType id,sizeType baseNr=0):ElemIter(mesh,id),_baseNr(baseNr) {}
  value_type operator*() const {
    if(dynamic_cast<const DHMConstrainedVertex*>(&(_mesh.getVert(_id))))
      return 2.0f;
    else if(_mesh.getVert(_id)._index < _baseNr)
      return 1.0f;
    else return 0.0f;
  }
  sizeType _baseNr;
};
struct CellIdIter : public ElemIter<scalar> {
  CellIdIter(const DHMMesh& mesh,sizeType id):ElemIter(mesh,id) {}
  value_type operator*() const {
    return (scalar)(_mesh.getCell(_id)._index);
  }
};
struct CellVidIter : public ElemIter<DHMTraits::CELLID> {
  CellVidIter(const DHMMesh& mesh,sizeType id):ElemIter(mesh,id) {}
  value_type operator*() const {
    value_type ret;
    for(char v=0; v<8; v++)
      ret[v]=_mesh.getCell(_id)._verts[v]->_index;
    return ret;
  }
};
struct CellLayerIter : public ElemIter<sizeType> {
  CellLayerIter(const DHMMesh& mesh,sizeType id):ElemIter(mesh,id) {}
  value_type operator*() const {
    return _mesh.getCell(_id)._layer;
  }
};
//comparator
struct LSSVertId {
  bool operator()(boost::shared_ptr<DHMVertex> a,boost::shared_ptr<DHMVertex> b) const {
    return a->_index<b->_index;
  }
};
struct LSSCellLayer {
  bool operator()(boost::shared_ptr<DHMCell> a,boost::shared_ptr<DHMCell> b) const {
    return a->_layer<b->_layer;
  }
};
struct LSSTCoordD : public DHMTraits {
  LSSTCoordD(const TCOORDS& tc,char d):_tc(tc),_d(d) {}
  sizeType operator()(sizeType i) const {
    return _tc[i][_d];
  }
  bool operator()(sizeType i,sizeType j) const {
    return operator()(i)<operator()(j);
  }
private:
  const TCOORDS& _tc;
  char _d;
};

PRJ_END

#endif
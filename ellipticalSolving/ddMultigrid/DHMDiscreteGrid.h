#ifndef DISCRETE_GRID_H
#define DISCRETE_GRID_H

#include "DHMMesh.h"
#include "ParticleSet.h"

PRJ_BEGIN

class DHMDiscreteGrid : public multi_unordered_map<sizeType,sizeType,false>
{
public:
  DHMDiscreteGrid(const ParticleSetN& pset,scalar rad) {
    _minC.setConstant(numeric_limits<scalar>::max());
    for(sizeType i=0; i<pset.size(); i++)
      _minC=compMin(_minC,pset[i]._pos);
    _eps0=rad;

    for(sizeType i=0; i<pset.size(); i++)
      addLeaf(getId(pset[i]._pos),i);
  }
  DHMDiscreteGrid(const DHMMesh& mesh) {
    //build grid parameter
    _eps0=0.0f;
    _minC.setConstant(numeric_limits<scalar>::max());
    sizeType nrLeaf=mesh.nrCell();
    _bbcs.resize(nrLeaf);
    for(sizeType i=0; i<nrLeaf; i++) {
      const DHMCell& cell=mesh.getCell(i);
      for(char v=0; v<8; v++)
        _bbcs[i].setUnion(cell._verts[v]->_pos);
      _eps0+=_bbcs[i].getExtent().cwiseAbs().sum()/3;
      _minC=compMin(_minC,_bbcs[i]._minC);
    }
    _eps0/=(scalar)nrLeaf;

    //build linked list
    for(sizeType i=0; i<nrLeaf; i++) {
      Vec3i ai=getId(_bbcs[i]._minC),bi=getId(_bbcs[i]._maxC),id;
      for(id[0]=ai[0]; id[0]<=bi[0]; id[0]++)
        for(id[1]=ai[1]; id[1]<=bi[1]; id[1]++)
          for(id[2]=ai[2]; id[2]<=bi[2]; id[2]++)
            addLeaf(id,i);
    }
  }
  template <typename DistCallback>
  void rangeQuery(const Vec3& pos,scalar rad,DistCallback& cb) const {
    scalar sqrDist;
    Vec3i ai=getId(pos-Vec3::Constant(rad)),bi=getId(pos+Vec3::Constant(rad)),id;
    boost::unordered_map<sizeType,sizeType>::const_iterator iter;
    for(id[0]=ai[0]; id[0]<=bi[0]; id[0]++)
      for(id[1]=ai[1]; id[1]<=bi[1]; id[1]++)
        for(id[2]=ai[2]; id[2]<=bi[2]; id[2]++)
          if((iter=_keys.find(getKey(id))) != _keys.end())
            updateDist(iter->second,pos,1.0f,cb,sqrDist);
  }
  template <typename DistCallback>
  void pointQuery(const Vec3& pos,scalar thres,DistCallback& cb,scalar& sqrDist) const {
    cb._cid=0;
    Vec3i ai=getId(pos),bi=ai,id=ai;
    sqrDist=numeric_limits<scalar>::max();
    boost::unordered_map<sizeType,sizeType>::const_iterator iter;

    if((iter=_keys.find(getKey(id))) != _keys.end() && updateDist(iter->second,pos,thres,cb,sqrDist))
      return;
    while(sqrDist == numeric_limits<scalar>::max()) {
      ai-=Vec3i::Ones();
      bi+=Vec3i::Ones();

#define ITER_FACE(D0,V0, D1,A1,B1, D2,A2,B2)	\
for(id[D0]=V0,id[D1]=(A1);id[D1]<=(B1);id[D1]++)	\
	for(id[D2]=(A2);id[D2]<=(B2);id[D2]++){	\
		if((iter=_keys.find(getKey(id))) != _keys.end() && updateDist(iter->second,pos,thres,cb,sqrDist))	\
			return;	\
	}
      ITER_FACE(0,ai[0], 1,ai[1]+1,bi[1]-1, 2,ai[2]+1,bi[2]-1)
      ITER_FACE(0,bi[0], 1,ai[1]+1,bi[1]-1, 2,ai[2]+1,bi[2]-1)

      ITER_FACE(1,ai[1], 0,ai[0],bi[0], 2,ai[2]+1,bi[2]-1)
      ITER_FACE(1,bi[1], 0,ai[0],bi[0], 2,ai[2]+1,bi[2]-1)

      ITER_FACE(2,ai[2], 0,ai[0],bi[0], 1,ai[1],bi[1])
      ITER_FACE(2,bi[2], 0,ai[0],bi[0], 1,ai[1],bi[1])
    }
  }
  scalar eps0() const {
    return _eps0;
  }
private:
  template <typename DistCallback>
  FORCE_INLINE bool updateDist(sizeType head,const Vec3& pos,scalar thres,DistCallback& cb,scalar& sqrDist) const {
    while(head != -1) {
      if(_bbcs.empty() || _bbcs[_lists[head].first].distToSqr(pos) < sqrDist) {
        cb.updateDist(_lists[head].first,sqrDist);
        if(sqrDist < thres)
          return true;
      }
      head=_lists[head].second;
    }
    return false;
  }
  FORCE_INLINE sizeType getKey(const Vec3i& id) const {
    if(id[0] < 0 || id[1] < 0 || id[2] < 0)
      return -1;
    else return (id[0]<<40)+(id[1]<<20)+id[2];
  }
  FORCE_INLINE Vec3i getId(Vec3 pt) const {
    return floor(Vec3((pt-_minC)/_eps0));
  }
  FORCE_INLINE void addLeaf(const Vec3i& id,sizeType i) {
    sizeType key=getKey(id);
    ASSERT(key >= 0)
    boost::unordered_map<sizeType,sizeType>::iterator iter=_keys.find(key);
    if(iter == _keys.end()) {
      _keys[key]=(sizeType)_lists.size();
      _lists.push_back(make_pair(i,-1));
    } else {
      _lists.push_back(make_pair(i,iter->second));
      iter->second=(sizeType)_lists.size()-1;
    }
  }
  //data
  vector<BBox<scalar> > _bbcs;
  scalar _eps0;
  Vec3 _minC;
};

PRJ_END

#endif
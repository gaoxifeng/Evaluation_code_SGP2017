#include "DHMParticleSet.h"
#include "DHMDiscreteGrid.h"
#include "DHMAdvection.h"

USE_PRJ_NAMESPACE

struct HasNeighCallback {
  HasNeighCallback(const DHMParticleSet& pset,const Vec3& pos,scalar rad):_pset(pset),_pos(pos),_rad(rad),_has(false) {}
  void updateDist(sizeType i,scalar& sqrDist) {
    if((_pset[i]._pos-_pos).norm() < _rad)
      _has=true;
    sqrDist=_has ? 0.0f : 1.0f;
  }
  const DHMParticleSet& _pset;
  const Vec3& _pos;
  scalar _rad;
  bool _has;
};
void DHMParticleSet::transferVel(DHMCVField& cv)
{

}
void DHMParticleSet::buildSpatialHash(scalar cellSz)
{
  _cells.reset(new DHMDiscreteGrid(*this,cellSz));
}
bool DHMParticleSet::hasNeigh(const Vec3& pos,scalar rad) const
{
  HasNeighCallback cb(*this,pos,rad);
  _cells->rangeQuery(pos,rad,cb);
  return cb._has;
}
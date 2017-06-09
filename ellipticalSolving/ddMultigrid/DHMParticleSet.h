#ifndef DHM_PARTICLE_SET_H
#define DHM_PARTICLE_SET_H

#include "ParticleSet.h"
#include "DHMField.h"

PRJ_BEGIN

class DHMAdvection;
class DHMDiscreteGrid;
class DHMParticleSet : public ParticleSetN
{
public:
  void transferVel(DHMCVField& cv);
  void buildSpatialHash(scalar cellSz);
  bool hasNeigh(const Vec3& pos,scalar rad) const;
  scalar _cRad;
private:
  boost::shared_ptr<DHMDiscreteGrid> _cells;
};

PRJ_END

#endif
#ifndef COLLISION_FUNC_H
#define COLLISION_FUNC_H

#include "ParticleSet.h"

PRJ_BEGIN

template<typename P_TYPE>
class CollisionFunc
{
public:
  virtual void operator()(const P_TYPE& p,const sizeType& id) =0;
};
template<typename VEC3>
struct ExtractPosDirect {
  static FORCE_INLINE const VEC3& extract(const VEC3& p) {
    return p;
  }
};
template<typename P_TYPE>
struct ExtractPosParticle {
  typedef typename ScalarUtil<typename P_TYPE::scalarType>::ScalarVec3 Vec3Type;
  static FORCE_INLINE const Vec3Type& extract(const P_TYPE& p) {
    return p._pos;
  }
};

PRJ_END

#endif
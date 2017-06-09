#ifndef DHM_COARSER_CONFORMAL_H
#define DHM_COARSER_CONFORMAL_H

#include "DHMCoarser.h"

PRJ_BEGIN

class DHMEnergyPool;
class DHMCoarserGeomConformal : public DHMCoarserGeom
{
public:
  DHMCoarserGeomConformal(DHMMesh& meshF,DHMMesh& meshC,const SMAT& P);
  void setParameter();
  bool updateWeight(DHMEnergyPool& EC) const;
  Mat3d findTransfer(const DHMCell& cellF,const DHMCell& cellC) const;
  virtual void coarsen();
};

PRJ_END

#endif

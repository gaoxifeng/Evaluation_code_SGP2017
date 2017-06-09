#ifndef DHM_SMOOTHER_H
#define DHM_SMOOTHER_H

#include "DHMMesh.h"
#include "solvers/LinearSolver.h"

PRJ_BEGIN

class DHMEnergyPool;
class DHMSmoother : public DHMTraits, public DHMEnergyTraits<scalarD>
{
public:
  DHMSmoother(DHMMesh& mesh);
  void setParameter();
  bool smooth();
protected:
  //make flip free
  bool updateWeight(DHMEnergyPool& EC) const;
  bool makeValid();
  //smooth surface/interior
  void buildK(SMAT& L) const;
  void smoothSurface();
  void smoothInterior(DHMEnergyPool& E);
  //data
  DHMMesh& _mesh;
  ObjMesh _smesh;
  boost::unordered_map<sizeType,sizeType> _surf;
  PCGSolver<scalarD,Kernel<scalarD>,NoPreconSolver<scalarD> > _sol;
};

PRJ_END

#endif
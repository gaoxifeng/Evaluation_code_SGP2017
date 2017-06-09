#ifndef DHM_WATER_H
#define DHM_WATER_H

#include "DHMLinearSolver.h"
#include "DHMAdvection.h"

PRJ_BEGIN

class DHMWater
{
public:
  DHMWater(Hierarchy& hier);
  void setDefaultParam(bool reset=true);
  void advance(scalar dt);
  //IO
  DHMParticleSet& getPSet();
  const DHMParticleSet& getPSet() const;
  void writeFrame(ostream& os) const;
  void readFrame(istream& is);
  void profileWrite() const;
  boost::property_tree::ptree _tree;
protected:
  //data
  boost::shared_ptr<DHMParticleSet> _pset;
  boost::shared_ptr<DHMCVField> _cv;
  boost::shared_ptr<DHMVSField> _lv;
  scalar _time;
  //solver
  Hierarchy _hier;
  boost::shared_ptr<DHMAdvection> _adv;
  boost::shared_ptr<DHMPressureSolver> _sol;
};

PRJ_END

#endif
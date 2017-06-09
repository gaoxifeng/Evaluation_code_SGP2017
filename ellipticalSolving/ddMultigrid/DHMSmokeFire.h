#ifndef DHM_SMOKE_FIRE_H
#define DHM_SMOKE_FIRE_H

#include "DHMLinearSolver.h"
#include "DHMAdvection.h"

PRJ_BEGIN

//Solve Homogeneous-Neumann-Laplacian Equation
class DHMLaplacianSolver : DHMEnergyTraits<scalarD>
{
public:
  DHMLaplacianSolver(const DHMMesh& mesh,DHMVSField& t);
  void advance(scalar dt,scalar kappa);
protected:
  boost::shared_ptr<DHMLinearSystem> _sys;
  DHMVSField& _t;
  scalar _kdt;
  SMAT _M;
};
//Solve Fire Equation
class DHMSmokeFire
{
public:
  DHMSmokeFire(Hierarchy& hier);
  void setDefaultParam(bool reset=true);
  void advance(scalar dt);
  //IO
  void writeFrame(ostream& os) const;
  void readFrame(istream& is);
  void profileWrite() const;
  const DHMVSField& getFuel() const;
  const DHMVSField& getSoot() const;
  const DHMVSField& getTemp() const;
  const DHMCVField& getVel() const;
  DHMVSField& getFuel();
  DHMVSField& getSoot();
  DHMVSField& getTemp();
  DHMCVField& getVel();
  boost::property_tree::ptree _tree;
  static const scalar ABSOLUTE;
protected:
  //data
  boost::shared_ptr<DHMCVField> _cv;
  boost::shared_ptr<DHMVSField> _f,_s,_t;
  scalar _time;
  //solver
  Hierarchy _hier;
  boost::shared_ptr<DHMAdvection> _adv;
  boost::shared_ptr<DHMPressureSolver> _sol;
  boost::shared_ptr<DHMLaplacianSolver> _scatter;
};

PRJ_END

#endif
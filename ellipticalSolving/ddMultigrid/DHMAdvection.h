#ifndef DHM_ADVECTION_H
#define DHM_ADVECTION_H

#include "DHMField.h"
#include "DHMParticleSet.h"
#include "DHMLinearSolver.h"
#include "CollisionDetection.h"

PRJ_BEGIN

//Source Setup
struct DHMSource {
  //data
  enum TYPE {
    SPHERE,
    BOX,
  };
  TYPE _type;
  Vec3 _vel;
  scalar _value;
  BBox<scalar> _bb;
  Sphere<scalar> _sphere;
  //parser
  static void readStamp(const boost::property_tree::ptree& tree,vector<scalar>& stamps);
  static bool isValidSource(const boost::property_tree::ptree& tree,scalar time);
  static void addSourcePoint(const std::vector<DHMSource>& s,DHMParticleSet& pset,scalar cRad);
  static bool checkSourcePoint(const Vec3& a,const Vec3& b,const std::vector<DHMSource>& s,scalarD& val);
  static void parseSingleSource(boost::property_tree::ptree& tree,scalar time,vector<DHMSource>& source,const std::string& name);
  static void parseSource(boost::property_tree::ptree& tree,scalar time,vector<DHMSource>& source,const std::string& name);
  static void writeSource(const boost::property_tree::ptree& tree);
  //creator
  static sizeType findNewId(boost::property_tree::ptree& tree);
  static sizeType addSphereSource(boost::property_tree::ptree& tree,const Vec3& ctr,const scalar rad,scalar val,const std::string& name);
  static sizeType addBoxSource(boost::property_tree::ptree& tree,const BBox<scalar>& bb,scalar val,const std::string& name);
  static void setSourceParticleVel(boost::property_tree::ptree& tree,sizeType id,const Vec3& vel);
  //modifier
  static void addSourceDuration(boost::property_tree::ptree& tree,sizeType i,scalar t0,scalar t1);
  static void setSourceLoop(boost::property_tree::ptree& tree,sizeType i,bool loop);
};
//Advection
enum DHMStepper {
  EULER,
  RK2,
  RK4,
};
class DHMDiscreteGrid;
class DHMAdvection : public DHMTraits, public DHMEnergyTraits<scalar>
{
public:
  typedef DHMCVField *DHMCVFieldPTR;
  typedef DHMVSField *DHMVSFieldPTR;
  typedef const DHMVSField *CDHMVSFieldPTR;
  struct AdvPoint {
    Vec3 _pos,_crd;
    sizeType _cid;
  };
  DHMAdvection(const DHMMesh& mesh);
  void setAdaptiveMesh(const DHMMesh* meshA,sizeType lv);
  void setUseAdaptiveBasis(bool use);
  void advectParticleSet(const DHMCVField& cv,const DHMVSField& lv,DHMParticleSet& pset,scalar dt,const DHMCVField* cvOld=NULL,scalar FLIPWeight=0.0f) const;
  void advectMacCormack(const DHMCVField& cv,DHMCVFieldPTR cvOut,CDHMVSFieldPTR* vs,DHMVSFieldPTR* vsOut,sizeType nr,scalar dt,DHMStepper stepper=EULER) const;
  void advectSemiLagrangian(const DHMCVField& cv,DHMCVFieldPTR cvOut,CDHMVSFieldPTR* vs,DHMVSFieldPTR* vsOut,sizeType nr,scalar dt,DHMStepper stepper=EULER,Matd* vsMINMAX=NULL) const;
  bool trace(const Vec3& pos,sizeType& cid,Vec3& crd) const;
  const vector<AdvPoint>& getSamplePoints() const;
  sizeType getNrSamplePerCell() const;
  scalar getCRad() const;
  void debug() const;
  void debug(const std::string& path,const Matd& vsMINMAX) const;
  vector<vector<DHMSource> > _source;
protected:
  //data
  const DHMMesh& _mesh;
  boost::shared_ptr<DHMLinearSystem> _sys;
  boost::shared_ptr<DHMDiscreteGrid> _cells;
  bool _useAdaptiveBasis;
  vector<AdvPoint> _pts;
  sizeType _dim;
};
//Mass Consistency Energy
class DHMScalarConsistencyEnergy : public DHMEnergy<8>
{
public:
  DHMScalarConsistencyEnergy(const DHMCell& cell,const ImplicitFunc<scalar>* adv=NULL,const DHMVSField* vsOut=NULL);
  GRAD grad(const Vec3d& pos) const;
  GRAD gradStencil(const Vec3d& pos) const;
  HESS hess(const Vec3d& pos) const;
  const ImplicitFunc<scalar>* _adv;
  const DHMVSField* _vsOut;
};

PRJ_END

#endif
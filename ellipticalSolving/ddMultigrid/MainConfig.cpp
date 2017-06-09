#include "MainConfig.h"
#include "DHMAdaptiveMesh.h"
#include "DHMSmoother.h"
#include "DHMCoarser.h"
#include "DHMEnergy.h"
#include "DHMUtil.h"

#include "DHMAdvection.h"
#include "DHMLinearSolver.h"
#include "DHMSmokeFire.h"
#include "DHMWater.h"

#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

class TestDirSFunc : public ImplicitFunc<scalar>
{
public:
  TestDirSFunc(const Vec3& dir):_dir(dir) {}
  virtual scalar operator()(const Vec3& pos) const {
    return _dir.dot(pos);
  }
  const Vec3 _dir;
};
class TestSphereSFunc : public ImplicitFunc<scalar>
{
public:
  TestSphereSFunc(const Vec3& ctr,scalar rad,bool binary)
    :_ctr(ctr),_rad(rad),_binary(binary) {}
  virtual scalar operator()(const Vec3& pos) const {
    if(_binary)
      return (pos-_ctr).norm() < _rad ? 1.0f : 0.0f;
    else return (pos-_ctr).norm()-_rad;
  }
  Vec3 _ctr;
  scalar _rad;
  bool _binary;
};
class TestCubeSFunc : public ImplicitFunc<scalar>
{
public:
  virtual scalar operator()(const Vec3& pos) const {
    return 1.0f;
  }
};

class TestRotVFunc : public VelFunc<scalar>
{
public:
  TestRotVFunc(const Vec3& ctr,const Vec3& dir):_ctr(ctr),_dir(dir) {}
  virtual Vec3 operator()(const Vec3& pos) const {
    return _dir.cross(pos-_ctr);
  }
  const Vec3 _ctr,_dir;
};
class TestCubeVFunc : public VelFunc<scalar>
{
public:
  virtual Vec3 operator()(const Vec3& pos) const {
    if(pos[0] > 0.501f)
      return -Vec3::Unit(0)*0.002f;
    else if(pos[0] < 0.499f)
      return Vec3::Unit(0)*0.002f;
    else return Vec3::Zero();
  }
};
class TestSphereVFunc : public VelFunc<scalar>
{
public:
  virtual Vec3 operator()(const Vec3& pos) const {
    if( pos[0] >= -0.375f && pos[0] <= 0.375f &&
        pos[1] >= -0.375f && pos[1] <= 0.375f &&
        pos[2] >= -0.375f && pos[2] <= 0.375f)
      return Vec3::Unit(0)*0.1f;
    else return Vec3::Zero();
    //return Vec3::Unit(0)*0.1f;
  }
};
class TestSinCosVFunc : public VelFunc<scalar>
{
public:
  virtual Vec3 operator()(const Vec3& pos) const {
    return Vec3::Unit(0)*std::sin(pos[0]*M_PI*16);
  }
};
class TestKnotVFunc : public VelFunc<scalar>
{
public:
  virtual Vec3 operator()(const Vec3& pos) const {
    if((pos-Vec3(-47.0f,0.0f,0.0f)).norm() < 6.0f)
      return Vec3::Unit(1)+Vec3::Unit(2);
    else return Vec3::Zero();

    //if( pos[0] <= -35.0f)
    //    return Vec3::Unit(2)*10.0f;
    //else return Vec3::Zero();
  }
};
class TestDancingVFunc : public VelFunc<scalar>
{
public:
  virtual Vec3 operator()(const Vec3& pos) const {
    if((pos-Vec3(30.0f,0.0f,15.0f)).norm() < 10.0f)
      return Vec3::Unit(1);
    else return Vec3::Zero();

    //if( pos[0] <= -35.0f)
    //    return Vec3::Unit(2)*10.0f;
    //else return Vec3::Zero();
  }
};
class ConstantVFunc : public VelFunc<scalar>
{
public:
  virtual Vec3 operator()(const Vec3& pos) const {
    //return Vec3(pos.z(),pos.x(),0.0f)+Vec3::Ones();
    return Vec3::Unit(0)*0.1f;
  }
};
class TwoSideVFunc : public VelFunc<scalar>
{
public:
  virtual Vec3 operator()(const Vec3& pos) const {
    //if(pos[0] < 0.0f)
    //	return Vec3::Unit(2)*(1.0f-abs(pos[0]))*0.1f;
    //else return Vec3::Unit(2)*(abs(pos[0])-1.0f)*0.1f;

    if(pos[0] < 0.0f)
      return Vec3::Unit(2)*0.1f;
    else return -Vec3::Unit(2)*0.1f;
  }
};

int MainConfig::mainBuildHierarchy(int argc,char** argv)
{
  if(argc < 1) {
    INFO("Usage BuildHierarchy: [pathIn]")
    system("pause");
    return 1;
  }

  std::string pathIn(argv[0]);
  std::string parent=boost::filesystem::path(pathIn).parent_path().string();
  //DHMEnergyPool::debugEnergy();

  Hierarchy mesh;
  mesh.readVTK(pathIn);
  readParams(*(mesh._mesh.front()->_tree),"mesh",argc,argv);
  DHMCoarser coarser(mesh);
  coarser.coarsen();
  mesh.writeVTK(parent+"/MG",false);
  mesh.write(boost::filesystem::ofstream(parent+"/MG.dat",ios::binary));
  system("pause");
  return 0;
}
int MainConfig::mainSmooth(int argc,char** argv)
{
  if(argc < 2) {
    INFO("Usage Smooth: [pathIn] [pathOut]")
    system("pause");
    return 1;
  }

  DHMMesh mesh;
  mesh.readVTK(argv[0]);
  readParams(*(mesh._tree),"mesh",argc,argv);
  DHMSmoother smoother(mesh);
  if(!smoother.smooth()) {
    INFO("Smoother Failed!")
    system("pause");
    return 1;
  }
  mesh.writeVTK(argv[1]);
  system("pause");
  return 0;
}
int MainConfig::mainSubdivide(int argc,char** argv)
{
  if(argc < 2) {
    INFO("Usage Subdivide: [pathIn] [pathOut] (nrPadding)")
    system("pause");
    return 1;
  }

  sizeType nrPadding=0;
  if(argc > 2)
    nrPadding=readCommand<sizeType>(argv[2]);

  DHMMesh mesh;
  mesh.readVTK(argv[0]);
  if(nrPadding == 0)
    mesh.subdivide();
  else if(nrPadding > 0) {
    mesh.addPadding(nrPadding,true);
    mesh.patchPOT(true);
  } else if(nrPadding < 0) {
    mesh.addPadding(-nrPadding,false);
    mesh.patchPOT(false);
  }
  mesh.writeVTK(argv[1]);
  system("pause");
  return 0;
}
int MainConfig::mainMakeMesh(int argc,char** argv)
{
  if(argc != 2) {
    INFO("Usage MakeMesh: [nrPad] [nrCell]")
    system("pause");
    return 1;
  }

  sizeType nrPad=readCommand<sizeType>(argv[0]);
  sizeType nrCell=readCommand<sizeType>(argv[1]);
  if(nrPad == 0)
    DHMMeshMaker::makeCube("./cube.vtk",BBox<scalar>(Vec3::Zero(),Vec3::Ones()),Vec3i::Constant(nrCell));
  else DHMMeshMaker::makeSphere("./sphere.vtk",1.0f,nrPad,Vec3i::Constant(nrCell));
  system("pause");
  return 0;
}
int MainConfig::mainAdaptive(int argc,char** argv)
{
  if(argc != 4) {
    INFO("Usage Adaptive: [pathIn] [minC] [maxC] [lv]")
    system("pause");
    return 1;
  }

  Vec3 minC=readVec3(argv[1]);
  Vec3 maxC=readVec3(argv[2]);
  sizeType lv=readCommand<sizeType>(argv[3]);

  Hierarchy hier;
  if(boost::filesystem::path(argv[0]).extension() == ".vtk")
    hier.readVTK(argv[0]);
  else hier.read(boost::filesystem::ifstream(argv[0],ios::binary));

  if(lv > 0) {
    DHMAdaptiveMesh mesh(hier);
    mesh.createAdaptiveMesh(BBox<scalar>(minC,maxC),lv);
    mesh.parityCheck();
    hier.write(boost::filesystem::ofstream("./MG.dat",ios::binary));
    for(sizeType l=0; l<lv; l++) {
      ostringstream oss;
      oss << "./level" << l << ".vtk";
      mesh.writeVTKLevel(oss.str(),l);
    }
  } else {
    lv=-lv;
    sizeType i=0;
    for(sizeType l=lv; l>0; l--) {
      ostringstream oss;
      oss << "./level" << (i++) << ".sedge";
      hier._mesh.front()->writeAdaptiveEdge(BBox<scalar>(minC,maxC),l,oss.str());
    }
    for(sizeType l=0; l<(sizeType)hier._mesh.size(); l++) {
      ostringstream oss;
      oss << "./level" << (i++) << ".sedge";
      hier._mesh[l]->writePov(oss.str());
    }
  }
  system("pause");
  return 0;
}
int MainConfig::mainCube(int argc,char** argv)
{
  Hierarchy hier;
#define WRITE
#ifdef WRITE
  DHMMeshMaker::makeCube("./cube.vtk",BBox<scalar>(Vec3::Zero(),Vec3::Ones()),Vec3i(32,32,32));
  hier.readVTK("./cube.vtk");
  DHMCoarser coarser(hier);
  coarser.coarsen();
  //DHMAdaptiveMesh adaptive(hier);
  //adaptive.createAdaptiveMesh(BBox<scalar>(Vec3::Zero(),Vec3(1.0f,1.0f,0.49f)),1);
  hier.writeVTK("./cubeMG",false);
  hier.write(boost::filesystem::ofstream("./cubeMG/cubeMG.dat",ios::binary));
#else
  hier.read(boost::filesystem::ifstream("./cubeMG/cubeMG.dat",ios::binary));
  readParams(*(hier._mesh.front()->_tree),"mesh",argc,argv);
#endif

  boost::filesystem::create_directory("./cubeMG/");
  DHMSmokeFire sol(hier);
  readParams(sol._tree,"solver",argc,argv);
  DHMSource::addSphereSource(sol._tree,Vec3(0.5f,0.5f,0.1f),0.03f,1.0f,"sourceS");
  DHMSource::addSphereSource(sol._tree,Vec3(0.5f,0.5f,0.1f),0.03f,320.0f,"sourceT");
  DHMSource::addSphereSource(sol._tree,Vec3(0.5f,0.5f,0.9f),0.03f,1.0f,"sourceS");
  DHMSource::addSphereSource(sol._tree,Vec3(0.5f,0.5f,0.9f),0.03f,DHMSmokeFire::ABSOLUTE*2.0f-320.0f,"sourceT");

  for(sizeType i=0; i<1000; i++) {
    sol.advance(0.01f);

    ostringstream oss;
    oss << "./cubeMG/frm" << i << ".vtk";
    sol.getSoot().writeSVTK(oss.str());

    ostringstream oss2;
    oss2 << "./cubeMG/frmT" << i << ".vtk";
    sol.getTemp().writeSVTK(oss2.str());
  }

  system("pause");
  return 0;
}
int MainConfig::mainKnot(int argc,char** argv)
{
  Hierarchy hier;
  hier.read(boost::filesystem::ifstream("./knot/MGAdaptive.dat",ios::binary));
  readParams(*(hier._mesh.front()->_tree),"mesh",argc,argv);
  boost::filesystem::create_directory("./knotMG/");

  DHMSmokeFire sol(hier);
  sol._tree.put<scalar>("kappa",0.0f);
  readParams(sol._tree,"solver",argc,argv);
  DHMSource::addSphereSource(sol._tree,Vec3(-47.0f,0.0f,0.0f),3.0f,1.0f,"sourceS");
  DHMSource::addSphereSource(sol._tree,Vec3(-47.0f,0.0f,0.0f),3.0f,500.0f,"sourceT");

  Vec3 g(0.0f,-1.0f,-1.0f);
  g=g.normalized()*9.81f;
  sol._tree.put<scalar>("gx",g[0]);
  sol._tree.put<scalar>("gy",g[1]);
  sol._tree.put<scalar>("gz",g[2]);

  for(sizeType i=0; i<1000; i++) {
    sol.advance(0.05f);
    ostringstream oss;
    oss << "./knotMG/frm" << i << ".vtk";
    sol.getSoot().writeSVTK(oss.str());
  }

  system("pause");
  return 0;
}
int MainConfig::mainSphere(int argc,char** argv)
{
  /*DHMMesh mesh;
  mesh.readVTK("./sphereMG/hexSmooth.vtk");

  vector<scalar> data(mesh.nrVert());
  for(sizeType i=0;i<mesh.nrVert();i++)
  	if(mesh.getVert(i)._pos[0] < 0)
  		data[i]=1.0f;
  	else data[i]=0.0f;
  mesh.writeVTK("e:/sphere.vtk",false,NULL,&data);*/

  Hierarchy hier;
//#define WRITE
#ifdef WRITE
  hier.readVTK("./sphereMG/hexSmooth.vtk");
  hier._mesh.front()->subdivide();
  DHMCoarser coarser(hier);
  coarser.coarsen();
  hier.writeVTK("./sphereMG",false);
  hier.write(boost::filesystem::ofstream("./sphereMG/sphereMG.dat",ios::binary));
#else
  hier.read(boost::filesystem::ifstream("./sphereMG/sphereMG.dat",ios::binary));
  readParams(*(hier._mesh.front()->_tree),"mesh",argc,argv);
#endif

  boost::filesystem::create_directory("./sphereMGResult/");
  DHMSmokeFire sol(hier);
  //sol._tree.put<bool>("interpAdv",true);
  readParams(sol._tree,"solver",argc,argv);
  DHMSource::addSphereSource(sol._tree,Vec3(0.0f,0.0f,-0.7f),0.1f,1.0f,"sourceS");
  DHMSource::addSphereSource(sol._tree,Vec3(0.0f,0.0f,-0.7f),0.1f,320.0f,"sourceT");

  for(sizeType i=0; i<1000; i++) {
    sol.advance(0.05f);

    ostringstream oss;
    oss << "./sphereMGResult/frm" << i << ".vtk";
    sol.getSoot().writeSVTK(oss.str());

    ostringstream oss2;
    oss2 << "./sphereMGResult/frmT" << i << ".vtk";
    sol.getTemp().writeSVTK(oss2.str());
  }

  system("pause");
  return 0;
}
int MainConfig::mainWaterSource(int argc,char** argv)
{
  DHMMesh mesh;
  DHMMeshMaker::makeCube("./cube.vtk",BBox<scalar>(Vec3::Zero(),Vec3::Ones()),Vec3i(16,16,16));
  mesh.readVTK("./cube.vtk");

  DHMCVField cv(mesh);
  DHMVSField lv(mesh);
  cv.clear(0.0f);
  lv.clear(-1.0f);
  DHMFieldConverter::convertAdaptive(cv,TestRotVFunc(Vec3::Zero(),Vec3::Ones()));

  boost::property_tree::ptree tree;
  DHMSource::addSphereSource(tree,Vec3(0.5f,0.75f,0.5f),0.25f,0.0f,"sourceW");
  DHMAdvection adv(mesh);
  adv._source.resize(1);
  DHMSource::setSourceParticleVel(tree,0,Vec3::Ones());
  DHMSource::parseSource(tree,0.01f,adv._source[0],"sourceW");

  DHMParticleSet pset;
  pset._cRad=1.0f/16.0f/2.0f;
  boost::filesystem::create_directory("./frm");
  for(sizeType i=0; i<100; i++) {
    ostringstream oss;
    oss << "./frm/pset" << i << ".vtk";
    adv.advectParticleSet(cv,lv,pset,0.01f);
    pset.writeVTK(oss.str());
  }

  system("pause");
  return 0;
}
int MainConfig::mainWaterCube(int argc,char** argv)
{
  Hierarchy hier;
#define WRITE
#ifdef WRITE
  DHMMeshMaker::makeCube("./cube.vtk",BBox<scalar>(Vec3::Zero(),Vec3::Ones()),Vec3i(32,32,32));
  hier.readVTK("./cube.vtk");
  DHMCoarser coarser(hier);
  coarser.coarsen();
  //DHMAdaptiveMesh adaptive(hier);
  //adaptive.createAdaptiveMesh(BBox<scalar>(Vec3::Zero(),Vec3(1.0f,1.0f,0.49f)),1);
  hier.writeVTK("./cubeMG",false);
  hier.write(boost::filesystem::ofstream("./cubeMG/cubeMG.dat",ios::binary));
#else
  hier.read(boost::filesystem::ifstream("./cubeMG/cubeMG.dat",ios::binary));
  readParams(*(hier._mesh.front()->_tree),"mesh",argc,argv);
#endif

  boost::filesystem::create_directory("./cubeMG/");
  DHMWater sol(hier);
  readParams(sol._tree,"solver",argc,argv);
  DHMSource::addBoxSource(sol._tree,BBox<scalar>(Vec3::Zero(),Vec3(0.25f,0.25f,0.75f)),0.0f,"sourceW");

  for(sizeType i=0; i<1000; i++) {
    sol.advance(0.01f);

    ostringstream oss;
    oss << "./cubeMG/frmW" << i << ".vtk";
    sol.getPSet().writeSVTK(oss.str());
  }

  system("pause");
  return 0;
}
int MainConfig::main(int argc,char** argv)
{
  if(argc < 2) {
    INFO("Usage: exe [option]")
    system("pause");
    return 1;
  }

  istringstream iss(argv[1]);
  sizeType option;
  iss >> option;

  if(option == 0)
    return mainBuildHierarchy(argc-2,argv+2);
  if(option == 1)
    return mainSmooth(argc-2,argv+2);
  if(option == 2)
    return mainSubdivide(argc-2,argv+2);
  if(option == 3)
    return mainMakeMesh(argc-2,argv+2);
  if(option == 4)
    return mainAdaptive(argc-2,argv+2);
  if(option == 5)
    return mainCube(argc-2,argv+2);
  if(option == 6)
    return mainKnot(argc-2,argv+2);
  if(option == 7)
    return mainSphere(argc-2,argv+2);
  if(option == 8)
    return mainWaterSource(argc-2,argv+2);
  if(option == 9)
    return mainCubeWater(argc-2,argv+2);

  INFO("Unknown Command!")
  system("pause");
  return 0;
}
void MainConfig::readParams(boost::property_tree::ptree& tre,std::string prefix,int argc,char** argv)
{
  vector<std::string> args;
  for(sizeType i=0; i<argc; i++)
    args.push_back(argv[i]);
  readParams(tre,prefix,args);
}
void MainConfig::readParams(boost::property_tree::ptree& tre,std::string prefix,const vector<std::string>& args)
{
  char buf1[1024],buf2[1024],buf3[1024];
  for(int i=0; i<args.size(); i++) {
    std::string str=args[i];
    std::replace(str.begin(),str.end(),',',' ');
    if(sscanf_s(str.c_str(),"%s %s %s",buf1,1024,buf2,1024,buf3,1024) == 3)
      if(std::string(buf1) == prefix) {
        if(buf3 == std::string("true") ||
           buf3 == std::string("false")) {
          tre.put<bool>(buf2,buf3 == std::string("true"));
          continue;
        }
        try {
          sizeType val=boost::lexical_cast<sizeType>(buf3);
          tre.put<sizeType>(buf2,val);
          continue;
        } catch(...) {}
        try {
          scalar val=boost::lexical_cast<scalar>(buf3);
          tre.put<scalar>(buf2,val);
          continue;
        } catch(...) {}
        tre.put<std::string>(buf2,buf3);
      }
  }
}
int main(int argc,char** argv)
{
  return MainConfig::main(argc,argv);
}
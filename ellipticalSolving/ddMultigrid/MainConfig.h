#ifndef MAIN_CONFIG_H
#define MAIN_CONFIG_H

#include "MathBasic.h"
#include <vector>
#include <string>
#include <boost/property_tree/ptree.hpp>

PRJ_BEGIN

class MainConfig
{
public:
  static int mainBuildHierarchy(int argc,char** argv);
  static int mainSmooth(int argc,char** argv);
  static int mainSubdivide(int argc,char** argv);
  static int mainMakeMesh(int argc,char** argv);
  static int mainAdaptive(int argc,char** argv);
  static int mainCube(int argc,char** argv);
  static int mainKnot(int argc,char** argv);
  static int mainSphere(int argc,char** argv);
  static int mainWaterSource(int argc,char** argv);
  static int mainWaterCube(int argc,char** argv);
  static int main(int argc,char** argv);
protected:
  template <typename T>
  static T readCommand(char* argv) {
    T ret;
    std::istringstream iss(argv);
    iss >> ret;
    return ret;
  }
  static Vec3 readVec3(char* argv) {
    Vec3 ret;
    sscanf_s(argv,"%f,%f,%f",ret.data(),ret.data()+1,ret.data()+2);
    return ret;
  }
  static void readParams(boost::property_tree::ptree& tre,std::string prefix,int argc,char** argv);
  static void readParams(boost::property_tree::ptree& tre,std::string prefix,const vector<std::string>& args);
};

PRJ_END

#endif
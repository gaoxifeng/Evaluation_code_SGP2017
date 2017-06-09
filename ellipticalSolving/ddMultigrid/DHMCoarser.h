#ifndef DHM_COARSER_H
#define DHM_COARSER_H

#include "DHMMesh.h"
#include "CommonFile/DisjointSet.h"

PRJ_BEGIN

//multigrid builder with both topological and geometric coarsening
struct Group {
  Group():_id(numeric_limits<sizeType>::max()) {}
  DHMTraits::POSSET _cells;
  Vec4i _vid[8];
  Vec4i _corner;
  Vec3i _crd;
  sizeType _id;
};
struct Hierarchy : public DHMTraits, public DHMEnergyTraits<scalarD>, public Serializable {
  Hierarchy();
  void readVTK(const std::string& path);
  void writeVTK(const std::string& path,bool TCoords=false) const;
  void writePVTK(const std::string& path,bool TCoords=false) const;
  virtual bool read(istream& is);
  virtual bool write(ostream& os) const;
  vector<boost::shared_ptr<DHMMesh> > _mesh;
  vector<SMAT> _prolong;
};
class DHMCoarser : public DHMTraits, public DHMEnergyTraits<scalarD>
{
public:
  enum GEOM_COARSEN_ALGOR {
    INJECT,
    AVERAGE,
    CONFORMAL,
  };
  DHMCoarser(Hierarchy& hier);
  //algor
  void setParameter();
  void coarsen();
  void coarsenGeom();
protected:
  Hierarchy& _hier;
};
class DHMCoarserGeom : public DHMEnergyTraits<scalarD>, public DHMTraits
{
public:
  DHMCoarserGeom(DHMMesh& meshF,DHMMesh& meshC,const SMAT& P);
  virtual void coarsen() =0;
protected:
  void buildParentMap(boost::unordered_map<sizeType,sizeType>& pMap) const;
  DHMMesh& _meshF;
  DHMMesh& _meshC;
  const SMAT& _P;
};
class DHMCoarserGeomAverage : public DHMCoarserGeom
{
public:
  DHMCoarserGeomAverage(DHMMesh& meshF,DHMMesh& meshC,const SMAT& P);
  virtual void coarsen();
};
class DHMCoarserGeomInject : public DHMCoarserGeom
{
public:
  DHMCoarserGeomInject(DHMMesh& meshF,DHMMesh& meshC,const SMAT& P);
  virtual void coarsen();
};

PRJ_END

#endif

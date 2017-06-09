#include "DHMCoarser.h"
#include "DHMCoarserConformal.h"
#include "DHMUtil.h"
#include <boost/filesystem/operations.hpp>

PRJ_BEGIN
template <>
struct Zero<Group> {
  static Group value() {
    return Group();
  }
};
PRJ_END

USE_PRJ_NAMESPACE

Hierarchy::Hierarchy():Serializable(-1) {}
void Hierarchy::readVTK(const std::string& path)
{
  _mesh.clear();
  _prolong.clear();

  _mesh.push_back(boost::shared_ptr<DHMMesh>(new DHMMesh));
  _mesh.back()->readVTK(path);
}
void Hierarchy::writeVTK(const std::string& path,bool TCoords) const
{
  boost::filesystem::create_directory(path);
  //delete old/write new hierarchy
  boost::filesystem::directory_iterator end;
  if(TCoords) {
    for(boost::filesystem::directory_iterator it(path); it!=end; it++)
      if(beginsWith(it->path().filename().string(),"lvTH"))
        boost::filesystem::remove(it->path());
    scalar coef=1.0f;
    for(sizeType i=0; i<(sizeType)_mesh.size(); i++) {
      ostringstream oss;
      oss << path << "/lvTH" << i << ".vtk";
      _mesh[i]->writeTVTK(oss.str(),coef);
      coef*=2.0f;
    }
  } else {
    for(boost::filesystem::directory_iterator it(path); it!=end; it++)
      if(beginsWith(it->path().filename().string(),"lvH"))
        boost::filesystem::remove(it->path());
    for(sizeType i=0; i<(sizeType)_mesh.size(); i++) {
      ostringstream oss;
      oss << path << "/lvH" << i << ".vtk";
      _mesh[i]->writeVTK(oss.str(),true);
    }
  }
  //delete old prolongation matrix
  for(boost::filesystem::directory_iterator it(path); it!=end; it++)
    if(beginsWith(it->path().filename().string(),"lvP"))
      boost::filesystem::remove(it->path());
  //write new prolongation matrix
  for(sizeType i=0; i<(sizeType)_prolong.size(); i++) {
    ostringstream oss;
    oss << path << "/lvP" << i << ".dat";
    boost::filesystem::ofstream os(oss.str(),ios::binary);
    _prolong[i].write(os);
  }
}
void Hierarchy::writePVTK(const std::string& path,bool TCoords) const
{
  scalar coef=1.0f;
  for(sizeType i=0; i<(sizeType)_prolong.size(); i++) {
    ostringstream oss;
    oss << path << "/P" << i;
    boost::filesystem::create_directory(oss.str());

    const DHMMesh& meshF=*(_mesh[i]);
    const DHMMesh& meshC=*(_mesh[i+1]);
    const SMAT& P=_prolong[i];

    TCOORDS tcoordsF,tcoordsC;
    if(TCoords) {
      meshF.buildTCoord(tcoordsF);
      meshC.buildTCoord(tcoordsC);
    }
    for(sizeType j=0; j<meshF.nrCell(); j++) {
      vector<scalar> color;
      POSES pos;
      for(char v=0; v<8; v++) {
        sizeType vid=meshF.getCell(j)._verts[v]->_index;
        for(ConstSMIterator<scalarD> beg=P.begin(vid),end=P.end(vid); beg!=end; ++beg) {
          color.push_back(1.0f);
          color.push_back(0.0f);
          if(TCoords) {
            pos.push_back(tcoordsF[beg.row()].cast<scalar>()*coef);
            pos.push_back(tcoordsC[beg.col()].cast<scalar>()*coef*2.0f);
          } else {
            pos.push_back(meshF.getVert(beg.row())._pos);
            pos.push_back(meshC.getVert(beg.col())._pos);
          }
        }
      }

      ostringstream ossI;
      ossI << oss.str() << "/cell" << j << ".vtk";
      VTKWriter<scalar> os("prolong",ossI.str(),true);
      os.appendPoints(pos.begin(),pos.end());
      os.appendCustomPointData("color",color.begin(),color.end());
      os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                     VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)pos.size()/2,2,0),
                     VTKWriter<scalar>::LINE);
      os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,1,0),
                     VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)pos.size(),1,0),
                     VTKWriter<scalar>::POINT);
    }
    coef*=2.0f;
  }
}
bool Hierarchy::read(istream& is)
{
  sizeType nrM;
  readBinaryData(nrM,is);
  _mesh.resize(nrM);
  for(sizeType i=0; i<nrM; i++) {
    _mesh[i].reset(new DHMMesh);
    _mesh[i]->read(is);
  }

  sizeType nrP;
  readBinaryData(nrP,is);
  _prolong.resize(nrP);
  for(sizeType i=0; i<nrP; i++)
    _prolong[i].read(is);

  return is.good();
}
bool Hierarchy::write(ostream& os) const
{
  sizeType nrM=(sizeType)_mesh.size();
  writeBinaryData(nrM,os);
  for(sizeType i=0; i<nrM; i++)
    _mesh[i]->write(os);

  sizeType nrP=(sizeType)_prolong.size();
  writeBinaryData(nrP,os);
  for(sizeType i=0; i<nrP; i++)
    _prolong[i].write(os);

  return os.good();
}
//topological coarsening
void mergeVertex(DisjointSet<sizeType>& vGroups,const Group& G0,const Group& G1)
{
  if(G0._id == numeric_limits<sizeType>::max() || G1._id == numeric_limits<sizeType>::max())
    return;
  Vec4i ret;
  sizeType nr=0;
  for(char vi=0; vi<8; vi++)
    for(char vj=0; vj<8; vj++)
      if(G0._vid[vi] == G1._vid[vj]) {
        vGroups.joinSafe(G0._id*8+vi,G1._id*8+vj);
        ret[nr++]=G0._id*8+vi;
      }
  ASSERT(nr == 4)
}
void deleteRedundantCell(DisjointSet<Group>& cGroups,DisjointSet<sizeType>& vGroups)
{
  sizeType nrE=(sizeType)cGroups._elts.size();
  DHMTraits::CMAP cells;
  for(sizeType c=0; c<nrE; c++) {
    //check for valid group
    Group& G=cGroups._elts[c]._Int;
    if(G._id == numeric_limits<sizeType>::max())
      continue;
    //find cell id
    DHMTraits::CELLID cid;
    for(char v=0; v<8; v++)
      cid[v]=vGroups.find(G._id*8+v);
    std::sort(cid.data(),cid.data()+8);
    //check if it is unique
    DHMTraits::CCITER it=cells.find(cid);
    if(it != cells.end()) {
      //merge with that cell
      Group& GO=cGroups._elts[it->second]._Int;
      //silly assert
      ASSERT(GO._crd == G._crd && GO._corner == G._corner);
      for(char v=0; v<8; v++)
        ASSERT(GO._vid[v] == G._vid[v]);
      //merge the group
      GO._cells.insert(G._cells.begin(),G._cells.end());
      G._cells.clear();
      //invalidate this
      G._id=numeric_limits<sizeType>::max();
    } else cells[cid]=c;
  }
}
void addPRow(const DHMTraits::CELLID& vid,const Vec3i voff,DHMTraits::CELLID& PRow)
{
  ASSERT(compGE(voff,Vec3i::Zero()) && compLE(voff,Vec3i::Constant(2)));
  if(PRow[0] != numeric_limits<sizeType>::max())
    return;

  DHMTraits::CELLID dists;
  sizeType minDist=3;
  for(char v=0; v<8; v++) {
    Vec3i off(v&1?2:0,v&2?2:0,v&4?2:0);
    dists[v]=(off-voff).cwiseAbs().sum();
    minDist=std::min(minDist,dists[v]);
  }
  for(char v=0,neigh=0; v<8; v++)
    if(dists[v] == minDist)
      PRow[neigh++]=vid[v];
}
void buildP(const DHMTraits::CELLIDS& trips,DHMEnergyTraits<scalarD>::SMAT& P)
{
  DHMEnergyTraits<scalarD>::TRIPS tripsP;
  for(sizeType i=0; i<(sizeType)trips.size(); i++) {
    char nr=0;
    for(char v=0; v<8; v++,nr++)
      if(trips[i][v] == numeric_limits<sizeType>::max())
        break;

    scalar coef=1.0f/(scalar)nr;
    for(char v=0; v<nr; v++)
      tripsP.push_back(Eigen::Triplet<scalarD,sizeType>(i,trips[i][v],coef));
  }
  P.buildFromTripletsDepulicate(tripsP,0.0f);
}
bool coarsenTopo(const DHMMesh& meshF,DHMMesh& meshC,DHMEnergyTraits<scalarD>::SMAT& P)
{
  //build texture coords
  DHMTraits::TCOORDS tcoords;
  sizeType nrPadding=meshF.getPadding();
  ASSERT(countBits<sizeType>(nrPadding) <= 1)
  meshF.buildTCoord(tcoords);

  //check for termination condition
  BBox<sizeType> bb;
  for(sizeType i=0; i<meshF.nrCell(); i++)
    bb.setUnion(tcoords[i]);
  sizeType stopThres=meshF._tree->get<sizeType>("stopThres");
  if(bb.getExtent().minCoeff() < stopThres)
    return false;

  //find cell index/vertex layer
  vector<sizeType> vLayers(meshF.nrVert(),nrPadding);
  vector<Vec4i> cindices(meshF.nrCell(),Vec4i::Zero());
  for(sizeType i=0; i<meshF.nrCell(); i++) {
    const DHMCell& c=meshF.getCell(i);
    //sum over texture coordinates
    Vec4i& cid=cindices[i];
    for(char v=0; v<8; v++) {
      sizeType vid=c._verts[v]->_index;
      cid.block<3,1>(0,0)+=tcoords[vid];
      vLayers[vid]=std::min<sizeType>(vLayers[vid],-c._layer);
    }
    //extra dimension for padding
    if(c._layer < 0) {
      cid[3]=-c._layer*8+12;
      if(cid[0]%8 == 0)
        cid[3]+=0*nrPadding;
      else if(cid[1]%8 == 0)
        cid[3]+=8*nrPadding;
      else if(cid[2]%8 == 0)
        cid[3]+=16*nrPadding;
    }
  }

  //merge cell into coarse cell
  DisjointSet<Group> cGroups(meshF.nrCell());
  for(sizeType i=0; i<meshF.nrFace(); i++) {
    const DHMFace& face=meshF.getFace(i);
    if(!face._cells[1])
      continue;
    sizeType cid0=face._cells[0]->_index;
    sizeType cid1=face._cells[1]->_index;
    if(cindices[cid0]/16 == cindices[cid1]/16)
      cGroups.joinSafe(cid0,cid1);
  }
  sizeType nrG=0;
  for(sizeType i=0; i<meshF.nrCell(); i++) {
    //if we have only one padding layer
    //we collapse instead of merge layer
    if(nrPadding == 1 && meshF.getCell(i)._layer < 0)
      continue;
    sizeType root=cGroups.find(i);
    cGroups._elts[root]._Int._cells.insert(i);
    if(root == i)
      cGroups._elts[root]._Int._id=nrG++;
  }

  //determine coarse cell coordinate
  for(sizeType i=0; i<meshF.nrCell(); i++) {
    Group& G=cGroups._elts[i]._Int;
    if(G._id == numeric_limits<sizeType>::max())
      continue;
    const DHMTraits::POSSET& cells=G._cells;
    //create parent map
    Vec4i& corner=G._corner;
    Vec3i& crd=G._crd;
    //ASSERT(cells.size() <= 8)
    ASSERT(nrPadding <= 1 || cells.size() == 8) {
      const DHMCell& cell0=meshF.getCell(*(cells.begin()));
      //determine corner coordinate
      Vec4i cid=cindices[cell0._index];
      corner=cid/16*2;
      corner[3]=std::max<sizeType>((-cell0._layer-1)/2*2,0);
      //determine coordinate frame
      crd=Vec3i(0,1,2);
      if(cid[0]%8 == 0)
        crd=Vec3i(3,1,2);
      else if(cid[1]%8 == 0)
        crd=Vec3i(0,3,2);
      else if(cid[2]%8 == 0)
        crd=Vec3i(0,1,3);
    }
    for(char v=0; v<8; v++) {
      Vec4i& vt=G._vid[v]=corner;
      if(v&1)vt[crd[0]]+=2;
      if(v&2)vt[crd[1]]+=2;
      if(v&4)vt[crd[2]]+=2;
    }
  }

  //merge vertex for coarse cells
  DisjointSet<sizeType> vGroups(nrG*8);
  for(sizeType i=0; i<meshF.nrFace(); i++) {
    const DHMFace& face=meshF.getFace(i);
    if(!face._cells[1])
      continue;
    sizeType r0=cGroups.find(face._cells[0]->_index);
    sizeType r1=cGroups.find(face._cells[1]->_index);
    if(r0 != r1)
      mergeVertex(vGroups,cGroups._elts[r0]._Int,cGroups._elts[r1]._Int);
  }
  //delete cells whose 8 vertices are same
  deleteRedundantCell(cGroups,vGroups);

  //assemble coarse mesh
  DHMTraits::VERTPTRS cverts;
  for(sizeType i=0,vid=0; i<(sizeType)vGroups._elts.size(); i++)
    if(vGroups._elts[i]._p == i) {
      vGroups._elts[i]._Int=vid;
      cverts.push_back(boost::shared_ptr<DHMVertex>(new DHMVertex(vid)));
      vid++;
    }
  DHMTraits::CELLIDS labels;
  DHMTraits::CELLPTRS ccells;
  for(sizeType i=0; i<meshF.nrCell(); i++) {
    const Group& G=cGroups._elts[i]._Int;
    if(G._id == numeric_limits<sizeType>::max())
      continue;
    //add cell
    boost::shared_ptr<DHMCell> cell(new DHMCell((sizeType)ccells.size()));
    for(char v=0; v<8; v++) {
      sizeType root=vGroups.find(G._id*8+v);
      cell->_verts[v]=cverts[vGroups._elts[root]._Int];
    }
    ccells.push_back(cell);
    //add label
    DHMTraits::CELLID cid;
    cid.block<4,1>(0,0)=G._corner;
    cid.block<3,1>(4,0)=G._crd;
    labels.push_back(cid);
  }

  //assemble prolongation matrix
  DHMTraits::CELLID invalid=DHMTraits::CELLID::Constant(numeric_limits<sizeType>::max());
  DHMTraits::CELLIDS tripsP(meshF.nrVert(),invalid);
  for(sizeType i=0; i<meshF.nrCell(); i++) {
    const Group& G=cGroups._elts[i]._Int;
    if(G._id == numeric_limits<sizeType>::max())
      continue;
    const DHMTraits::POSSET& cells=G._cells;
    const Vec4i& corner=G._corner;
    const Vec3i& crd=G._crd;
    DHMTraits::CELLID vid;
    for(char v=0; v<8; v++) {
      sizeType root=vGroups.find(G._id*8+v);
      vid[v]=vGroups._elts[root]._Int;
    }
    for(DHMTraits::CPOSITER cbeg=cells.begin(),cend=cells.end(); cbeg!=cend; cbeg++) {
      const DHMCell& cell=meshF.getCell(*cbeg);
      for(char v=0; v<8; v++) {
        sizeType id=cell._verts[v]->_index;
        Vec4i voff=Vec4i(tcoords[id][0],tcoords[id][1],tcoords[id][2],vLayers[id])-corner;
        addPRow(vid,Vec3i(voff[crd[0]],voff[crd[1]],voff[crd[2]]),tripsP[id]);
      }
    }
  }
  //collapse padding if nrPadding == 1
  if(nrPadding == 1)
    for(sizeType i=0; i<meshF.nrCell(); i++) {
      const DHMCell& cell=meshF.getCell(i);
      if(cell._layer < 0)
        for(char e=0; e<12; e++) {
          Vec2i eid=cell.getEdge(e);
          if(tripsP[eid[0]] == invalid && tripsP[eid[1]] != invalid)
            tripsP[eid[0]]=tripsP[eid[1]];
          else if(tripsP[eid[1]] == invalid && tripsP[eid[0]] != invalid)
            tripsP[eid[1]]=tripsP[eid[0]];
        }
    }

  //build coarse mesh and prolongation matrix
  meshC.reset(cverts,ccells,&labels,nrPadding/2);
  meshC._tree=meshF._tree;
  P.resize(meshF.nrVert(),meshC.nrVert());
  buildP(tripsP,P);
  return true;
}
//DHMCoarser
DHMCoarser::DHMCoarser(Hierarchy& hier):_hier(hier)
{
  setParameter();
}
void DHMCoarser::setParameter()
{
  DHMMesh& mesh=*(_hier._mesh.front());
  putNoOverwrite<sizeType>(*(mesh._tree),"geomCoarsenAlgor",0);
  putNoOverwrite<sizeType>(*(mesh._tree),"stopThres",3);
}
void DHMCoarser::coarsen()
{
  _hier._mesh.resize(1);
  _hier._prolong.clear();
  INFO("Patch singularities to POWER-OF-TWO boundary!")
  _hier._mesh[0]->patchPOT();
  INFO("Start Topological Coarsening!")
  while(true) {
    SMAT P;
    boost::shared_ptr<DHMMesh> meshF=_hier._mesh.back();
    boost::shared_ptr<DHMMesh> meshC(new DHMMesh);
    if(!coarsenTopo(*meshF,*meshC,P))
      break;
    _hier._mesh.push_back(meshC);
    _hier._prolong.push_back(P);
  }
  //_hier.writeVTK("./MG",true);
  //_hier.writePVTK("./MG",true);
  INFO("Start Geometrical Coarsening!")
  coarsenGeom();
}
void DHMCoarser::coarsenGeom()
{
  vector<boost::shared_ptr<DHMMesh> >& mesh=_hier._mesh;
  vector<SMAT>& prolong=_hier._prolong;
  sizeType nrM=(sizeType)mesh.size();
  for(sizeType i=0; i<nrM-1; i++) {
    switch(_hier._mesh.front()->_tree->get<sizeType>("geomCoarsenAlgor")) {
    case INJECT:
      DHMCoarserGeomAverage(*(mesh[i]),*(mesh[i+1]),prolong[i]).coarsen();
      break;
    case AVERAGE:
      DHMCoarserGeomInject(*(mesh[i]),*(mesh[i+1]),prolong[i]).coarsen();
      break;
    case CONFORMAL:
      DHMCoarserGeomConformal(*(mesh[i]),*(mesh[i+1]),prolong[i]).coarsen();
      break;
    default:
      ASSERT_MSG(false,"Unknown Coarsening Algorithm!")
      break;
    }
  }
}
//geometrical coarsening
DHMCoarserGeom::DHMCoarserGeom(DHMMesh& meshF,DHMMesh& meshC,const SMAT& P)
  :_meshF(meshF),_meshC(meshC),_P(P) {}
void DHMCoarserGeom::buildParentMap(boost::unordered_map<sizeType,sizeType>& pMap) const
{
  //build coarse vertex-cell map
  CMAP cells;
  for(sizeType c=0; c<_meshC.nrCell(); c++) {
    CELLID cid;
    const DHMCell& cell=_meshC.getCell(c);
    for(char v=0; v<8; v++)
      cid[v]=cell._verts[v]->_index;
    std::sort(cid.data(),cid.data()+8);
    if(cells.find(cid) != cells.end()) {
      INFO("ABC")
      system("pause");
    }
    ASSERT(cells.find(cid) == cells.end())
    cells[cid]=c;
  }

  //build parent map
  for(sizeType c=0; c<_meshF.nrCell(); c++) {
    CELLID cid;
    const DHMCell& cell=_meshF.getCell(c);
    for(char v=0; v<8; v++) {
      sizeType vid=cell._verts[v]->_index;
      if(_P.numElement(vid) == 8) {
        char d=0;
        for(ConstSMIterator<scalarD> beg=_P.begin(vid),end=_P.end(vid); beg!=end; ++beg)
          cid[d++]=beg.col();
        std::sort(cid.data(),cid.data()+8);
        pMap[c]=cells.find(cid)->second;
        break;
      }
    }
  }
}
DHMCoarserGeomAverage::DHMCoarserGeomAverage(DHMMesh& meshF,DHMMesh& meshC,const SMAT& P)
  :DHMCoarserGeom(meshF,meshC,P) {}
void DHMCoarserGeomAverage::coarsen()
{
  //simple coarsening methods, just average the related cells
  SMAT R;
  _P.transpose(R);
  for(sizeType r=0; r<R.rows(); r++) {
    Vec3& pos=_meshC.getVert(r)._pos;
    pos.setZero();

    scalarD total=0.0f;
    for(ConstSMIterator<scalarD> beg=R.begin(r),end=R.end(r); beg!=end; ++beg)
      total+=*beg;
    for(ConstSMIterator<scalarD> beg=R.begin(r),end=R.end(r); beg!=end; ++beg)
      pos+=_meshF.getVert(beg.col())._pos*(scalar)(*beg/total);
  }
}
DHMCoarserGeomInject::DHMCoarserGeomInject(DHMMesh& meshF,DHMMesh& meshC,const SMAT& P)
  :DHMCoarserGeom(meshF,meshC,P) {}
void DHMCoarserGeomInject::coarsen()
{
  vector<bool> visited(_meshC.nrVert(),false);
  for(sizeType r=0; r<_P.rows(); r++)
    if(_P.numElement(r) == 1) {
      sizeType c=_P.begin(r).col();
      visited[c]=true;
      _meshC.getVert(c)._pos=
        _meshF.getVert(r)._pos;
    }

  SMAT R;
  _P.transpose(R);
  for(sizeType r=0; r<R.rows(); r++)
    if(!visited[r]) {
      Vec3& pos=_meshC.getVert(r)._pos;
      pos.setZero();

      scalarD total=0.0f;
      for(ConstSMIterator<scalarD> beg=R.begin(r),end=R.end(r); beg!=end; ++beg)
        total+=*beg;
      for(ConstSMIterator<scalarD> beg=R.begin(r),end=R.end(r); beg!=end; ++beg)
        pos+=_meshF.getVert(beg.col())._pos*(scalar)(*beg/total);
    }
}

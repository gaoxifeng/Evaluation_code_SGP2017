#include "DHMAdaptiveMesh.h"
#include "DHMIter.h"
#include "DHMUtil.h"
#include "CommonFile/IO.h"
#include <deque>
#include <boost/filesystem/fstream.hpp>

USE_PRJ_NAMESPACE

//DHMVertex
DHMVertex::DHMVertex():Serializable(1),_index(numeric_limits<sizeType>::max()) {}
DHMVertex::DHMVertex(sizeType index):Serializable(1),_index(index) {}
bool DHMVertex::read(istream& is,IOData* dat)
{
  readBinaryData(_pos,is);
  readBinaryData(_index,is);
  readBinaryData(_isSurface,is);
  return is.good();
}
bool DHMVertex::write(ostream& os,IOData* dat) const
{
  writeBinaryData(_pos,os);
  writeBinaryData(_index,os);
  writeBinaryData(_isSurface,os);
  return os.good();
}
boost::shared_ptr<Serializable> DHMVertex::copy() const
{
  return boost::shared_ptr<Serializable>(new DHMVertex);
}
DHMVertex::DHMVertex(sizeType type,sizeType index):Serializable(type),_index(index) {}
//DHMCell
DHMCell::DHMCell():Serializable(2),_index(numeric_limits<sizeType>::max()),_layer(numeric_limits<sizeType>::max()) {}
DHMCell::DHMCell(sizeType index):Serializable(2),_index(index),_layer(numeric_limits<sizeType>::max()) {}
void DHMCell::writeVTK(const std::string& path) const
{
  VTKWriter<scalar> os("cell",path,true);
  DHMTraits::POSES pts;
  for(char v=0; v<8; v++)
    pts.push_back(_verts[v]->_pos);
  os.appendPoints(pts.begin(),pts.end());

  DHMTraits::CELLID id;
  DHMTraits::CELLIDS ids;
  id << 0,1,2,3,4,5,6,7;
  ids.push_back(id);
  os.appendCells(ids.begin(),ids.end(),VTKWriter<scalar>::VOXEL);
}
bool DHMCell::read(istream& is,IOData* dat)
{
  for(char v=0; v<8; v++)
    readBinaryData(_verts[v],is,dat);
  readBinaryData(_index,is);
  readBinaryData(_layer,is);
  return is.good();
}
bool DHMCell::write(ostream& os,IOData* dat) const
{
  for(char v=0; v<8; v++)
    writeBinaryData(_verts[v],os,dat);
  writeBinaryData(_index,os);
  writeBinaryData(_layer,os);
  return os.good();
}
boost::shared_ptr<Serializable> DHMCell::copy() const
{
  return boost::shared_ptr<Serializable>(new DHMCell);
}
Vec2i DHMCell::getEdgeId(sizeType i)
{
  char dime=1<<(i/4);
  char dim1=((i%4)/2)<<((i/4+1)%3);
  char dim2=((i%4)%2)<<((i/4+2)%3);
  return Vec2i(dim1+dim2,dim1+dim2+dime);
}
Vec4i DHMCell::getFaceId(sizeType i)
{
  char dimf=(i%2)<<(i/2);
  char dim1=1<<((i/2+1)%3);
  char dim2=1<<((i/2+2)%3);
  if(dim1>dim2)std::swap(dim1,dim2);
  return Vec4i(dimf,dimf+dim1,dimf+dim2,dimf+dim1+dim2);
}
Vec2i DHMCell::getEdge(sizeType i) const
{
  Vec2i ret=getEdgeId(i);
  for(char v=0; v<2; v++)
    ret[v]=_verts[ret[v]]->_index;
  sort2(ret[0],ret[1]);
  return ret;
}
Vec4i DHMCell::getFace(sizeType i) const
{
  Vec4i ret=getFaceId(i);
  for(char v=0; v<4; v++)
    ret[v]=_verts[ret[v]]->_index;
  sort4(ret[0],ret[1],ret[2],ret[3]);
  return ret;
}
Vec4i DHMCell::getFace(const TCOORDS& tcoords,sizeType tid,char d) const
{
  Vec4i ret;
  for(char f=0; f<6; f++) {
    ret=getFace(f);
    if(tcoords[ret[0]][d] == tid && tcoords[ret[1]][d] == tid &&
       tcoords[ret[2]][d] == tid && tcoords[ret[3]][d] == tid)
      return ret;
  }
  ASSERT(false);
  return ret;
}
//DHMFace
DHMFace::DHMFace():Serializable(3) {}
Vec4i DHMFace::getId() const
{
  Vec4i id;
  for(char v=0; v<4; v++)
    id[v]=_verts[v]->_index;
  sort4(id[0],id[1],id[2],id[3]);
  return id;
}
bool DHMFace::read(istream& is,IOData* dat)
{
  readBinaryData(_verts[0],is,dat);
  readBinaryData(_verts[1],is,dat);
  readBinaryData(_verts[2],is,dat);
  readBinaryData(_verts[3],is,dat);
  readBinaryData(_cells[0],is,dat);
  readBinaryData(_cells[1],is,dat);
  return is.good();
}
bool DHMFace::write(ostream& os,IOData* dat) const
{
  writeBinaryData(_verts[0],os,dat);
  writeBinaryData(_verts[1],os,dat);
  writeBinaryData(_verts[2],os,dat);
  writeBinaryData(_verts[3],os,dat);
  writeBinaryData(_cells[0],os,dat);
  writeBinaryData(_cells[1],os,dat);
  return os.good();
}
boost::shared_ptr<Serializable> DHMFace::copy() const
{
  return boost::shared_ptr<Serializable>(new DHMFace);
}
//DHMMesh
//getter
DHMMesh::DHMMesh():Serializable(-1)
{
  _tree.reset(new boost::property_tree::ptree);
}
sizeType DHMMesh::nrVertNC() const
{
  sizeType ret=nrVert()-1;
  while(boost::dynamic_pointer_cast<DHMConstrainedVertex>(_verts[ret]))
    ret--;
  return ret+1;
}
sizeType DHMMesh::nrVert() const
{
  return (sizeType)_verts.size();
}
sizeType DHMMesh::nrCell() const
{
  return (sizeType)_cells.size();
}
sizeType DHMMesh::nrFace() const
{
  return (sizeType)_faces.size();
}
DHMVertex& DHMMesh::getVert(sizeType i)
{
  return *(_verts[i]);
}
const DHMVertex& DHMMesh::getVert(sizeType i) const
{
  return *(_verts[i]);
}
DHMCell& DHMMesh::getCell(sizeType i)
{
  return *(_cells[i]);
}
const DHMCell& DHMMesh::getCell(sizeType i) const
{
  return *(_cells[i]);
}
DHMFace& DHMMesh::getFace(sizeType i)
{
  return *(_faces[i]);
}
const DHMFace& DHMMesh::getFace(sizeType i) const
{
  return *(_faces[i]);
}
const vector<sizeType>& DHMMesh::getAdaptive() const
{
  return _adaptive;
}
void DHMMesh::getPos(COL& pos) const
{
  sizeType nrP=nrVert();
  pos.resize(nrP*3);
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrP; i++)
    pos.block<3,1>(i*3,0)=_verts[i]->_pos.cast<scalarD>();
}
void DHMMesh::setPos(const COL& pos)
{
  sizeType nrP=nrVert();
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<nrP; i++)
    _verts[i]->_pos=pos.block<3,1>(i*3,0).cast<scalar>();
}
//IO
void DHMMesh::assembleVIndex() const
{
  sizeType nrP=nrVert();
  for(sizeType i=0; i<nrP; i++)
    _verts[i]->_index=i;
}
void DHMMesh::reset(const VERTPTRS &verts,const CELLPTRS &cells,const CELLIDS *labels,sizeType nrPadding)
{
  _verts=verts;
  _cells=cells;
  assembleIndex();
  //eliminate hanging vertex
  TBEG("Eliminate Hanging")
  boost::unordered_map<sizeType,boost::shared_ptr<DHMVertex> > validVerts;
  for(sizeType i=0; i<nrCell(); i++)
    for(char v=0; v<8; v++) {
      boost::shared_ptr<DHMVertex> vert=_cells[i]->_verts[v];
      validVerts[vert->_index]=vert;
    }
  _verts.clear();
  for(boost::unordered_map<sizeType,boost::shared_ptr<DHMVertex> >::const_iterator
      beg=validVerts.begin(),end=validVerts.end(); beg!=end; beg++)
    _verts.push_back(beg->second);
  std::sort(_verts.begin(),_verts.end(),LSSVertId());
  assembleIndex();
  TEND
  //build face
  TBEG("Build Face")
  buildFace(labels);
  TEND
  //build layer
  TBEG("Build Layer")
  buildLayer();
  TEND
  TBEG("Build Padding")
  if(nrPadding == numeric_limits<sizeType>::max())
    nrPadding=buildPadding();
  sizeType nrC=nrCell();
  for(sizeType i=0; i<nrC; i++) {
    sizeType& l=_cells[i]->_layer;
    l=std::min<sizeType>(0,l-nrPadding);
  }
  TEND
  //sort according to layer
  TBEG("Ensure RHS")
  assembleIndex();
  //ensureRHS();
  TEND
  INFOV("Mesh has %ld padding layers",getPadding())
}
bool DHMMesh::readVTK(const std::string& path)
{
  Vec3 posV;
  char buf[1024];
  int nrP,nrC,nrI,I[8],type;
  boost::filesystem::ifstream is(path);
  while(is.getline(buf,1024).good()) {
    if(sscanf(buf,"POINTS %d float",&nrP) == 1 ||
       sscanf(buf,"POINTS %d double",&nrP) == 1) {
      _verts.resize(nrP);
      for(sizeType i=0; i<nrP; i++) {
        is >> posV[0] >> posV[1] >> posV[2];
        _verts[i].reset(new DHMVertex(i));
        _verts[i]->_pos=posV;
      }
    } else if(sscanf(buf,"CELLS %d %d",&nrC,&nrI) == 2) {
      _cells.resize(nrC);
      for(sizeType i=0; i<nrC; i++) {
        is.getline(buf,1024);
        sscanf(buf,"8 %d %d %d %d %d %d %d %d",
               I+0,I+1,I+2,I+3,I+4,I+5,I+6,I+7);
        _cells[i].reset(new DHMCell(i));
        for(char v=0; v<8; v++)
          _cells[i]->_verts[v]=_verts[I[v]];
      }
    } else if(sscanf(buf,"CELL_TYPES %d",&nrC) == 1) {
      for(sizeType i=0; i<nrC; i++) {
        is.getline(buf,1024);
        sscanf(buf,"%d",&type);
        if(type == 12) {
          std::swap(_cells[i]->_verts[2],_cells[i]->_verts[3]);
          std::swap(_cells[i]->_verts[6],_cells[i]->_verts[7]);
        }
      }
    } else if(is.eof())ASSERT(false)
    }
  reset(_verts,_cells,NULL);
  return is.good();
}
void DHMMesh::writeVTK(const std::string& path,bool TCoords,const vector<scalar>* data,const vector<scalar>* vdata) const
{
  assembleVIndex();
  VTKWriter<scalar> os("DHMMesh",path,false);
  os.appendPoints(VertPosIter(*this,0),VertPosIter(*this,nrVert()));
  os.appendCells(CellVidIter(*this,0),CellVidIter(*this,nrCell()),VTKWriter<scalar>::VOXEL);
  os.appendCustomData("layer",CellLayerIter(*this,0),CellLayerIter(*this,nrCell()));
  //build TCoords
  if(TCoords && nrVert() == nrVertNC()) {
    //build TCoord
    TCOORDS tcoords;
    buildTCoord(tcoords);
    os.appendCustomPointVectorData("TCoords",tcoords.begin(),tcoords.end());
  }
  //add custom data
  if(data)
    os.appendCustomData("custom",data->begin(),data->end());
  if(vdata)
    os.appendCustomPointData("customP",vdata->begin(),vdata->end());
}
void DHMMesh::writeTVTK(const std::string& path,scalar coef,const vector<scalar>* data) const
{
  assembleVIndex();
  //build TCoord
  TCOORDS tcoords;
  buildTCoord(tcoords);
  POSES verts(tcoords.size());
  for(sizeType i=0; i<(sizeType)tcoords.size(); i++)
    verts[i]=tcoords[i].cast<scalar>()*coef;
  //VTKWriter
  VTKWriter<scalar> os("DHMMesh",path,true);
  os.appendPoints(verts.begin(),verts.end());
  os.appendCells(CellVidIter(*this,0),CellVidIter(*this,nrCell()),VTKWriter<scalar>::VOXEL);
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)tcoords.size(),0,0),
                 VTKWriter<scalar>::POINT);
  os.appendCustomPointVectorData("TCoords",tcoords.begin(),tcoords.end());
  //add duplicity
  {
    boost::unordered_map<Vec3i,sizeType,Hash> dup;
    for(sizeType i=0; i<(sizeType)tcoords.size(); i++)
      dup[tcoords[i]]++;
    vector<sizeType> dups(nrVert());
    for(sizeType i=0; i<nrVert(); i++)
      dups[i]=dup[tcoords[i]];
    os.appendCustomPointData("duplicity",dups.begin(),dups.end());
  }
  //add custom data
  if(data)
    os.appendCustomData("custom",data->begin(),data->end());
}
bool DHMMesh::read(istream& is)
{
  IOData dat;
  dat.registerType<DHMVertex>();
  dat.registerType<DHMConstrainedVertex>();
  dat.registerType<DHMCell>();
  dat.registerType<DHMFace>();
  {
    readVector(_verts,is,&dat);
    readVector(_cells,is,&dat);
    readVector(_faces,is,&dat);
    readVector(_adaptive,is,&dat);
  }
  return is.good();
}
bool DHMMesh::write(ostream& os) const
{
  IOData dat;
  dat.registerType<DHMVertex>();
  dat.registerType<DHMConstrainedVertex>();
  dat.registerType<DHMCell>();
  dat.registerType<DHMFace>();
  {
    writeVector(_verts,os,&dat);
    writeVector(_cells,os,&dat);
    writeVector(_faces,os,&dat);
    writeVector(_adaptive,os,&dat);
  }
  return os.good();
}
void DHMMesh::writePov(const std::string& path) const
{
  boost::unordered_set<Vec2i,Hash> se;
  ObjMesh mesh=getSurface(NULL,&se);

  boost::filesystem::path p(path);
  boost::filesystem::ofstream osPov(p.replace_extension(".pov"));
  mesh.writePov(osPov,false);
  boost::filesystem::ofstream ossPov(p.replace_extension(".spov"));
  mesh.writePov(ossPov,true);

  boost::filesystem::ofstream ose(p.replace_extension(".sedge"));
  ose << se.size() << "," << std::endl;
  for(boost::unordered_set<Vec2i,Hash>::const_iterator
      beg=se.begin(),end=se.end(); beg!=end; beg++) {
    const Vec3& v0=mesh.getV()[(*beg)[0]];
    const Vec3& v1=mesh.getV()[(*beg)[1]];
    ose << v0[0] << "," << v0[1] << "," << v0[2] << ","
        << v1[0] << "," << v1[1] << "," << v1[2] << ",1," << std::endl;
  }
}
void DHMMesh::writeAdaptiveEdge(const BBox<scalar>& bb,sizeType lv,const std::string& path) const
{
  //build subdivision tag
  boost::unordered_set<sizeType> subd;
  for(sizeType i=0; i<nrCell(); i++) {
    const DHMCell& cell=*(_cells[i]);
    for(char v=0; v<8; v++)
      if(bb.contain(cell._verts[v]->_pos)) {
        subd.insert(cell._index);
        break;
      }
  }
  //build face map
  boost::unordered_set<Vec4i,Hash> fids;
  sizeType nrF=nrFace();
  for(sizeType i=0; i<nrF; i++)
    fids.insert(_faces[i]->getId());
  //map surface vert
  vector<Vec4i,Eigen::aligned_allocator<Vec4i> > sf;
  boost::unordered_set<Vec2i,Hash> se;
  for(sizeType i=0; i<nrCell(); i++)
    for(char f=0; f<6; f++) {
      const DHMCell& cell=*(_cells[i]);
      Vec4i cf=cell.getFace(f);
      Vec4i cfid=cell.getFaceId(f);
      if(fids.find(cf) == fids.end())
        if(subd.find(cell._index) != subd.end()) {
          sf.push_back(Vec4i(cell._verts[cfid[0]]->_index,
                             cell._verts[cfid[1]]->_index,
                             cell._verts[cfid[2]]->_index,
                             cell._verts[cfid[3]]->_index));
        } else {
#define ADD_SEDGE(A,B)	\
eid=Vec2i(cell._verts[cfid[A]]->_index,cell._verts[cfid[B]]->_index);	\
sort2(eid[0],eid[1]);	\
se.insert(eid);
          Vec2i eid;
          ADD_SEDGE(0,1)
          ADD_SEDGE(0,2)
          ADD_SEDGE(3,1)
          ADD_SEDGE(3,2)
        }
#undef ADD_SEDGE
    }

  //write edge
  POSES epos;
  sizeType nr=1<<lv;
  for(sizeType i=0; i<(sizeType)sf.size(); i++) {
    Vec3 pos[4]= {_verts[sf[i][0]]->_pos,
                  _verts[sf[i][1]]->_pos,
                  _verts[sf[i][2]]->_pos,
                  _verts[sf[i][3]]->_pos
                 };
#define DEL_SEDGE(A,B)	\
eid=Vec2i(sf[i][A],sf[i][B]);	\
sort2(eid[0],eid[1]);	\
se.erase(eid);
    Vec2i eid;
    DEL_SEDGE(0,1)
    DEL_SEDGE(0,2)
    DEL_SEDGE(3,1)
    DEL_SEDGE(3,2)
#undef DEL_SEDGE
    for(sizeType x=0; x<=nr; x++)
      for(sizeType y=0; y<nr; y++) {
        epos.push_back(interp2D(pos[0],pos[1],pos[2],pos[3],(scalar)x/(scalar)nr,(scalar)y/(scalar)nr));
        epos.push_back(interp2D(pos[0],pos[1],pos[2],pos[3],(scalar)x/(scalar)nr,(scalar)(y+1)/(scalar)nr));
      }
    for(sizeType x=0; x<nr; x++)
      for(sizeType y=0; y<=nr; y++) {
        epos.push_back(interp2D(pos[0],pos[1],pos[2],pos[3],(scalar)x/(scalar)nr,(scalar)y/(scalar)nr));
        epos.push_back(interp2D(pos[0],pos[1],pos[2],pos[3],(scalar)(x+1)/(scalar)nr,(scalar)y/(scalar)nr));
      }
  }
  sizeType nrSE=epos.size();
  for(boost::unordered_set<Vec2i,Hash>::const_iterator
      beg=se.begin(),end=se.end(); beg!=end; beg++) {
    epos.push_back(_verts[(*beg)[0]]->_pos);
    epos.push_back(_verts[(*beg)[1]]->_pos);
  }

  //write edge
  boost::filesystem::ofstream ose(path);
  ose << epos.size()/2 << "," << std::endl;
  for(sizeType i=0; i<(sizeType)epos.size(); i+=2) {
    ose << epos[i][0] << "," << epos[i][1] << "," << epos[i][2] << ","
        << epos[i+1][0] << "," << epos[i+1][1] << "," << epos[i+1][2] << "," << (i<nrSE?0:1) << "," << std::endl;
  }
}
ObjMesh DHMMesh::getSurface(boost::unordered_map<sizeType,sizeType>* vertMap,boost::unordered_set<Vec2i,Hash>* edgeMap) const
{
  //build face map
  boost::unordered_set<Vec4i,Hash> fids;
  sizeType nrF=nrFace();
  for(sizeType i=0; i<nrF; i++)
    fids.insert(_faces[i]->getId());
  //map surface vert
  ObjMesh ret;
  sizeType nrV=0;
  boost::unordered_map<sizeType,sizeType> vmap;
  for(sizeType i=0; i<nrCell(); i++)
    for(char f=0; f<6; f++) {
      const DHMCell& cell=*(_cells[i]);
      Vec4i cf=cell.getFace(f);
      if(fids.find(cf) == fids.end()) {
        //map vertex
        for(char v=0; v<4; v++)
          if(vmap.find(_verts[cf[v]]->_index) == vmap.end()) {
            ret.getV().push_back(_verts[cf[v]]->_pos);
            vmap[cf[v]]=nrV++;
          }

        //map index
        Vec4i fid=cell.getFaceId(f);
#define VID(I) vmap[cell._verts[fid[I]]->_index]
        ret.getI().push_back(Vec3i(VID(0),VID(1),VID(3)));
        ret.getI().push_back(Vec3i(VID(0),VID(3),VID(2)));
        if(edgeMap) {
          Vec2i eid;
#define ADD_SEDGE(A,B)	\
eid=Vec2i(VID(A),VID(B));	\
sort2(eid[0],eid[1]);	\
edgeMap->insert(eid);
          ADD_SEDGE(0,1)
          ADD_SEDGE(0,2)
          ADD_SEDGE(3,1)
          ADD_SEDGE(3,2)
#undef ADD_SEDGE
        }
#undef VID
      }
    }
  if(vertMap)
    *vertMap=vmap;
  ret.smooth();
  ret.makeUniform();
  ret.smooth();
  return ret;
}
//add padding by subdivision
sizeType DHMMesh::getPadding() const
{
  return -_cells[0]->_layer;
}
void DHMMesh::addPadding(sizeType nrL,bool postResample)
{
  //subdivide
  boost::unordered_set<Vec4i,Hash> fids;
  sizeType nrF=nrFace();
  for(sizeType i=0; i<nrF; i++)
    fids.insert(_faces[i]->getId());
  //find surface cell
  PMAP padMap;
  CELLPTRS newCells;
  sizeType nrPadding=getPadding();
  for(sizeType i=0; i<nrCell(); i++)
    for(char f=0; f<6; f++) {
      Vec4i cf=_cells[i]->getFace(f);
      if(fids.find(cf) == fids.end()) {
        //subdivide
        sizeType nrNew=(sizeType)newCells.size();
        subdivide(*(_cells[i]),cf,nrL+1,padMap,newCells);
        //assign layer for new paddings
        for(sizeType l=0; l<=nrL; l++)
          newCells[nrNew+l]->_layer=-nrL-nrPadding+l;
        _cells[i]->_layer=numeric_limits<sizeType>::max();
        break;
      }
    }
  //erase old cells
  deleteRedundant();
  //insert new cells
  _cells.insert(_cells.end(),newCells.begin(),newCells.end());
  //reassemble index
  assembleIndex();
  buildFace();
  ensureRHS();
  //resample
  if(postResample) {
    TCOORDS tcoords;
    buildTCoord(tcoords);
    resample(tcoords,-1);
  }
}
void DHMMesh::subdivide()
{
  //double padding
  addPadding(getPadding(),true);

  //build TCOORDS
  TCOORDS tcoords;
  buildTCoord(tcoords);

  //find min/max TCOORDS
  BBox<sizeType> bb;
  for(sizeType i=0; i<(sizeType)tcoords.size(); i++)
    bb.setUnion(tcoords[i]);
  for(char d=0; d<3; d++)
    for(sizeType t=bb._minC[d],tt=t; t<bb._maxC[d]; t++,tt+=2)
      addLayer(tcoords,tt,d,1);

  //reassemble
  assembleIndex();
  buildFace();
  ensureRHS();
}
//add layers to ensure multgrid coarsening safety
void DHMMesh::patchPOT(bool postResample)
{
  //finding minimal padding boundary
  sizeType nrPadding=getPadding();
  if(!nrPadding)
    return;
  sizeType pot=1;
  for(; pot < nrPadding; pot<<=1);
  INFOV("Need to align to %ld boundary!",pot)
  if(pot > nrPadding)
    addPadding(pot-nrPadding);
  //build singularity
  POSSET svert;
  buildPadding(&svert);
  //build tcoords
  TCOORDS tcoords;
  buildTCoord(tcoords);
  //sort singularity along XYZ
  vector<sizeType> tcoordsD;
  for(CPOSITER beg=svert.begin(),end=svert.end(); beg!=end; beg++)
    tcoordsD.push_back(*beg);
  for(char d=0; d<3; d++) {
    //find different of two coords, force them to be POT by subdivision
    LSSTCoordD lss(tcoords,d);
    std::sort(tcoordsD.begin(),tcoordsD.end(),lss);
    for(sizeType i=1; i<(sizeType)tcoordsD.size(); i++) {
      sizeType t1=lss(tcoordsD[i]);
      sizeType t0=lss(tcoordsD[i-1]);
      if((t1-t0)%pot != 0) {
        sizeType nrL=pot-(t1-t0)%pot;
        addLayer(tcoords,t0,d,nrL);
        if(postResample)
          resample(tcoords,d,t0,t1+nrL);
        //check still sorted
        for(sizeType t=1; t<(sizeType)tcoordsD.size(); t++)
          ASSERT(!lss(tcoordsD[t],tcoordsD[t-1]));
      }
    }
  }

  //reassemble
  assembleIndex();
  buildFace();
  ensureRHS();
}
void assignTCoord(const DHMCell& cell,const DHMCell& cellO,const Vec4i& fid,DHMTraits::TCOORDS& tcoords)
{
  Vec3i delta;
  for(char d=0; d<3; d++) {
    if(cell.getFace(d*2+0) == fid) {
      delta=-Vec3i::Unit(d);
      break;
    } else if(cell.getFace(d*2+1) == fid) {
      delta=Vec3i::Unit(d);
      break;
    }
  }
  Vec4i oppV0;
  Vec4i oppV1;
  for(char d=0; d<3; d++) {
    if(cellO.getFace(d*2+0) == fid) {
      oppV0=DHMCell::getFaceId(d*2+0);
      oppV1=DHMCell::getFaceId(d*2+1);
      break;
    } else if(cellO.getFace(d*2+1) == fid) {
      oppV0=DHMCell::getFaceId(d*2+1);
      oppV1=DHMCell::getFaceId(d*2+0);
      break;
    }
  }
  for(char v=0; v<4; v++)
    tcoords[cellO._verts[oppV1[v]]->_index]=
      tcoords[cellO._verts[oppV0[v]]->_index]+delta;
}
void copyTCoord(const DHMCell& cellO,const Vec4i& fid,DHMTraits::TCOORDS& tcoords)
{
  Vec4i oppV0;
  Vec4i oppV1;
  for(char d=0; d<3; d++) {
    if(cellO.getFace(d*2+0) == fid) {
      oppV0=DHMCell::getFaceId(d*2+0);
      oppV1=DHMCell::getFaceId(d*2+1);
      break;
    } else if(cellO.getFace(d*2+1) == fid) {
      oppV0=DHMCell::getFaceId(d*2+1);
      oppV1=DHMCell::getFaceId(d*2+0);
      break;
    }
  }
  for(char v=0; v<4; v++)
    tcoords[cellO._verts[oppV1[v]]->_index]=
      tcoords[cellO._verts[oppV0[v]]->_index];
}
void DHMMesh::buildTCoord(TCOORDS& tcoords) const
{
  sizeType nrP=nrVert();
  sizeType nrC=nrCell();
  tcoords.assign(nrP,Vec3i::Constant(numeric_limits<sizeType>::max()));

  //search for an inner cell
  std::deque<sizeType> que;
  vector<bool> visited(nrP);
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& c=*(_cells[i]);
    if(c._layer >= 0) {
      visited[i]=true;
      que.push_back(i);
      for(char v=0; v<8; v++)
        tcoords[c._verts[v]->_index]=Vec3i(v&1?1:0,v&2?1:0,v&4?1:0);
      break;
    }
  }

  //fast march other cells
  FNMAP faceMap;
  faceMap.reserve(nrFace()*2);
  for(sizeType i=0; i<nrFace(); i++) {
    const DHMFace& face=*(_faces[i]);
    faceMap.insert(face._cells[0]->_index,_faces[i]);
    faceMap.insert(face._cells[1]->_index,_faces[i]);
  }
  FNMAP::SET neighF;
  while(!que.empty()) {
    const DHMCell& c=*(_cells[que.front()]);
    que.pop_front();
    faceMap.get(c._index,neighF);
    for(FNMAP::SET::const_iterator beg=neighF.begin(),end=neighF.end(); beg!=end; beg++) {
      const DHMFace& face=**beg;
      if(!face._cells[1])
        continue;
      const DHMCell& co=(face._cells[0]->_index == c._index) ? *(face._cells[1]) : *(face._cells[0]);
      if(co._layer >= 0 && !visited[co._index]) {
        visited[co._index]=1;
        que.push_back(co._index);
        assignTCoord(c,co,face.getId(),tcoords);
      }
    }
  }

  //copy to padded cells
  sizeType nrPadding=getPadding();
  for(sizeType l=0; l<nrPadding; l++)
    for(sizeType i=0; i<nrFace(); i++)
      if(_faces[i]->_cells[1]) {
        const DHMCell& c0=*(_faces[i]->_cells[0]);
        const DHMCell& c1=*(_faces[i]->_cells[1]);
        if(c0._layer==-l && c1._layer==-l-1)
          copyTCoord(c1,_faces[i]->getId(),tcoords);
        else if(c0._layer==-l-1 && c1._layer==-l)
          copyTCoord(c0,_faces[i]->getId(),tcoords);
      }

  //force each tcoord to be positive
  Vec3i tc=Vec3i::Constant(numeric_limits<sizeType>::max());
  for(sizeType i=0; i<nrP; i++) {
    ASSERT(tcoords[i][0] != numeric_limits<sizeType>::max())
    ASSERT(tcoords[i][1] != numeric_limits<sizeType>::max())
    ASSERT(tcoords[i][2] != numeric_limits<sizeType>::max())
    tc=compMin(tc,tcoords[i]);
  }
  for(sizeType i=0; i<nrP; i++)
    tcoords[i]-=tc;
}
//assembler
void DHMMesh::ensureRHS()
{
  //build faceMap and tcoords
  TCOORDS tcoords;
  buildTCoord(tcoords);

  //build cindices
  sizeType nrC=nrCell();
  sizeType nrV=nrVert();
  sizeType nrF=nrFace();
  vector<sizeType> vLayers(nrV,getPadding());
  vector<Vec3i,Eigen::aligned_allocator<Vec3i> > cindices(nrC,Vec3i::Zero());
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& cell=*(_cells[i]);
    for(char v=0; v<8; v++) {
      sizeType vid=cell._verts[v]->_index;
      cindices[i]+=tcoords[vid];
      vLayers[vid]=std::min<sizeType>(vLayers[vid],-cell._layer);
    }
  }

  //build cell type
  vector<sizeType> cType(nrC,-1);
  for(sizeType l=0; l>-getPadding(); l--)
    for(sizeType i=0; i<nrF; i++) {
      const DHMFace& face=*(_faces[i]);
      sizeType cs0=face._cells[0]->_index;
      sizeType cs1=face._cells[1]->_index;
      sizeType l0=face._cells[0]->_layer;
      sizeType l1=face._cells[1]->_layer;
      if(l1 > l0) {
        std::swap(l0,l1);
        std::swap(cs0,cs1);
      }
      if(l0 == l && l1 == l-1)
        if(l0 == 0) {
          for(char d=0; d<3; d++)
            if(cindices[cs0][d] > cindices[cs1][d]) {
              cType[cs1]=d*2;
              break;
            } else if(cindices[cs0][d] < cindices[cs1][d]) {
              cType[cs1]=d*2+1;
              break;
            }
        } else cType[cs1]=cType[cs0];
    }

  //build coordinate frame
  for(sizeType i=0; i<nrC; i++) {
    Vec3i crd;
    switch(cType[i]) {
    case 0:
      crd=Vec3i(1,3,2);
      break;//-X
    case 1:
      crd=Vec3i(3,1,2);
      break;//+X
    case 2:
      crd=Vec3i(3,0,2);
      break;//-Y
    case 3:
      crd=Vec3i(0,3,2);
      break;//+Y
    case 4:
      crd=Vec3i(1,0,3);
      break;//-Z
    case 5:
      crd=Vec3i(0,1,3);
      break;//+Z
    default:
      crd=Vec3i(0,1,2);
      break;
    }

    Vec4i evid[8],evidMin=Vec4i::Constant(numeric_limits<sizeType>::max());
    boost::shared_ptr<DHMVertex> vNew[8];
    bool debug[8];
    DHMCell& cell=*(_cells[i]);
    for(char v=0; v<8; v++) {
      sizeType vid=cell._verts[v]->_index;
      evid[v]=Vec4i(tcoords[vid][0],tcoords[vid][1],tcoords[vid][2],vLayers[vid]);
      evidMin=compMin(evidMin,evid[v]);
      debug[v]=false;
    }
    for(char v=0; v<8; v++) {
      Vec4i dt=evid[v]-evidMin;
      sizeType off=Vec3i(dt[crd[0]],dt[crd[1]],dt[crd[2]]).dot(Vec3i(1,2,4));
      ASSERT(off >= 0 && off < 8)
      vNew[off]=cell._verts[v];
      debug[off]=true;
    }
    for(char v=0; v<8; v++) {
      ASSERT(debug[v])
      cell._verts[v]=vNew[v];
    }
  }
}
void DHMMesh::assembleIndex()
{
  assembleVIndex();
  std::sort(_cells.begin(),_cells.end(),LSSCellLayer());
  sizeType nrC=nrCell();
  for(sizeType i=0; i<nrC; i++)
    _cells[i]->_index=i;
}
void DHMMesh::buildFace(const CELLIDS *labels)
{
  //build face map
  typedef multi_unordered_map<Vec4i,sizeType,true> FMAP;
  typedef FMAP::const_iterator CFITER;
  FMAP faceMap;
  sizeType nrC=nrCell();
  faceMap.reserve(nrC*6);
  for(sizeType i=0; i<nrC; i++)
    for(char f=0; f<6; f++)
      faceMap.insert(_cells[i]->getFace(f),i);

  _faces.clear();
  //assemble faces
  FMAP::SET vss;
  bool nonManifold=false;
  for(CFITER beg=faceMap.begin(),end=faceMap.end(); beg!=end; beg++) {
    //build neighbor cells
    faceMap.get(beg->first,vss);
    if(vss.size() > 2)
      nonManifold=true;

    if(vss.size() == 1)
    {
      _verts[beg->first[0]]->_isSurface=true;
      _verts[beg->first[1]]->_isSurface=true;
      _verts[beg->first[2]]->_isSurface=true;
      _verts[beg->first[3]]->_isSurface=true;
    }

    //add faces
    sizeType nrV=(sizeType)vss.size();
    for(sizeType i=0; i<nrV; i++)
      for(sizeType j=i+1; j<nrV; j++)
        if(!labels || (*labels)[vss[i]] != (*labels)[vss[j]]) {
          boost::shared_ptr<DHMFace> face(new DHMFace);
          for(char v=0; v<4; v++)
            face->_verts[v]=_verts[beg->first[v]];
          face->_cells[0]=_cells[vss[i]];
          face->_cells[1]=_cells[vss[j]];
          _faces.push_back(face);
        }
  }
  if(nonManifold) {
    ASSERT(labels)
  }
}
void DHMMesh::buildLayer()
{
  //initialize
  sizeType nrC=nrCell();
  for(sizeType i=0; i<nrC; i++)
    _cells[i]->_layer=numeric_limits<sizeType>::max();

  //build outtermost cells
  FNMAP faceMap;
  FIDSET fids;
  std::deque<sizeType> que;
  sizeType nrF=nrFace();
  faceMap.reserve(nrF*2);
  fids.reserve(nrF);
  for(sizeType i=0; i<nrF; i++) {
    const DHMFace& face=*(_faces[i]);
    faceMap.insert(face._cells[0]->_index,_faces[i]);
    faceMap.insert(face._cells[1]->_index,_faces[i]);
    fids.insert(face.getId());
  }
  for(sizeType i=0; i<nrC; i++)
    for(char f=0; f<6; f++) {
      Vec4i cf=_cells[i]->getFace(f);
      if(fids.find(cf) == fids.end()) {
        _cells[i]->_layer=0;
        que.push_back(_cells[i]->_index);
        break;
      }
    }

  //fast march inner cells
  FNMAP::SET neighF;
  for(sizeType l=1; !que.empty(); l++)
    for(sizeType q=0,nrQ=(sizeType)que.size(); q<nrQ; q++) {
      const DHMCell& c=*(_cells[que.front()]);
      que.pop_front();
      faceMap.get(c._index,neighF);
      for(FNMAP::SET::const_iterator beg=neighF.begin(),end=neighF.end(); beg!=end; beg++) {
        const DHMFace& face=**beg;
        if(!face._cells[1])
          continue;
        DHMCell& co=(face._cells[0]->_index == c._index) ? *(face._cells[1]) : *(face._cells[0]);
        if(co._layer == numeric_limits<sizeType>::max()) {
          co._layer=l;
          que.push_back(co._index);
        }
      }
    }
}
sizeType DHMMesh::buildPadding(POSSET* svert) const
{
  typedef boost::fast_pool_allocator<std::pair<const Vec2i,Vec2i> > ALLOC_EMAP;
  typedef boost::unordered_map<Vec2i,Vec2i,Hash,std::equal_to<Vec2i>,ALLOC_EMAP> EMAP;
  typedef EMAP::const_iterator CEITER;
  typedef EMAP::iterator EITER;
  //build edge info map
  //if edgeMap::second[1]=numeric_limits<sizeType>::max(), this means edge cross layers
  //we don't consider such edge
  Vec2i edge;
  EMAP edgeMap;
  sizeType nrC=nrCell();
  edgeMap.reserve(nrC*12);
  for(sizeType i=0; i<nrC; i++) {
    const DHMCell& c=*(_cells[i]);
    for(char e=0; e<12; e++) {
      edge=c.getEdge(e);
      EITER it=edgeMap.find(edge);
      if(it == edgeMap.end())
        edgeMap[edge]=Vec2i(1,c._layer);
      else if(it->second[1] != numeric_limits<sizeType>::max()) {
        if(it->second[1] == c._layer)
          it->second[0]++;
        else it->second[1]=numeric_limits<sizeType>::max();
      }
    }
  }

  //exclude surface edge
  FIDSET fids;
  sizeType nrF=nrFace();
  fids.reserve(nrF);
  for(sizeType i=0; i<nrF; i++)
    fids.insert(_faces[i]->getId());
  for(sizeType i=0; i<nrC; i++)
    for(char f=0; f<6; f++) {
      Vec4i cf=_cells[i]->getFace(f);
      if(fids.find(cf) == fids.end())
        for(char a=0; a<4; a++)
          for(char b=a; b<4; b++) {
            Vec2i eid(cf[a],cf[b]);
            sort2(cf[0],cf[1]);
            EITER it=edgeMap.find(eid);
            if(it != edgeMap.end())
              it->second[1]=numeric_limits<sizeType>::max();
          }
    }

  //build padding according to singularity edge
  sizeType padding=0;
  for(CEITER beg=edgeMap.begin(),end=edgeMap.end(); beg!=end; beg++)
    if(beg->second[0] != 4 && beg->second[1] != numeric_limits<sizeType>::max()) {
      if(svert) {
        svert->insert(beg->first[0]);
        svert->insert(beg->first[1]);
      }
      padding=std::max(padding,beg->second[1]+1);
    }
  return padding;
}
//subdivider
Vec3i getMinTC(const DHMCell& c,const DHMCell::TCOORDS& tc)
{
  Vec3i ret=tc[c._verts[0]->_index];
  for(char v=1; v<8; v++)
    ret=compMin(ret,tc[c._verts[v]->_index]);
  return ret;
}
Vec3i getMaxTC(const DHMCell& c,const DHMCell::TCOORDS& tc)
{
  Vec3i ret=tc[c._verts[0]->_index];
  for(char v=1; v<8; v++)
    ret=compMax(ret,tc[c._verts[v]->_index]);
  return ret;
}
void DHMMesh::subdivide(const DHMCell& cell,const Vec4i& face,sizeType nrSlice,PMAP& padMap,CELLPTRS& cells)
{
  //subdivide [cell] along [face] into [nrSlice] cells
  //put new vertices into [padMap]
  //put new cells into [cells]
#define GETV(VID)_verts[VID]->_pos
  Vec2i sedge[4],sedgeu;//find 4 edges connecting the opposite face
  vector<sizeType> svert[4];//subdivided vertex offset
  char sedgeMap[8];//map each vertex to an edge
  for(char d=0; d<3; d++) {
    bool isD0=face == cell.getFace(d*2);
    bool isD1=face == cell.getFace(d*2+1);
    if(isD0 || isD1) {
      for(char i=0; i<4; i++) {
        //build map
        sedge[i]=DHMCell::getEdgeId(d*4+i);
        sedgeMap[sedge[i][0]]=i;
        sedgeMap[sedge[i][1]]=i;
        //build edge tag
        sedge[i][0]=cell._verts[sedge[i][0]]->_index;
        sedge[i][1]=cell._verts[sedge[i][1]]->_index;
        if(isD1)std::swap(sedge[i][0],sedge[i][1]);
        //insert subdivide vertices
        sedgeu=sedge[i];	//build a unique index
        sort2(sedgeu[0],sedgeu[1]);
        CPITER it=padMap.find(sedgeu);
        if(it == padMap.end()) {
          svert[i].push_back(sedge[i][0]);
          for(sizeType v=1; v<nrSlice; v++) {
            boost::shared_ptr<DHMVertex> vert(new DHMVertex);
            vert->_pos=(GETV(sedge[i][0])*(scalar)(nrSlice-v)+
                        GETV(sedge[i][1])*(scalar)(v))/(scalar)nrSlice;
            vert->_index=nrVert();
            svert[i].push_back(vert->_index);
            _verts.push_back(vert);
          }
          svert[i].push_back(sedge[i][1]);
          padMap.set(sedgeu,svert[i].begin(),svert[i].end());
        } else padMap.get(it->first,svert[i]);
      }
      break;
    }
  }
#undef GETV
  //now subdivide
  for(sizeType i=0; i<nrSlice; i++) {
    boost::shared_ptr<DHMCell> c(new DHMCell);
    for(char v=0; v<8; v++) {
      const Vec2i& sedgev=sedge[sedgeMap[v]];
      const vector<sizeType>& svertv=svert[sedgeMap[v]];
      if(cell._verts[v]->_index == sedgev[0]) {
        c->_verts[v]=_verts[svertv[i]];
      } else {
        ASSERT(cell._verts[v]->_index == sedgev[1])
        c->_verts[v]=_verts[svertv[i+1]];
      }
    }
    cells.push_back(c);
  }
}
void DHMMesh::addLayer(TCOORDS& tcoords,sizeType tid,char d,sizeType nrL)
{
  //increase tcoords
  for(sizeType i=0; i<(sizeType)tcoords.size(); i++)
    if(tcoords[i][d] > tid)
      tcoords[i][d]+=nrL;
  //find cells crossing multiple tcoords and subdivide them
  PMAP padMap;
  CELLPTRS newCells;
  sizeType nrC=nrCell();
  for(sizeType i=0; i<nrC; i++) {
    DHMCell& c=*(_cells[i]);
    sizeType minD=getMinTC(c,tcoords)[d];
    sizeType maxD=getMaxTC(c,tcoords)[d];
    if(maxD == minD+nrL+1) {
      sizeType nrNew=(sizeType)newCells.size();
      subdivide(c,c.getFace(tcoords,minD,d),nrL+1,padMap,newCells);
      for(; nrNew<(sizeType)newCells.size(); nrNew++)
        newCells[nrNew]->_layer=c._layer;
      c._layer=numeric_limits<sizeType>::max();
    }
  }
  //assign new tcoords
  tcoords.resize(nrVert());
  vector<sizeType> vss;
  for(CPITER beg=padMap.begin(),end=padMap.end(); beg!=end; beg++) {
    padMap.get(beg->first,vss);
    for(sizeType i=1; i<=nrL; i++)
      tcoords[vss[i]]=tcoords[vss[0]]+Vec3i::Unit(d)*i;
    ASSERT(tcoords[vss[0]]+Vec3i::Unit(d)*(nrL+1) == tcoords[vss[nrL+1]])
  }
  //delete redundent
  deleteRedundant();
  //insert new cells
  _cells.insert(_cells.end(),newCells.begin(),newCells.end());
  //reassemble index (deprecated for internal operation)
  assembleIndex();
}
void DHMMesh::deleteRedundant()
{
  for(sizeType i=0; i<nrCell();)
    if(_cells[i]->_layer == numeric_limits<sizeType>::max()) {
      _cells[i]=_cells.back();
      _cells.pop_back();
    } else i++;
}
//resample vertices according to arc-length
void DHMMesh::resample(const TCOORDS& tcoords,char d,sizeType f,sizeType t)
{
  //find all edges satisfying condition
  typedef multi_unordered_map<sizeType,sizeType,true> CONNMAP;
  typedef CONNMAP::const_iterator CCONNITER;
  CONNMAP connMap;
  sizeType nrC=nrCell();
  for(sizeType i=0; i<nrC; i++)
    for(char e=0; e<12; e++) {
      const Vec2i eid=_cells[i]->getEdge(e);
      const Vec3i& t0=tcoords[eid[0]];
      const Vec3i& t1=tcoords[eid[1]];
      if(d == -1) {
        if(t0 == t1) {
          connMap.insert(eid[0],eid[1]);
          connMap.insert(eid[1],eid[0]);
        }
      } else if(t0[d] != t1[d] && t0[d] >= f && t0[d] <= t && t1[d] >= f && t1[d] <= t) {
        connMap.insert(eid[0],eid[1]);
        connMap.insert(eid[1],eid[0]);
      }
    }

  //resample each edge
  typedef boost::fast_pool_allocator<sizeType> ALLOC_NPMAP;
  typedef boost::unordered_map<sizeType,Vec3,Hash,std::equal_to<sizeType>,ALLOC_NPMAP> NPMAP;
  typedef NPMAP::const_iterator CNPITER;
  NPMAP newPos;
  CONNMAP::SET conn;
  vector<bool> visited(nrVert(),false);
  for(CCONNITER beg=connMap.begin(),end=connMap.end(); beg!=end; beg++) {
    connMap.get(beg->first,conn);
    sizeType nrConn=(sizeType)conn.size();
    ASSERT(nrConn == 1 || nrConn == 2)
    if(nrConn == 1 && !visited[beg->first]) {
      //initialize
      vector<sizeType> que;
      visited[beg->first]=true;
      que.push_back(beg->first);
      //search
      bool finish=false;
      while(!finish) {
        connMap.get(que.back(),conn);
        finish=true;
        for(vector<sizeType>::const_iterator beg=conn.begin(),end=conn.end(); beg!=end; beg++)
          if(!visited[*beg]) {
            visited[*beg]=true;
            que.push_back(*beg);
            finish=false;
            break;
          }
      }
      //resample
      scalar len=0.0f;
      for(sizeType i=1; i<(sizeType)que.size(); i++)
        len+=(_verts[que[i-1]]->_pos-_verts[que[i]]->_pos).norm();
      for(sizeType i=1; i<(sizeType)que.size()-1; i++)
        newPos[que[i]]=resample(que,len*(scalar)i/(scalar)(que.size()-1));
    }
  }

  //assign new pos
  for(CNPITER beg=newPos.begin(),end=newPos.end(); beg!=end; beg++)
    _verts[beg->first]->_pos=beg->second;
}
Vec3 DHMMesh::resample(const vector<sizeType>& que,scalar len) const
{
  scalar curr=0.0f;
  for(sizeType i=1; i<(sizeType)que.size(); i++) {
    const Vec3& p0=_verts[que[i-1]]->_pos;
    const Vec3& p1=_verts[que[i]]->_pos;
    scalar lenSeg=(p0-p1).norm();
    if(curr+lenSeg < len)
      curr+=lenSeg;
    else return (p0*(lenSeg+curr-len)+p1*(len-curr))/lenSeg;
  }

  ASSERT(false)
  return Vec3::Zero();
}

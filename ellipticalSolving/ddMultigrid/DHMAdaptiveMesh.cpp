#include "DHMAdaptiveMesh.h"
#include "DHMCoarser.h"
#include "DHMUtil.h"
#include "DHMIter.h"

USE_PRJ_NAMESPACE

//DHMConstrainedVertex
DHMConstrainedVertex::DHMConstrainedVertex():DHMVertex(4,numeric_limits<sizeType>::max()) {}
DHMConstrainedVertex::DHMConstrainedVertex(sizeType index):DHMVertex(4,index) {}
bool DHMConstrainedVertex::read(istream& is,IOData* dat)
{
  DHMVertex::read(is,dat);
  readBinaryData(_verts[0],is,dat);
  readBinaryData(_verts[1],is,dat);
  readBinaryData(_verts[2],is,dat);
  readBinaryData(_verts[3],is,dat);
  readBinaryData(_weights,is);
  return is.good();
}
bool DHMConstrainedVertex::write(ostream& os,IOData* dat) const
{
  DHMVertex::write(os,dat);
  writeBinaryData(_verts[0],os,dat);
  writeBinaryData(_verts[1],os,dat);
  writeBinaryData(_verts[2],os,dat);
  writeBinaryData(_verts[3],os,dat);
  writeBinaryData(_weights,os);
  return os.good();
}
boost::shared_ptr<Serializable> DHMConstrainedVertex::copy() const
{
  return boost::shared_ptr<Serializable>(new DHMConstrainedVertex);
}
//DHMAdaptiveMesh
DHMAdaptiveMesh::DHMAdaptiveMesh(Hierarchy& hier):_hier(hier) {}
const DHMMesh& DHMAdaptiveMesh::getBaseMesh() const
{
  for(sizeType i=0;; i++)
    if(i >= (sizeType)_hier._prolong.size() || _hier._prolong[i].rows() == _hier._mesh[i]->nrVert())
      return *(_hier._mesh[i]);
  return *(DHMMesh*)NULL;
}
void DHMAdaptiveMesh::writeVTKLevel(const std::string& path,sizeType lv)
{
  //reorder
  DHMMesh& mesh=*(_hier._mesh[lv]);
  mesh.assembleVIndex();
  VERTPTRS& verts=mesh._verts;
  CELLPTRS& cells=mesh._cells;
  //write
  VTKWriter<scalar> os("DHMMesh",path,false);
  os.appendPoints(VertPosIter(mesh,0),VertPosIter(mesh,mesh.nrVert()));
  os.appendCells(CellVidIter(mesh,0),CellVidIter(mesh,mesh.nrCell()),VTKWriter<scalar>::VOXEL);
  os.appendCustomData("layer",CellLayerIter(mesh,0),CellLayerIter(mesh,mesh.nrCell()));
  os.appendCustomData("index",CellIdIter(mesh,0),CellIdIter(mesh,mesh.nrCell()));
  os.appendCustomPointData("cons",VertTypeIter(mesh,0,getBaseMesh().nrVert()),VertTypeIter(mesh,mesh.nrVert()));
}
void DHMAdaptiveMesh::clearAdaptiveMesh()
{
  while(!_hier._prolong.empty() && _hier._prolong.front().rows() < _hier._mesh.front()->nrVert()) {
    _hier._mesh.erase(_hier._mesh.begin());
    _hier._prolong.erase(_hier._prolong.begin());
  }
}
void DHMAdaptiveMesh::createAdaptiveMesh(const BBox<scalar>& bb,sizeType lv)
{
  clearAdaptiveMesh();
  boost::shared_ptr<DHMMesh> baseMesh=_hier._mesh.front();
  vector<sizeType> cells;
  for(sizeType i=0; i<baseMesh->nrCell(); i++) {
    for(char v=0; v<8; v++)
      if(bb.contain(baseMesh->getCell(i)._verts[v]->_pos)) {
        cells.push_back(i);
        break;
      }
  }
  createAdaptiveMesh(cells,lv);
}
void addPRow(sizeType id,const Vec3i& off,boost::shared_ptr<DHMVertex> vertP[8],DHMEnergyTraits<scalarD>::TRIPS& trips)
{
  ASSERT(compGE(off,Vec3i::Zero()) && compLE(off,Vec3i::Constant(2)));
  Eigen::Matrix<scalar,8,1> Ms;
  DHMMapping::MStencil(Ms,(off-Vec3i::Ones()).cast<scalar>());
  for(char v=0; v<8; v++)
    if(Ms[v] > 0.0f) {
      boost::shared_ptr<DHMConstrainedVertex> vert=
        boost::dynamic_pointer_cast<DHMConstrainedVertex>(vertP[v]);
      if(vert) {
        for(char cv=0; cv<4; cv++)
          if(vert->_verts[cv])
            trips.push_back(Eigen::Triplet<scalarD,sizeType>(id,vert->_verts[cv]->_index,Ms[v]*vert->_weights[cv]));
      } else trips.push_back(Eigen::Triplet<scalarD,sizeType>(id,vertP[v]->_index,Ms[v]));
    }
}
boost::shared_ptr<DHMVertex> getSubdVert
(const Vec3& tc,scalar coef,const DHMCell& cell,
 const DHMTraits::VERTPTRS& baseVerts,
 const boost::unordered_set<Vec4i,Hash>& fMap,
 boost::unordered_map<DHMTraits::CELLID,boost::shared_ptr<DHMVertex>,Hash>& fvertMap,
 boost::unordered_map<Vec4i,boost::shared_ptr<DHMVertex>,Hash>& ivertMap)
{
  Eigen::Matrix<scalar,8,1> Ms;
  DHMMapping::MStencil(Ms,tc*2.0f-Vec3::Ones());
  if(tc[0] == 0.0f || tc[0] == 1.0f ||
     tc[1] == 0.0f || tc[1] == 1.0f ||
     tc[2] == 0.0f || tc[2] == 1.0f) {
    DHMTraits::CELLID cid=DHMTraits::CELLID::Constant(numeric_limits<sizeType>::max());
    char d=0;
    for(char v=0; v<8; v++)
      if(Ms[v] > 0.0f) {
        cid[d]=cell._verts[v]->_index;
        cid[d+4]=(sizeType)(Ms[v]*coef);
        d++;
      }
    //in this case subdivided vertex coincide with coarse vertex
    if(d == 1)
      return baseVerts[cid[0]];
    //otherwise the subdivided vertex is an internal face vertex
    sort4Map(cid[0],cid[1],cid[2],cid[3],
             cid[4],cid[5],cid[6],cid[7]);
    if(fvertMap.find(cid) == fvertMap.end()) {
      //we need to create this vertex
      boost::shared_ptr<DHMVertex> vert;
      if(fMap.find(cid.block<4,1>(0,0)) == fMap.end()) {
        vert.reset(new DHMVertex);
      } else {
        vert.reset(new DHMConstrainedVertex);
        for(char v=0; v<d; v++)
          if(cid[v] != numeric_limits<sizeType>::max()) {
            boost::dynamic_pointer_cast<DHMConstrainedVertex>(vert)->_verts[v]=baseVerts[cid[v]];
            boost::dynamic_pointer_cast<DHMConstrainedVertex>(vert)->_weights[v]=(scalar)cid[v+4]/(scalar)coef;
          }
      }
      vert->_pos.setZero();
      for(char v=0; v<d; v++)
        vert->_pos+=baseVerts[cid[v]]->_pos*(scalar)cid[v+4]/(scalar)coef;
      fvertMap[cid]=vert;
      return vert;
    } else //this vertex has been created
      return fvertMap[cid];
  } else {
    Vec4i cid;
    cid[3]=cell._index;
    cid.block<3,1>(0,0)=(tc*coef).cast<sizeType>();
    //we need to create this vertex
    if(ivertMap.find(cid) == ivertMap.end()) {
      boost::shared_ptr<DHMVertex> vert(new DHMVertex);
      vert->_pos.setZero();
      for(char v=0; v<8; v++)
        vert->_pos+=cell._verts[v]->_pos*Ms[v];
      ivertMap[cid]=vert;
      return vert;
    } else //this vertex has been created
      return ivertMap[cid];
  }
}
void DHMAdaptiveMesh::createAdaptiveMesh(const vector<sizeType>& cells,sizeType lv)
{
  clearAdaptiveMesh();
  if(cells.empty())
    return;

  //initialize
  boost::unordered_set<Vec4i,Hash> fMap;
  boost::unordered_map<CELLID,boost::shared_ptr<DHMVertex>,Hash> fvertMap;
  boost::unordered_map<Vec4i,boost::shared_ptr<DHMVertex>,Hash> ivertMap;

  //build subdivision tag
  sizeType nrSC=(sizeType)cells.size();
  boost::shared_ptr<DHMMesh> baseMesh=_hier._mesh.front();
  vector<sizeType> subd(baseMesh->nrCell(),numeric_limits<sizeType>::max());
  for(sizeType i=0; i<nrSC; i++)
    subd[cells[i]]=0;
  for(sizeType i=0; i<baseMesh->nrFace(); i++) {
    const DHMFace& face=baseMesh->getFace(i);
    if(face._cells[1] && subd[face._cells[0]->_index] != subd[face._cells[1]->_index]) {
      Vec4i fid=face.getId();
      fMap.insert(fid);
      for(char vi=0; vi<4; vi++)
        for(char vj=vi+1; vj<4; vj++) {
          Vec4i eid=Vec4i::Constant(numeric_limits<sizeType>::max());
          eid[0]=fid[vi];
          eid[1]=fid[vj];
          fMap.insert(eid);
        }
    }
  }
  CELLPTRS otherCells;
  for(sizeType i=0; i<baseMesh->nrCell(); i++)
    if(subd[i] == numeric_limits<sizeType>::max()) {
      subd[i]=(sizeType)otherCells.size();
      otherCells.push_back(baseMesh->_cells[i]);
    }

  //build adaptive mesh
  sizeType nrSubd=1<<lv;
  sizeType offset=(nrSubd+1)*(nrSubd+1)*(nrSubd+1);
  Vec3i stride((nrSubd+1)*(nrSubd+1),(nrSubd+1),1);
  scalar coef=(scalar)(nrSubd*nrSubd*nrSubd);
  VERTPTRS vertMap(nrSC*offset);
  VERTPTRS verts=baseMesh->_verts;
  VERTPTRS cverts;
  for(sizeType step=nrSubd>>1,l=0; step>=1; step>>=1,l++) {
    //build vertex and cells
    _hier._mesh.insert(_hier._mesh.begin(),boost::shared_ptr<DHMMesh>(new DHMMesh));
    DHMMesh& lv=*(_hier._mesh.front());
    lv._cells=otherCells;
    sizeType lastVerts=(sizeType)verts.size();
    for(sizeType i=0,off=0; i<nrSC; i++,off+=offset) {
      const DHMCell& cell=baseMesh->getCell(cells[i]);
      subd[cell._index]=-1-(sizeType)lv._cells.size();
      for(sizeType x=0; x<nrSubd; x+=step)
        for(sizeType y=0; y<nrSubd; y+=step)
          for(sizeType z=0; z<nrSubd; z+=step) {
            //we reuse DHMCell::_index to store BaseMesh cell index
            boost::shared_ptr<DHMCell> scell(new DHMCell(cell._index));
            scell->_layer=cell._layer;
            lv._cells.push_back(scell);
            for(char v=0; v<8; v++) {
              Vec3i tci(x+(v&1?step:0),y+(v&2?step:0),z+(v&4?step:0));
              boost::shared_ptr<DHMVertex>& svert=vertMap[off+tci.dot(stride)];
              if(!svert) {
                svert=getSubdVert(tci.cast<scalar>()/(scalar)nrSubd,coef,cell,verts,fMap,fvertMap,ivertMap);
                if(svert->_index == numeric_limits<sizeType>::max())
                  if(boost::dynamic_pointer_cast<DHMConstrainedVertex>(svert)) {
                    svert->_index=(sizeType)cverts.size();
                    cverts.push_back(svert);
                  } else {
                    svert->_index=(sizeType)verts.size();
                    verts.push_back(svert);
                  }
              }
              scell->_verts[v]=svert;
            }
          }
    }
    //put cverts in the end
    lv._verts.insert(lv._verts.end(),verts.begin(),verts.end());
    lv._verts.insert(lv._verts.end(),cverts.begin(),cverts.end());
    lv._adaptive=subd;

    //build prolongation
    TRIPS tripsP;
    sizeType ccellSz=step*2;
    boost::shared_ptr<DHMVertex> vertP[8];
    vector<bool> visited(lv._verts.size(),false);
    for(sizeType i=0; i<baseMesh->nrVert(); i++) {
      visited[i]=true;
      tripsP.push_back(Eigen::Triplet<scalarD,sizeType>(i,i,1.0f));
    }
    for(sizeType i=0,off=0; i<nrSC; i++,off+=offset) {
      const DHMCell& cell=baseMesh->getCell(cells[i]);
      for(sizeType x=0; x<=nrSubd; x+=step)
        for(sizeType y=0; y<=nrSubd; y+=step)
          for(sizeType z=0; z<=nrSubd; z+=step) {
            ASSERT(vertMap[off+Vec3i(x,y,z).dot(stride)])
            boost::shared_ptr<DHMVertex> vert=vertMap[off+Vec3i(x,y,z).dot(stride)];
            if(!visited[vert->_index] && !boost::dynamic_pointer_cast<DHMConstrainedVertex>(vert)) {
              visited[vert->_index]=true;
              Vec3i base=(Vec3i(x,y,z)/ccellSz*ccellSz).cwiseMin(Vec3i::Constant(nrSubd-ccellSz));
              for(char v=0; v<8; v++) {
                Vec3i tcv=base+Vec3i(v&1?ccellSz:0,v&2?ccellSz:0,v&4?ccellSz:0);
                vertP[v]=vertMap[off+tcv.dot(stride)];
              }
              addPRow(vert->_index,(Vec3i(x,y,z)-base)/step,vertP,tripsP);
            }
          }
    }
    _hier._prolong.insert(_hier._prolong.begin(),SMAT());
    _hier._prolong.front().resize(verts.size(),lastVerts);
    _hier._prolong.front().buildFromTripletsDepulicate(tripsP,0.0f);
  }
}
void DHMAdaptiveMesh::parityCheck() const
{
  for(sizeType i=0; i<(sizeType)_hier._prolong.size(); i++)
    if(_hier._prolong[i].rows() < _hier._mesh[i]->nrVert()) {
      //check that DHMConstrainedVertex are inserted to the back
      const DHMMesh& lv=*(_hier._mesh[i]);
      const SMAT& lvP=_hier._prolong[i];
      for(sizeType v=0; v<lvP.rows(); v++)
        ASSERT(!boost::dynamic_pointer_cast<DHMConstrainedVertex>(lv._verts[v]))
        for(sizeType v=lvP.rows(); v<(sizeType)lv.nrVert(); v++) {
          ASSERT(boost::dynamic_pointer_cast<DHMConstrainedVertex>(lv._verts[v]))
        }

      //check that prolongation is exact
      const VERTPTRS& last=_hier._mesh[i+1]->_verts;
      for(sizeType v=0; v<lvP.rows(); v++) {
        Vec3 pos=Vec3::Zero();
        for(ConstSMIterator<scalarD> beg=lvP.begin(v),end=lvP.end(v); beg!=end; ++beg)
          pos+=last[beg.col()]->_pos*(scalar)(*beg);
        ASSERT((pos-lv._verts[v]->_pos).norm() < 1E-4f)
      }
    }
}

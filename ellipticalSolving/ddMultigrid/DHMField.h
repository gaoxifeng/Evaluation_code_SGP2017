#ifndef DHM_FIELD_H
#define DHM_FIELD_H

#include "CommonFile/ImplicitFuncInterface.h"
#include "CommonFile/Zero.h"
#include "DHMAdaptiveMesh.h"
#include "DHMIter.h"
#include "DHMUtil.h"

PRJ_BEGIN

//DHM data field
template <typename T,bool CELL,char c>
struct DHMField : public DHMTraits, public DHMEnergyTraits<T> {
  using typename DHMEnergyTraits<T>::COL;
  typedef typename Eigen::Matrix<T,c,1> COMP;
  typedef typename Eigen::Block<COL,c,1> COMPREF;
  typedef typename Eigen::Matrix<T,3,1> PT3;
  DHMField(const DHMMesh& mesh)
    :_mesh(mesh) {
    clear(Zero<T>::value());
  }
  //IO
  bool read(std::istream& is) {
    return readBinaryData(_data,is).good();
  }
  bool write(std::ostream& os) const {
    return writeBinaryData(_data,os).good();
  }
  void writeSVTK(const std::string& path) const {
    _mesh.assembleVIndex();
    VTKWriter<T> os("DHMMesh",path,true);
    os.appendPoints(VertPosIter(_mesh,0),VertPosIter(_mesh,_mesh.nrVert()));
    os.appendCells(CellVidIter(_mesh,0),CellVidIter(_mesh,_mesh.nrCell()),VTKWriter<T>::VOXEL);
    for(char i=0; i<c; i++) {
      ostringstream oss;
      oss << "comp" << (sizeType)i;
      if(CELL) {
        vector<T> vals;
        for(sizeType ci=0; ci<_mesh.nrCell(); ci++)
          vals.push_back(getComp(ci)[i]);
        os.appendCustomData(oss.str(),vals.begin(),vals.end());
      } else {
        vector<T> vals;
        for(sizeType vi=0; vi<_mesh.nrVert(); vi++)
          vals.push_back(getVScalar<T>(i,_mesh.getVert(vi)));
        os.appendCustomPointData(oss.str(),vals.begin(),vals.end());
      }
    }
  }
  void writeVVTK(const std::string& path) const {
    _mesh.assembleVIndex();
    VTKWriter<T> os("DHMMesh",path,true);
    os.appendPoints(VertPosIter(_mesh,0),VertPosIter(_mesh,_mesh.nrVert()));
    os.appendCells(CellVidIter(_mesh,0),CellVidIter(_mesh,_mesh.nrCell()),VTKWriter<T>::VOXEL);
    if(CELL) {
      for(char i=0; i<c/8; i++) {
        ostringstream oss;
        oss << "comp" << (sizeType)i;

        vector<PT3,Eigen::aligned_allocator<PT3> > vss;
        for(sizeType ci=0; ci<_mesh.nrCell(); ci++)
          vss.push_back(getVel(i,ci,PT3::Zero()));
        os.appendCustomVectorData(oss.str(),vss.begin(),vss.end());
      }
    } else {
      for(char i=0; i<c/3; i++) {
        ostringstream oss;
        oss << "comp" << (sizeType)i;

        vector<PT3,Eigen::aligned_allocator<PT3> > vss;
        for(sizeType vi=0; vi<_mesh.nrVert(); vi++)
          vss.push_back(PT3(getVScalar<T>(i*3+0,_mesh.getVert(vi)),getVScalar<T>(i*3+1,_mesh.getVert(vi)),getVScalar<T>(i*3+2,_mesh.getVert(vi))));
        os.appendCustomPointVectorData(oss.str(),vss.begin(),vss.end());
      }
    }
  }
  void writeDVVTK(const std::string& path) const {
    _mesh.assembleVIndex();
    VTKWriter<T> os("DHMMesh",path,true);
    POSES poses;
    for(sizeType i=0; i<_mesh.nrCell(); i++) {
      const DHMCell& cell=_mesh.getCell(i);
      for(char v=0; v<8; v++)
        poses.push_back(cell._verts[v]->_pos);
    }
    os.appendPoints(poses.begin(),poses.end());
    os.appendCells(typename VTKWriter<T>::template IteratorIndex<CELLID>(0,8,0),
                   typename VTKWriter<T>::template IteratorIndex<CELLID>(_mesh.nrCell(),8,0),
                   VTKWriter<T>::VOXEL);
    for(char i=0; i<c/8; i++) {
      ostringstream oss;
      oss << "comp" << (sizeType)i;
      vector<T> data;
      for(sizeType ci=0; ci<_mesh.nrCell(); ci++)
        for(char v=0; v<8; v++)
          data.push_back(getComp(ci)[i*8+v]);
      os.appendCustomPointData(oss.str(),data.begin(),data.end());
    }
  }
  template <typename VEC>
  void writeSVVTK(const std::string& path,const VEC* sv=NULL,sizeType C=0,Vec3c nrP=Vec3c::Constant(2)) const {
    vector<sizeType> ids;
    if(sv) {
      std::copy(sv->begin(),sv->end(),std::back_inserter(ids));
    } else {
      for(sizeType i=0; i<_mesh.nrCell(); i++)
        ids.push_back(i);
    }

    ASSERT(CELL)
    POSES poses;
    vector<T> css;
    for(sizeType i=0; i<(sizeType)ids.size(); i++) {
      const DHMCell& cell=_mesh.getCell(ids[i]);
      for(char x=0; x<nrP[0]; x++)
        for(char y=0; y<nrP[1]; y++)
          for(char z=0; z<nrP[2]; z++) {
            Vec3 crd=(PT3(x,y,z)+PT3::Constant(0.5f))*2.0f;
            crd.array()/=nrP.array().cast<scalar>();
            crd-=PT3::Ones();
            PT3 pos=DHMMapping::M(cell,crd);
            poses.push_back(pos);
            poses.push_back(pos+getVel(C,ids[i],crd));

            css.push_back(0.0f);
            css.push_back(1.0f);
          }
    }

    VTKWriter<T> os("DHMMesh",path,true);
    os.appendPoints(poses.begin(),poses.end());
    os.appendCells(typename VTKWriter<T>::template IteratorIndex<Vec3i>(0,2,0),
                   typename VTKWriter<T>::template IteratorIndex<Vec3i>((sizeType)poses.size()/2,2,0),
                   VTKWriter<T>::LINE);
    os.appendCustomPointData("Color",css.begin(),css.end());
  }
  //utility
  PT3 getVel(sizeType comp,sizeType cid,const PT3& crd) const {
    if(CELL) {
      Eigen::Matrix<T,3,3> J;
      Eigen::Matrix<T,3,8> Js;
      DHMMapping::JStencil(Js,crd);
      DHMMapping::J(Js,_mesh.getCell(cid),J);
      return J.inverse().transpose()*(Js*_data.template block<8,1>(cid*c+comp*8,0)).eval();
    } else {
      return PT3(getScalar<T>(0,cid,crd),getScalar<T>(1,cid,crd),getScalar<T>(2,cid,crd));
    }
  }
  template <typename T2>
  T2 getScalar(sizeType comp,sizeType cid,const PT3& crd,T2* MIN=NULL,T2* MAX=NULL) const {
    ASSERT(!CELL)
    const DHMCell& cell=_mesh.getCell(cid);
    Eigen::Matrix<T,8,1> Ms,cellVal;
    DHMMapping::MStencil(Ms,crd);

    T2 ret=0.0f,val;
    for(char v=0; v<8; v++) {
      val=getVScalar<T2>(comp,*(cell._verts[v]));
      ret+=Ms[v]*val;
      if(MIN)
        *MIN=std::min(*MIN,val);
      if(MAX)
        *MAX=std::max(*MAX,val);
    }
    return ret;
  }
  template <typename T2>
  T2 getVScalar(sizeType comp,const DHMVertex& vert) const {
    const DHMConstrainedVertex* cv=
      dynamic_cast<const DHMConstrainedVertex*>(&vert);
    T2 val=0.0f;
    if(cv) {
      for(char cvi=0; cvi<4; cvi++)
        if(cv->_verts[cvi])
          val+=(T2)(cv->_weights[cvi]*_data[cv->_verts[cvi]->_index*c+comp]);
    } else val=(T2)_data[vert._index*c+comp];
    return val;
  }
  COMPREF getCompRef(sizeType i) {
    return _data.template block<c,1>(i*c,0);
  }
  COMP getComp(sizeType i) const {
    return _data.template block<c,1>(i*c,0);
  }
  sizeType nrComp() const {
    if(CELL)
      return _mesh.nrCell();
    else return _mesh.nrVertNC();
  }
  void clear(const T& val) {
    if(CELL)
      _data.setConstant(_mesh.nrCell()*c,val);
    else _data.setConstant(_mesh.nrVertNC()*c,val);
  }
  void swap(DHMField& other) {
    ASSERT(&(other._mesh) == &_mesh);
    ASSERT(other._data.size() == _data.size());
    other._data.swap(_data);
  }
  //data
  const DHMMesh& _mesh;
  COL _data;
};
typedef DHMField<scalar,false,3> DHMVVField;
typedef DHMField<scalar,true,8> DHMCVField;
typedef DHMField<scalar,false,1> DHMVSField;
typedef DHMField<scalar,true,1> DHMCSField;
//DHM field converter
struct DHMFieldConverter : public DHMTraits, public DHMEnergyTraits<scalar> {
  //converter
  static void convert(DHMVVField& field,const DHMCVField& func);
  static void convert(DHMVVField& field,const VelFunc<scalar>& func);
  static void convert(DHMVSField& field,const ImplicitFunc<scalar>& func);
  static void convert(DHMCSField& field,const ImplicitFunc<scalar>& func);
  static void convertAdaptive(DHMCVField& field,const COL& samples,bool useAdaptiveBasis=true);
  static void convertAdaptive(DHMCVField& field,const VelFunc<scalar>& func,bool useAdaptiveBasis=true);
  static void precomputeConvertAdaptive(DHMCVField& field);
  static void convertAdaptiveDirect(DHMCVField& field,const COL& samples);
  static void convertAdaptiveDirect(DHMCVField& field,const VelFunc<scalar>& func);
};

PRJ_END

#endif

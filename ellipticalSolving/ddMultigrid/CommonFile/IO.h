#ifndef IO_H
#define IO_H

#include "IOBasic.h"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/interprocess/streams/vectorstream.hpp>

PRJ_BEGIN

using namespace std;

//io for pair
template <typename T1,typename T2>
FORCE_INLINE ostream& writeBinaryData(const std::pair<T1,T2>& m,ostream& os,IOData* dat=NULL)
{
  writeBinaryData(m.first,os,dat);
  writeBinaryData(m.second,os,dat);
  return os;
}
template <typename T1,typename T2>
FORCE_INLINE istream& readBinaryData(std::pair<T1,T2>& m,istream& is,IOData* dat=NULL)
{
  readBinaryData(m.first,is,dat);
  readBinaryData(m.second,is,dat);
  return is;
}

//io for string
istream& readBinaryData(string& str,istream& is,IOData* dat=NULL);
ostream& writeBinaryData(const string& str,ostream& os,IOData* dat=NULL);

//io for vector
template <typename T,typename ALLOC>FORCE_INLINE istream& readVector(vector<T,ALLOC>& v,istream& is,IOData* dat=NULL)
{
  sizeType size;
  readBinaryData(size,is,dat);
  v.resize(size);
  for(sizeType i=0; i<(sizeType)v.size(); i++)
    readBinaryData(v[i],is,dat);
  return is;
}
template <typename T,typename ALLOC>FORCE_INLINE ostream& writeVector(const vector<T,ALLOC>& v,ostream& os,IOData* dat=NULL)
{
  sizeType size=(sizeType)v.size();
  writeBinaryData(size,os,dat);
  for(sizeType i=0; i<(sizeType)v.size(); i++)
    writeBinaryData(v[i],os,dat);
  return os;
}

//io for float vector is double vector
istream& readVector(vector<scalarD,Eigen::aligned_allocator<scalarD> >& v,istream& is,IOData* dat=NULL);
ostream& writeVector(const vector<scalarD,Eigen::aligned_allocator<scalarD> >& v,ostream& os,IOData* dat=NULL);
istream& readVector(vector<scalarF,Eigen::aligned_allocator<scalarF> >& v,istream& is,IOData* dat=NULL);
ostream& writeVector(const vector<scalarF,Eigen::aligned_allocator<scalarF> >& v,ostream& os,IOData* dat=NULL);

//io for vector of Eigen::Matrix
#define IO_VEC_FIXED_DECL(NAME)	\
istream& readVector(vector<NAME,Eigen::aligned_allocator<NAME> >& v,istream& is,IOData* dat=NULL);	\
ostream& writeVector(const vector<NAME,Eigen::aligned_allocator<NAME> >& v,ostream& os,IOData* dat=NULL);
IO_VEC_FIXED_DECL(Vec2f)
IO_VEC_FIXED_DECL(Vec3f)
IO_VEC_FIXED_DECL(Vec4f)
IO_VEC_FIXED_DECL(Mat2f)
IO_VEC_FIXED_DECL(Mat3f)
IO_VEC_FIXED_DECL(Mat4f)
IO_VEC_FIXED_DECL(Quatf)
IO_VEC_FIXED_DECL(Transf)
IO_VEC_FIXED_DECL(Affinef)
IO_VEC_FIXED_DECL(Vec2d)
IO_VEC_FIXED_DECL(Vec3d)
IO_VEC_FIXED_DECL(Vec4d)
IO_VEC_FIXED_DECL(Mat2d)
IO_VEC_FIXED_DECL(Mat3d)
IO_VEC_FIXED_DECL(Mat4d)
IO_VEC_FIXED_DECL(Quatd)
IO_VEC_FIXED_DECL(Transd)
IO_VEC_FIXED_DECL(Affined)
IO_VEC_FIXED_DECL(Vec2i)
IO_VEC_FIXED_DECL(Vec3i)
IO_VEC_FIXED_DECL(Vec4i)
#undef IO_VEC_FIXED_DECL

//utility
FORCE_INLINE bool beginsWith(const string& name,const string& regex)
{
  return name.size() >= regex.size() && name.substr(0,regex.size()) == regex;
}
FORCE_INLINE bool endsWith(const string& name,const string& regex)
{
  return name.size() >= regex.size() && name.substr(name.size()-regex.size(),regex.size()) == regex;
}
FORCE_INLINE string replaceDot(const string& name)
{
  string ret=name;
  for(sizeType i=0,nr=(sizeType)ret.size(); i<nr; i++)
    if(ret[i]=='.')
      ret[i]='_';
  return ret;
}

//for file format definition
class HasMagic
{
public:
  HasMagic(const sizeType& MAGIC);
  virtual ~HasMagic() {}
  bool readMagic(istream& is);
  bool writeMagic(ostream& os) const;
  HasMagic& operator=(const HasMagic& other);
  const sizeType _MAGIC;
};

//for VTK IO
void vtkWrite(ostream& oss,scalarF val);
void vtkWrite(ostream& oss,scalarD val);
void vtkWrite(ostream& oss,int val);
void vtkWrite(ostream& oss,unsigned char val);
template<typename T>
struct VTKWriter {
  enum VTK_DATA_TYPE {
    UNSTRUCTURED_GRID,
    STRUCTURED_POINTS,
  };
  enum VTK_CELL_TYPE {
    POINT=1,
    LINE=3,
    TRIANGLE=5,
    TETRA=10,
    VOXEL=11,
    HEX=12,
    POLYLINE=4,
    QUAD=9,
    QUADRATIC_LINE=21,
  };
public:
  struct Data {
    Data():_nr(0) {}
    string _str;
    sizeType _nr;
  };
  template <typename ITER>
  struct ValueTraits {
    typedef typename ITER::value_type value_type;
  };
  template <typename POINTED_TYPE>
  struct ValueTraits<POINTED_TYPE*> {
    typedef POINTED_TYPE value_type;
  };
  template <typename ID>
  struct IteratorIndex {
    typedef ID value_type;
    IteratorIndex(const sizeType& id,const sizeType stride,const sizeType& off)
      :_id(id),_stride(stride),_off(off) {}
    void operator++() {
      _id++;
    }
    bool operator!=(const IteratorIndex& other) const {
      return _id != other._id;
    }
    virtual ID operator*() const {
      ID ret;
      for(sizeType i=0; i<ret.size(); i++)ret(i)=(_stride == 0) ? _id+_off*i : _id*_stride+i;
      return ret;
    }
    sizeType _id;
    sizeType _stride;
    sizeType _off;
  };
  template <typename ID>
  struct IteratorRepeat {
    typedef ID value_type;
    IteratorRepeat(const sizeType& id,const ID& val)
      :_id(id),_val(val) {}
    void operator++() {
      _id++;
    }
    bool operator!=(const IteratorRepeat& other) const {
      return _id != other._id;
    }
    virtual ID operator*() const {
      return _val;
    }
    sizeType _id;
    ID _val;
  };
  template <typename ITER>
  struct IteratorAdd {
    typedef typename ValueTraits<ITER>::value_type value_type;
    IteratorAdd(ITER beg0,ITER beg1):_beg0(beg0),_beg1(beg1) {}
    void operator++() {
      _beg0++;
      _beg1++;
    }
    bool operator!=(const IteratorAdd& other) const {
      return _beg0 != other._beg0;
    }
    virtual value_type operator*() const {
      return (*_beg0)+(*_beg1);
    }
    ITER _beg0,_beg1;
  };
  template <typename ITER,typename SCALAR>
  struct IteratorAddMult {
    typedef typename ValueTraits<ITER>::value_type value_type;
    IteratorAddMult(ITER beg0,ITER beg1,SCALAR mult):_beg0(beg0),_beg1(beg1),_mult(mult) {}
    void operator++() {
      _beg0++;
      _beg1++;
    }
    bool operator!=(const IteratorAddMult& other) const {
      return _beg0 != other._beg0;
    }
    virtual value_type operator*() const {
      return (*_beg0)+(*_beg1)*_mult;
    }
    ITER _beg0,_beg1;
    SCALAR _mult;
  };
public:
  VTKWriter(const string& name,const boost::filesystem::path& path,bool binary)
    :_os(path,binary ? ios_base::binary : ios_base::out),
     _points(binary ? ios_base::binary : ios_base::out),
     _cells(binary ? ios_base::binary : ios_base::out),
     _cellTypes(binary ? ios_base::binary : ios_base::out),
     _cellDatas(binary ? ios_base::binary : ios_base::out),
     _nrPoint(0),_nrCell(0),_nrIndex(0),_nrData(0),_vdt(UNSTRUCTURED_GRID),
     _binary(binary) {
    _os << "# vtk DataFile Version 1.0" << endl;
    _os << name << endl;
    _os << (binary ? "BINARY" : "ASCII") << endl;
    _os << "DATASET " << "UNSTRUCTURED_GRID" << endl;
  }
  VTKWriter(const string& name,const boost::filesystem::path& path,bool binary,const BBox<T>& bb,const Vec3i& nrCell,bool center)
    :_os(path,binary ? ios_base::binary : ios_base::out),
     _points(binary ? ios_base::binary : ios_base::out),
     _cells(binary ? ios_base::binary : ios_base::out),
     _cellTypes(binary ? ios_base::binary : ios_base::out),
     _cellDatas(binary ? ios_base::binary : ios_base::out),
     _nrPoint(0),_nrCell(0),_nrIndex(0),_nrData(0),_vdt(STRUCTURED_POINTS),
     _binary(binary) {
    typename BBox<T>::PT ext=bb.getExtent();
    typename BBox<T>::PT spacing(ext.x()/nrCell.x(),ext.y()/nrCell.y(),ext.z()/nrCell.z());
    _os << "# vtk DataFile Version 1.0" << endl;
    _os << name << endl;
    _os << (binary ? "BINARY" : "ASCII") << endl;
    _os << "DATASET " << "STRUCTURED_POINTS" << endl;
    if(center) {
      typename BBox<T>::PT origin=bb._minC+spacing*0.5f;
      _os << "DIMENSIONS " << nrCell.x() << " " << nrCell.y() << " " << nrCell.z() << endl;
      _os << "ORIGIN " << origin.x() << " " << origin.y() << " " << origin.z() << endl;
      _os << "SPACING " << spacing.x() << " " << spacing.y() << " " << spacing.z() << endl;
    } else {
      typename BBox<T>::PT origin=bb._minC;
      _os << "DIMENSIONS " << (nrCell.x()+1) << " " << (nrCell.y()+1) << " " << (nrCell.z()+1) << endl;
      _os << "ORIGIN " << origin.x() << " " << origin.y() << " " << origin.z() << endl;
      _os << "SPACING " << spacing.x() << " " << spacing.y() << " " << spacing.z() << endl;
    }
  }
  virtual ~VTKWriter() {
    bool first;
    switch(_vdt) {
    case UNSTRUCTURED_GRID:
      _os << "POINTS " << _nrPoint << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
      _os << _points.str();
      _os << "CELLS " << _nrCell << " " << _nrIndex << endl;
      _os << _cells.str();
      _os << "CELL_TYPES " << _nrCell << endl;
      _os << _cellTypes.str();
      first=false;
      for(typename boost::unordered_map<string,Data>::const_iterator beg=_customData.begin(),end=_customData.end(); beg!=end; beg++) {
        if(!first)
          _os << "CELL_DATA " << beg->second._nr << endl;
        first=true;
        _os << beg->second._str << endl;
      }
      first=false;
      for(typename boost::unordered_map<string,Data>::const_iterator beg=_customPointData.begin(),end=_customPointData.end(); beg!=end; beg++) {
        if(!first)
          _os << "POINT_DATA " << beg->second._nr << endl;
        first=true;
        _os << beg->second._str << endl;
      }
      break;
    case STRUCTURED_POINTS:
      ASSERT(_nrData == 1)
      _os << "POINT_DATA " << _nrPoint << std::endl;
      _os << "SCALARS data " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
      _os << "LOOKUP_TABLE default" << endl;
      _os << _cellDatas.str();
      break;
    default:
      ASSERT_MSG(false,"Unsupported!")
    }
  }
  template <typename ITER> VTKWriter& appendPoints(ITER beg,ITER end) {
    typedef typename ValueTraits<ITER>::value_type value_type;
    sizeType nr=0;
    if(_binary) {
      for(; beg != end; ++beg) {
        const value_type& val=*beg;
        vtkWrite(_points,val(0));
        vtkWrite(_points,val(1));
        vtkWrite(_points,val(2));
        nr++;
      }
    } else {
      for(; beg != end; ++beg) {
        const value_type& val=*beg;
        _points << (T)val(0) << " " << (T)val(1) << " " << (T)val(2) << endl;
        nr++;
      }
    }
    _nrPoint+=nr;
    return *this;
  }
  template <typename ITER> VTKWriter& appendVoxels(ITER beg,ITER end,bool hex) {
    typedef typename Eigen::Matrix<sizeType,8,1> IDS;
    typedef typename ValueTraits<ITER>::value_type value_type;
    std::vector<value_type,Eigen::aligned_allocator<value_type> > points;
    std::vector<IDS,Eigen::aligned_allocator<IDS> > cells;
    for(; beg!=end;) {
      IDS ids;
      value_type minC=*beg++;
      value_type maxC=*beg++;
      value_type ext=maxC-minC;

      if(hex) ids << 0,1,3,2,4,5,7,6;
      else ids << 0,1,2,3,4,5,6,7;
      ids.array()+=points.size();
      cells.push_back(ids);

      points.push_back(minC+value_type(0.0f   ,0.0f   ,0.0f   ));
      points.push_back(minC+value_type(ext.x(),   0.0f,0.0f   ));
      points.push_back(minC+value_type(0.0f   ,ext.y(),0.0f   ));
      points.push_back(minC+value_type(ext.x(),ext.y(),0.0f   ));

      points.push_back(minC+value_type(0.0f   ,0.0f   ,ext.z()));
      points.push_back(minC+value_type(ext.x(),   0.0f,ext.z()));
      points.push_back(minC+value_type(0.0f   ,ext.y(),ext.z()));
      points.push_back(minC+value_type(ext.x(),ext.y(),ext.z()));
    }
    setRelativeIndex();
    appendPoints(points.begin(),points.end());
    appendCells(cells.begin(),cells.end(),hex ? HEX : VOXEL,true);
    return *this;
  }
  template <typename ITER> VTKWriter& appendCells(ITER beg,ITER end,VTK_CELL_TYPE ct,bool relativeIndex=false) {
    if(relativeIndex)
      ASSERT(_relativeCellIndex >= -1)
      int base=relativeIndex ? (int)_relativeCellIndex : 0;

    typedef typename ValueTraits<ITER>::value_type value_type;
    sizeType nr=0;
    sizeType nrIndex=0;
    if(_binary) {
      for(; beg != end; ++beg) {
        const value_type& val=*beg;
        switch(ct) {
        case POINT:
          nrIndex+=2;
          vtkWrite(_cells,1);
          vtkWrite(_cells,base+(int)val(0));
          break;
        case LINE:
          nrIndex+=3;
          vtkWrite(_cells,2);
          vtkWrite(_cells,base+(int)val(0));
          vtkWrite(_cells,base+(int)val(1));
          break;
        case TRIANGLE:
        case QUADRATIC_LINE:
          nrIndex+=4;
          vtkWrite(_cells,3);
          vtkWrite(_cells,base+(int)val(0));
          vtkWrite(_cells,base+(int)val(1));
          vtkWrite(_cells,base+(int)val(2));
          break;
        case TETRA:
        case QUAD:
          nrIndex+=5;
          vtkWrite(_cells,4);
          vtkWrite(_cells,base+(int)val(0));
          vtkWrite(_cells,base+(int)val(1));
          vtkWrite(_cells,base+(int)val(2));
          vtkWrite(_cells,base+(int)val(3));
          break;
        case VOXEL:
        case HEX:
          nrIndex+=9;
          vtkWrite(_cells,8);
          vtkWrite(_cells,base+(int)val(0));
          vtkWrite(_cells,base+(int)val(1));
          vtkWrite(_cells,base+(int)val(2));
          vtkWrite(_cells,base+(int)val(3));
          vtkWrite(_cells,base+(int)val(4));
          vtkWrite(_cells,base+(int)val(5));
          vtkWrite(_cells,base+(int)val(6));
          vtkWrite(_cells,base+(int)val(7));
          break;
        case POLYLINE:
          nrIndex+=val.rows()+1;
          vtkWrite(_cells,(int)val.rows());
          for(sizeType i=0; i<(int)val.rows(); i++)
            vtkWrite(_cells,base+(int)val[i]);
          break;
        }
        vtkWrite(_cellTypes,(int)ct);
        nr++;
      }
    } else {
      for(; beg != end; ++beg) {
        const value_type& val=*beg;
        switch(ct) {
        case POINT:
          nrIndex+=2;
          _cells << "1 " << (base+(int)val(0)) << endl;
          break;
        case LINE:
          nrIndex+=3;
          _cells << "2 " << (base+(int)val(0)) << " " << (base+(int)val(1)) << endl;
          break;
        case TRIANGLE:
        case QUADRATIC_LINE:
          nrIndex+=4;
          _cells << "3 " << (base+(int)val(0)) << " " << (base+(int)val(1)) << " " << (base+(int)val(2)) << endl;
          break;
        case TETRA:
        case QUAD:
          nrIndex+=5;
          _cells << "4 " << (base+(int)val(0)) << " " << (base+(int)val(1)) << " " << (base+(int)val(2)) << " " << (base+(int)val(3)) << endl;
          break;
        case VOXEL:
        case HEX:
          nrIndex+=9;
          _cells << "8 " << (base+(int)val(0)) << " " << (base+(int)val(1)) << " " << (base+(int)val(2)) << " " << (base+(int)val(3)) << " "
                 << (base+(int)val(4)) << " " << (base+(int)val(5)) << " " << (base+(int)val(6)) << " " << (base+(int)val(7)) << endl;
          break;
        case POLYLINE:
          nrIndex+=val.rows()+1;
          _cells << val.rows() << " ";
          for(sizeType i=0; i<(int)val.rows(); i++)
            _cells << (base+(int)val[i]) << " ";
          _cells << endl;
          break;
        }
        _cellTypes << ct << endl;
        nr++;
      }
    }
    _nrCell+=nr;
    _nrIndex+=nrIndex;
    return *this;
  }
  template <typename ITER> VTKWriter& appendDatas(const string name,ITER beg,ITER end) {
    if(_binary)
      for(; beg != end; ++beg,++_nrPoint)
        vtkWrite(_cellDatas,*beg);
    else
      for(; beg != end; ++beg,++_nrPoint)
        _cellDatas << (T)*beg << endl;
    _nrData++;
    return *this;
  }
  template <typename ITER> VTKWriter& appendCustomData(const string name,ITER beg,ITER end) {
    ostringstream os;
    if(_customData.find(name) == _customData.end()) {
      os << "SCALARS " << name << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
      os << "LOOKUP_TABLE default" << endl;
    }

    Data& dat=_customData[name];
    if(_binary)
      for(; beg != end; ++beg,++dat._nr)
        vtkWrite(os,(T)*beg);
    else
      for(; beg != end; ++beg,++dat._nr)
        os << (T)*beg << endl;
    dat._str+=os.str();
    ASSERT(dat._nr == _nrCell)
    return *this;
  }
  template <typename ITER> VTKWriter& appendCustomPointData(const string name,ITER beg,ITER end) {
    ostringstream os;
    if(_customPointData.find(name) == _customPointData.end()) {
      //Data& dat=_customPointData[name];
      os << "SCALARS " << name << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
      os << "LOOKUP_TABLE default" << endl;
    }

    Data& dat=_customPointData[name];
    if(_binary)
      for(; beg != end; ++beg,++dat._nr)
        vtkWrite(os,(T)*beg);
    else
      for(; beg != end; ++beg,++dat._nr)
        os << (T)*beg << endl;
    dat._str+=os.str();
    ASSERT(dat._nr == _nrPoint)
    return *this;
  }
  template <typename ITER> VTKWriter& appendCustomVectorData(const string name,ITER beg,ITER end) {
    ostringstream os;
    if(_customData.find(name) == _customData.end()) {
      os << "VECTORS " << name << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
    }

    Data& dat=_customData[name];
    if(_binary)
      for(; beg != end; ++beg,++dat._nr) {
        vtkWrite(os,(T)(*beg)[0]);
        vtkWrite(os,(T)(*beg)[1]);
        vtkWrite(os,(T)(*beg)[2]);
      }
    else
      for(; beg != end; ++beg,++dat._nr)
        os << (T)(*beg)[0] << " " << (T)(*beg)[1] << " " << (T)(*beg)[2] << endl;
    dat._str+=os.str();
    ASSERT(dat._nr == _nrCell)
    return *this;
  }
  template <typename ITER> VTKWriter& appendCustomPointVectorData(const string name,ITER beg,ITER end) {
    ostringstream os;
    if(_customPointData.find(name) == _customPointData.end()) {
      //Data& dat=_customPointData[name];
      os << "VECTORS " << name << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
    }

    Data& dat=_customPointData[name];
    if(_binary)
      for(; beg != end; ++beg,++dat._nr) {
        vtkWrite(os,(T)(*beg)[0]);
        vtkWrite(os,(T)(*beg)[1]);
        vtkWrite(os,(T)(*beg)[2]);
      }
    else
      for(; beg != end; ++beg,++dat._nr)
        os << (T)(*beg)[0] << " " << (T)(*beg)[1] << " " << (T)(*beg)[2] << endl;
    dat._str+=os.str();
    ASSERT(dat._nr == _nrPoint)
    return *this;
  }
  template <typename ITER> VTKWriter& appendPointsByAdd(ITER beg0,ITER beg1,ITER end0) {
    appendPoints(IteratorAdd<ITER>(beg0,beg1),IteratorAdd<ITER>(end0,end0));
    return *this;
  }
  template <typename ITER> VTKWriter& appendCustomPointColorData(const string name,ITER begC,ITER endC) {
    ostringstream os;
    if(_customPointData.find(name) == _customPointData.end()) {
      //Data& dat=_customPointData[name];
      os << "COLOR_SCALARS " << name << " " << (*begC).size() << endl;
    }

    int sz=(*begC).size();
    Data& dat=_customPointData[name];
    if(_binary) {
      for(; begC != endC; ++begC,++dat._nr)
        for(int d=0; d<sz; d++)
          vtkWrite(os,(unsigned char)((*begC)[d]*255));
    } else {
      for(; begC != endC; ++begC,++dat._nr) {
        for(int d=0; d<sz; d++)
          os << (T)(*begC)[d] << " ";
        os << endl;
      }
    }
    dat._str+=os.str();
    return *this;
  }
  void setRelativeIndex(sizeType rel=-1) {
    if(rel == -1)
      _relativeCellIndex=_nrPoint;
    else _relativeCellIndex=rel;
  }
private:
  boost::filesystem::ofstream _os;
  ostringstream _points,_cells,_cellTypes,_cellDatas;
  boost::unordered_map<string,Data> _customData;
  boost::unordered_map<string,Data> _customPointData;
  sizeType _nrPoint,_nrCell,_nrIndex,_nrData;
  VTK_DATA_TYPE _vdt;
  bool _binary;
  sizeType _relativeCellIndex;
};

PRJ_END

#endif

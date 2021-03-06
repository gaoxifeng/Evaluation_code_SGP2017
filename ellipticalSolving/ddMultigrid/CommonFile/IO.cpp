#include "IO.h"

PRJ_BEGIN

istream& readBinaryData(string& str,istream& is,IOData* dat)
{
  sizeType len;
  readBinaryData(len,is,dat);
  str.assign(len,' ');
  is.read(&(str[0]),len);
  return is;
}
ostream& writeBinaryData(const string& str,ostream& os,IOData* dat)
{
  const sizeType len=(sizeType)str.length();
  writeBinaryData(len,os,dat);
  return os.write(str.c_str(),len);
}

istream& readVector(vector<scalarD,Eigen::aligned_allocator<scalarD> >& v,istream& is,IOData* dat)
{
  sizeType size;
  readBinaryData(size,is,dat);
  v.resize(size);
  if(size > 0)
    is.read((char*)&(v[0]),sizeof(scalarD)*size);
  return is;
}
ostream& writeVector(const vector<scalarD,Eigen::aligned_allocator<scalarD> >& v,ostream& os,IOData* dat)
{
  sizeType size=(sizeType)v.size();
  writeBinaryData(size,os,dat);
  if(size > 0)
    os.write((char*)&(v[0]),sizeof(scalarD)*size);
  return os;
}
istream& readVector(vector<scalarF,Eigen::aligned_allocator<scalarF> >& v,istream& is,IOData* dat)
{
  std::vector<scalarD,Eigen::aligned_allocator<scalarD> > tmp;
  readVector(tmp,is,dat);

  v.resize(tmp.size());
  for(sizeType i=0; i<(sizeType)v.size(); i++)
    v[i]=convert<scalarF>()(tmp[i]);
  return is;
}
ostream& writeVector(const vector<scalarF,Eigen::aligned_allocator<scalarF> >& v,ostream& os,IOData* dat)
{
  std::vector<scalarD,Eigen::aligned_allocator<scalarD> > tmp(v.size());
  for(sizeType i=0; i<(sizeType)v.size(); i++)
    tmp[i]=v[i];

  writeVector(tmp,os,dat);
  return os;
}

#define IO_VEC_FIXED(NAME)																									\
istream& readVector(vector<NAME,Eigen::aligned_allocator<NAME> >& v,istream& is,IOData* dat)								\
{																															\
	sizeType size;																											\
	readBinaryData(size,is,dat);																							\
	v.resize(size);																											\
	for(sizeType i=0;i<(sizeType)v.size();i++)																				\
		readBinaryData(v[i],is,dat);																						\
	return is;																												\
}																															\
ostream& writeVector(const vector<NAME,Eigen::aligned_allocator<NAME> >& v,ostream& os,IOData* dat)							\
{																															\
	sizeType size=(sizeType)v.size();																						\
	writeBinaryData(size,os,dat);																							\
	for(sizeType i=0;i<(sizeType)v.size();i++)																				\
		writeBinaryData(v[i],os,dat);																						\
	return os;																												\
}
IO_VEC_FIXED(Vec2f)
IO_VEC_FIXED(Vec3f)
IO_VEC_FIXED(Vec4f)
IO_VEC_FIXED(Mat2f)
IO_VEC_FIXED(Mat3f)
IO_VEC_FIXED(Mat4f)
//IO_VEC_FIXED(Quatf)
IO_VEC_FIXED(Vec2d)
IO_VEC_FIXED(Vec3d)
IO_VEC_FIXED(Vec4d)
IO_VEC_FIXED(Mat2d)
IO_VEC_FIXED(Mat3d)
IO_VEC_FIXED(Mat4d)
//IO_VEC_FIXED(Quatd)
IO_VEC_FIXED(Vec2i)
IO_VEC_FIXED(Vec3i)
IO_VEC_FIXED(Vec4i)
#undef IO_VEC_FIXED

//for file format definition
HasMagic::HasMagic(const sizeType& MAGIC):_MAGIC(MAGIC) {}
bool HasMagic::readMagic(istream& is)
{
  if(!is.good())
    return false;

  sizeType MAGIC_Test;
  readBinaryData(MAGIC_Test,is);
  if(MAGIC_Test != _MAGIC)
    return false;

  return true;
}
bool HasMagic::writeMagic(ostream& os) const
{
  if(!os.good())
    return false;

  writeBinaryData(_MAGIC,os);
  return true;
}
HasMagic& HasMagic::operator=(const HasMagic& other)
{
  ASSERT(other._MAGIC == _MAGIC)return *this;
}

//VTK compatible atomic IO
class Endianness
{
public:
  static bool isLittleEndian() {
    union u {
      unsigned long l;
      unsigned char c[sizeof(unsigned long)];
    };
    u dummy;
    dummy.l = 1;
    return dummy.c[0] == 1;
  }
  static void swap2Bytes(unsigned char* &ptr) {
    unsigned char tmp;
    tmp = ptr[0];
    ptr[0] = ptr[1];
    ptr[1] = tmp;
  }
  static void swap4Bytes(unsigned char* &ptr) {
    unsigned char tmp;
    tmp = ptr[0];
    ptr[0] = ptr[3];
    ptr[3] = tmp;
    tmp = ptr[1];
    ptr[1] = ptr[2];
    ptr[2] = tmp;
  }
  static void swap8Bytes(unsigned char* &ptr) {
    unsigned char tmp;
    tmp = ptr[0];
    ptr[0] = ptr[7];
    ptr[7] = tmp;
    tmp = ptr[1];
    ptr[1] = ptr[6];
    ptr[6] = tmp;
    tmp = ptr[2];
    ptr[2] = ptr[5];
    ptr[5] = tmp;
    tmp = ptr[3];
    ptr[3] = ptr[4];
    ptr[4] = tmp;
  }
};
template<class T2>
FORCE_INLINE void byteSwap(T2& val)
{
  sizeType n=sizeof(T2);
  unsigned char *p=(unsigned char*)&val;
  switch( n ) {
  case 1:
    return;
  case 2:
    Endianness::swap2Bytes(p);
    break;
  case 4:
    Endianness::swap4Bytes(p);
    break;
  case 8:
    Endianness::swap8Bytes(p);
    break;
  default:
    break;
  }
}
template <typename T2>
FORCE_INLINE void vtkWriteTpl(ostream& oss,T2 val)
{
  if(typeid(val) == typeid(scalarD)) {
    double valD=convert<double>()(val);
    if(Endianness::isLittleEndian())
      byteSwap(valD);
    oss.write((const char*)&valD,sizeof(double));
  } else {
    if(Endianness::isLittleEndian())
      byteSwap(val);
    oss.write((const char*)&val,sizeof(T2));
  }
}
void vtkWrite(ostream& oss,scalarF val)
{
  vtkWriteTpl(oss,val);
}
void vtkWrite(ostream& oss,scalarD val)
{
  vtkWriteTpl(oss,val);
}
void vtkWrite(ostream& oss,int val)
{
  vtkWriteTpl(oss,val);
}
void vtkWrite(ostream& oss,unsigned char val)
{
  oss.write((const char*)&val,1);
}

PRJ_END

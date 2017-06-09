#include "IOBasic.h"

PRJ_BEGIN

#define IO_BASIC(T)                                                                                             \
ostream& writeBinaryData(const T& val,ostream& os,IOData* dat){os.write((char*)&val,sizeof(T));return os;}	\
istream& readBinaryData(T& val,istream& is,IOData* dat){is.read((char*)&val,sizeof(T));return is;}
IO_BASIC(char)
IO_BASIC(unsigned char)
IO_BASIC(short)
IO_BASIC(unsigned short)
IO_BASIC(int)
IO_BASIC(unsigned int)
//IO_BASIC(scalarD)
IO_BASIC(bool)
IO_BASIC(sizeType)
#undef IO_BASIC

ostream& writeBinaryData(scalarF val,ostream& os,IOData* dat)
{
  double valD=val;
  os.write((char*)&valD,sizeof(double));
  return os;
}
istream& readBinaryData(scalarF& val,istream& is,IOData* dat)
{
  double valD;
  is.read((char*)&valD,sizeof(double));
  val=(scalarF)valD;
  return is;
}
ostream& writeBinaryData(scalarD val,ostream& os,IOData* dat)
{
  double valD=convert<double>()(val);
  os.write((char*)&valD,sizeof(double));
  return os;
}
istream& readBinaryData(scalarD& val,istream& is,IOData* dat)
{
  double valD;
  is.read((char*)&valD,sizeof(double));
  val=scalarD(valD);
  return is;
}

#define IO_FIXED(NAME,TO_TYPE)									\
ostream& writeBinaryData(const NAME& v,ostream& os,IOData* dat)					\
{												\
	sizeType d0=NAME::RowsAtCompileTime;os.write((char*)&d0,sizeof(sizeType));		\
	sizeType d1=NAME::ColsAtCompileTime;os.write((char*)&d1,sizeof(sizeType));		\
	for(sizeType r=0;r<NAME::RowsAtCompileTime;r++)						\
	for(sizeType c=0;c<NAME::ColsAtCompileTime;c++)						\
	{											\
		TO_TYPE val=convert<TO_TYPE>()(v(r,c));						\
		os.write((char*)&val,sizeof(TO_TYPE));						\
	}											\
	return os;										\
}												\
istream& readBinaryData(NAME& v,istream& is,IOData* dat)					\
{												\
	sizeType d0;is.read((char*)&d0,sizeof(sizeType));ASSERT(d0 == NAME::RowsAtCompileTime)	\
	sizeType d1;is.read((char*)&d1,sizeof(sizeType));ASSERT(d1 == NAME::ColsAtCompileTime)	\
	for(sizeType r=0;r<NAME::RowsAtCompileTime;r++)						\
	for(sizeType c=0;c<NAME::ColsAtCompileTime;c++)						\
	{											\
		TO_TYPE val;									\
		is.read((char*)&val,sizeof(TO_TYPE));						\
		v(r,c)=NAME::Scalar(val);							\
	}											\
	return is;										\
}
//double
IO_FIXED(Vec2d,double)
IO_FIXED(Vec3d,double)
IO_FIXED(Vec4d,double)
IO_FIXED(Vec6d,double)
IO_FIXED(Mat2d,double)
IO_FIXED(Mat3d,double)
IO_FIXED(Mat4d,double)
IO_FIXED(Mat6d,double)
//float is double
IO_FIXED(Vec2f,double)
IO_FIXED(Vec3f,double)
IO_FIXED(Vec4f,double)
IO_FIXED(Vec6f,double)
IO_FIXED(Mat2f,double)
IO_FIXED(Mat3f,double)
IO_FIXED(Mat4f,double)
IO_FIXED(Mat6f,double)
//sizeType
IO_FIXED(Vec2i,sizeType)
IO_FIXED(Vec3i,sizeType)
IO_FIXED(Vec4i,sizeType)
//char
IO_FIXED(Vec2c,char)
IO_FIXED(Vec3c,char)
IO_FIXED(Vec4c,char)
#undef IO_FIXED

#define IO_NON_FIXED(NAME,TO_TYPE)									\
ostream& writeBinaryData(const NAME& v,ostream& os,IOData* dat)						\
{													\
	sizeType d0=v.rows();										\
	sizeType d1=v.cols();										\
	os.write((char*)&d0,sizeof(sizeType));								\
	os.write((char*)&d1,sizeof(sizeType));								\
	if(d0 <= 1 || d1 <= 1)										\
	{												\
		if(d0*d1 > 0)os.write((char*)(v.cast<TO_TYPE>().eval().data()),sizeof(TO_TYPE)*d0*d1);	\
		return os;										\
	}												\
	Eigen::Matrix<TO_TYPE,-1,1> sub;sub.setZero(d1);						\
	for(sizeType r=0;r<d0;r++)									\
	{												\
		sub=v.row(r).cast<TO_TYPE>();								\
		if(d1 > 0)os.write((char*)(sub.data()),sizeof(TO_TYPE)*d1);				\
	}												\
	return os;											\
}													\
istream& readBinaryData(NAME& v,istream& is,IOData* dat)						\
{													\
typedef Eigen::Matrix<TO_TYPE,NAME::RowsAtCompileTime,NAME::ColsAtCompileTime> TONAME;			\
	sizeType d0;is.read((char*)&d0,sizeof(sizeType));						\
	sizeType d1;is.read((char*)&d1,sizeof(sizeType));						\
	v.resize(d0,d1);										\
	if(d0 <= 1 || d1 <= 1)										\
	{												\
		TONAME vcast(d0,d1);									\
		if(d0*d1 > 0)is.read((char*)(vcast.data()),sizeof(TO_TYPE)*d0*d1);			\
		v=vcast.cast<NAME::Scalar>();								\
		return is;										\
	}												\
	Eigen::Matrix<TO_TYPE,-1,1> sub;sub.setZero(d1);						\
	for(sizeType r=0;r<d0;r++)									\
	{												\
		if(d1 > 0)is.read((char*)(sub.data()),sizeof(TO_TYPE)*d1);				\
		v.row(r)=sub.cast<NAME::Scalar>();							\
	}												\
	return is;											\
}
IO_NON_FIXED(Rowd,double)
IO_NON_FIXED(Cold,double)
IO_NON_FIXED(Matd,double)
IO_NON_FIXED(Rowf,double)
IO_NON_FIXED(Colf,double)
IO_NON_FIXED(Matf,double)
IO_NON_FIXED(Rowi,sizeType)
IO_NON_FIXED(Coli,sizeType)
IO_NON_FIXED(Mati,sizeType)
#undef IO_NON_FIXED

#define IO_FIXED_QUAT(NAMEQ,NAMET,NAMEA,TO_TYPE)			\
ostream& writeBinaryData(const NAMEQ& v,ostream& os,IOData* dat)	\
{									\
	TO_TYPE val=convert<TO_TYPE>()(v.w());				\
	os.write((char*)&val,sizeof(TO_TYPE));				\
	val=convert<TO_TYPE>()(v.x());					\
	os.write((char*)&val,sizeof(TO_TYPE));				\
	val=convert<TO_TYPE>()(v.y());					\
	os.write((char*)&val,sizeof(TO_TYPE));				\
	val=convert<TO_TYPE>()(v.z());					\
	os.write((char*)&val,sizeof(TO_TYPE));				\
	return os;							\
}									\
istream& readBinaryData(NAMEQ& v,istream& is,IOData* dat)		\
{									\
	TO_TYPE val;							\
	is.read((char*)&val,sizeof(TO_TYPE));				\
	v.w()=NAMEQ::Scalar(val);					\
	is.read((char*)&val,sizeof(TO_TYPE));				\
	v.x()=NAMEQ::Scalar(val);					\
	is.read((char*)&val,sizeof(TO_TYPE));				\
	v.y()=NAMEQ::Scalar(val);					\
	is.read((char*)&val,sizeof(TO_TYPE));				\
	v.z()=NAMEQ::Scalar(val);					\
	return is;							\
}									\
ostream& writeBinaryData(const NAMET& v,ostream& os,IOData* dat)	\
{									\
	TO_TYPE val=convert<TO_TYPE>()(v.x());				\
	os.write((char*)&val,sizeof(TO_TYPE));				\
	val=convert<TO_TYPE>()(v.y());					\
	os.write((char*)&val,sizeof(TO_TYPE));				\
	val=convert<TO_TYPE>()(v.z());					\
	os.write((char*)&val,sizeof(TO_TYPE));				\
	return os;							\
}									\
istream& readBinaryData(NAMET& v,istream& is,IOData* dat)		\
{									\
	TO_TYPE val;							\
	is.read((char*)&val,sizeof(TO_TYPE));				\
	v.x()=NAMEQ::Scalar(val);					\
	is.read((char*)&val,sizeof(TO_TYPE));				\
	v.y()=NAMEQ::Scalar(val);					\
	is.read((char*)&val,sizeof(TO_TYPE));				\
	v.z()=NAMEQ::Scalar(val);					\
	return is;							\
}									\
ostream& writeBinaryData(const NAMEA& v,ostream& os,IOData* dat)	\
{                                                                       \
  Eigen::Matrix<NAMEA::Scalar,3,3> l=v.linear();                        \
  Eigen::Matrix<NAMEA::Scalar,3,1> t=v.translation();                   \
  writeBinaryData(l,os,dat);                                            \
  writeBinaryData(t,os,dat);                                            \
  return os;                                                            \
}                                                                       \
istream& readBinaryData(NAMEA& v,istream& is,IOData* dat)		\
{									\
  Eigen::Matrix<NAMEA::Scalar,3,3> l;                                   \
  Eigen::Matrix<NAMEA::Scalar,3,1> t;                                   \
  readBinaryData(l,is,dat);                                             \
  readBinaryData(t,is,dat);                                             \
  v.linear()=l;                                                         \
  v.translation()=t;                                                    \
  return is;                                                            \
}
IO_FIXED_QUAT(Quatd,Transd,Affined,double)
IO_FIXED_QUAT(Quatf,Transf,Affinef,double)
#undef IO_FIXED_QUAT

#define IO_BB(NAME)							\
ostream& writeBinaryData(const BBox<NAME>& b,ostream& os,IOData* dat)	\
{									\
	writeBinaryData(b._minC,os,dat);				\
	writeBinaryData(b._maxC,os,dat);				\
	return os;							\
}									\
istream& readBinaryData(BBox<NAME>& b,istream& is,IOData* dat)		\
{									\
	readBinaryData(b._minC,is,dat);					\
	readBinaryData(b._maxC,is,dat);					\
	return is;							\
}
IO_BB(scalarD)
IO_BB(scalarF)
#undef IO_BB

PRJ_END

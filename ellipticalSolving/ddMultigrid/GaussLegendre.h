#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H

#include "CommonFile/MathBasic.h"
#include "CommonFile/GridOp.h"
#include "CommonFile/CollisionDetection.h"

PRJ_BEGIN

class GaussLegendre
{
public:
  static const scalarD WEIGHT[64][64];
  static const scalarD POINT[64][64];
};
class DunavantTriangle
{
public:
  static const scalarD DUNAVANT_WEIGHT[20][79];
  static const scalarD DUNAVANT_POINT[20][158];
  static const sizeType DUNAVANT_NR[20];
};
template <typename RET_TYPE,typename T>
class GaussLegendreIntegral : private GaussLegendre, private DunavantTriangle
{
public:
  typedef typename ScalarUtil<T>::ScalarVec3 PT3;
  template <typename FUNC_TYPE>
  static RET_TYPE integrate1D(FUNC_TYPE& f,sizeType deg) {
    ASSERT(deg >= 0 && deg < 64)
    const scalarD* weights=WEIGHT[deg];
    const scalarD* points=POINT[deg];

    RET_TYPE ret=Zero<RET_TYPE>::value();
    for(sizeType i=0; i<=deg; i++)
      ret+=f(PT3(T(points[i]),0.0f,0.0f))*T(weights[i]);
    return ret;
  }
  template <typename FUNC_TYPE>
  static RET_TYPE integrate2D(FUNC_TYPE& f,sizeType deg) {
    ASSERT(deg >= 0 && deg < 64)
    const scalarD* weights=WEIGHT[deg];
    const scalarD* points=POINT[deg];

    RET_TYPE ret=Zero<RET_TYPE>::value();
    for(sizeType i=0; i<=deg; i++)
      for(sizeType j=0; j<=deg; j++)
        ret+=f(PT3(T(points[i]),T(points[j]),0.0f))*T(weights[i]*weights[j]);
    return ret;
  }
  template <typename FUNC_TYPE>
  static RET_TYPE integrate3D(FUNC_TYPE& f,sizeType deg) {
    ASSERT(deg >= 0 && deg < 64)
    const scalarD* weights=WEIGHT[deg];
    const scalarD* points=POINT[deg];

    RET_TYPE ret=Zero<RET_TYPE>::value();
    for(sizeType i=0; i<=deg; i++)
      for(sizeType j=0; j<=deg; j++)
        for(sizeType k=0; k<=deg; k++)
          ret+=f(PT3(T(points[i]),T(points[j]),T(points[k])))*T(weights[i]*weights[j]*weights[k]);
    return ret;
  }
  template <typename FUNC_TYPE>
  static RET_TYPE integrate1DBF(FUNC_TYPE& f,const PT3& dx) {
    sizeType nrS=(sizeType)std::ceil(2.0f/dx[0]);
    T lenS=(scalarD)2.0f/(scalarD)nrS;

    RET_TYPE ret=Zero<RET_TYPE>::value();
    for(sizeType i=0; i<nrS; i++)
      ret+=f(PT3(scalarD(i+0.5f)/scalarD(nrS)-1.0f,0.0f,0.0f))*lenS;
    return ret;
  }
  template <typename FUNC_TYPE>
  static RET_TYPE integrate2DBF(FUNC_TYPE& f,const PT3& dx) {
    sizeType nrS=(sizeType)std::ceil(2.0f/dx[0]);
    T lenS=(scalarD)2.0f/(scalarD)nrS;
    lenS=lenS*lenS;

    RET_TYPE ret=Zero<RET_TYPE>::value();
    for(sizeType i=0; i<nrS; i++)
      for(sizeType j=0; j<nrS; j++)
        ret+=f(PT3(scalarD(i+0.5f)/scalarD(nrS)-1.0f,
                   scalarD(j+0.5f)/scalarD(nrS)-1.0f,0.0f))*lenS;
    return ret;
  }
  template <typename FUNC_TYPE>
  static RET_TYPE integrate3DBF(FUNC_TYPE& f,const PT3& dx) {
    sizeType nrS=(sizeType)std::ceil(2.0f/dx[0]);
    T lenS=(scalarD)2.0f/(scalarD)nrS;
    lenS=lenS*lenS*lenS;

    RET_TYPE ret=Zero<RET_TYPE>::value();
    for(sizeType i=0; i<nrS; i++)
      for(sizeType j=0; j<nrS; j++)
        for(sizeType k=0; k<nrS; k++)
          ret+=f(PT3(scalarD(i+0.5f)/scalarD(nrS)-1.0f,
                     scalarD(j+0.5f)/scalarD(nrS)-1.0f,
                     scalarD(k+0.5f)/scalarD(nrS)-1.0f))*lenS;
    return ret;
  }
  static void getLineSample(const PT3& a,const PT3& b,sizeType deg,std::vector<PT3,Eigen::aligned_allocator<PT3> >& eps,std::vector<T>& ews) {
    const PT3 ctr=(b+a)/2.0f;
    const PT3 v0=(b-a)/2.0f;
    const T w=v0.norm();

    ASSERT(deg >= 0 && deg < 64)
    const scalarD* weights=WEIGHT[deg];
    const scalarD* points=POINT[deg];

    for(sizeType i=0; i<=deg; i++) {
      PT3 ep=v0*T(points[i])+ctr;
      T ew=T(weights[i])*w;
      eps.push_back(ep);
      ews.push_back(ew);
    }
  }
  template <typename FUNC_TYPE>
  static RET_TYPE integrateLine(const PT3& a,const PT3& b,FUNC_TYPE& f,sizeType deg) {
    std::vector<PT3,Eigen::aligned_allocator<PT3> > eps;
    std::vector<T> ews;
    getLineSample(a,b,deg,eps,ews);

    RET_TYPE ret=Zero<RET_TYPE>::value();
    for(sizeType i=0; i<(sizeType)eps.size(); i++)
      ret+=f(eps[i])*ews[i];
    return ret;
  }
  //for this transformation see: "Gauss Legendre quadrature over a triangle" by H. T. RATHOD1 K. V. NAGARAJA2 B. VENKATESUDU3 AND N. L. RAMESH
  static void getTriSampleSimple(const PT3& a,const PT3& b,const PT3& c,sizeType deg,std::vector<PT3,Eigen::aligned_allocator<PT3> >& eps,std::vector<T>& ews) {
    const PT3 v0=b-a;
    const PT3 v1=c-a;
    const T w=v1.cross(v0).norm();

    ASSERT(deg >= 0 && deg < 64)
    const scalarD* weights=WEIGHT[deg];
    const scalarD* points=POINT[deg];

    for(sizeType i=0; i<=deg; i++)
      for(sizeType j=0; j<=deg; j++) {
        scalar c0=T(1.0f+points[i])/2.0f;
        scalar c1=T(1.0f-points[i])*T(1.0f+points[j])/4.0f;
        PT3 ep=v0*c0+v1*c1+a;
        T ew=T(weights[i]*weights[j])*w*(1.0f-points[i])/8.0f;
        eps.push_back(ep);
        ews.push_back(ew);
      }
  }
  //the famous Dunavant Quadrature rule see: "High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle" by David Dunavant
  static void getTriSample(const PT3& a,const PT3& b,const PT3& c,sizeType deg,std::vector<PT3,Eigen::aligned_allocator<PT3> >& eps,std::vector<T>& ews) {
    const PT3 v0=b-a;
    const PT3 v1=c-a;
    const T w=v1.cross(v0).norm()/2.0f;

    ASSERT(deg >= 0 && deg < 20)
    const scalarD* weights=DUNAVANT_WEIGHT[deg];
    const scalarD* points=DUNAVANT_POINT[deg];
    sizeType nrPoint=DUNAVANT_NR[deg];

    for(sizeType i=0; i<nrPoint; i++) {
      T c0=T(points[i*2]);
      T c1=T(points[i*2+1]);
      PT3 ep=a*c0+b*c1+c*(1.0f-c0-c1);
      T ew=T(weights[i])*w;
      eps.push_back(ep);
      ews.push_back(ew);
    }
  }
  template <typename FUNC_TYPE>
  static RET_TYPE integrateTri(const PT3& a,const PT3& b,const PT3& c,FUNC_TYPE& f,sizeType deg) {
    std::vector<PT3,Eigen::aligned_allocator<PT3> > eps;
    std::vector<T> ews;
    getTriSample(a,b,c,deg,eps,ews);

    RET_TYPE ret=Zero<RET_TYPE>::value();
    for(sizeType i=0; i<(sizeType)eps.size(); i++)
      ret+=f(eps[i])*ews[i];
    return ret;
  }
  template <typename FUNC_TYPE>
  static RET_TYPE integrateLineBF(const PT3& a,const PT3& b,FUNC_TYPE& f,const PT3& dx) {
    LineSegTpl<T> l(a,b);
    T len=l.length();
    if(len < dx[0])
      return len*f(l.masscenter());
    else {
      PT3 mid=l.masscenter();
      return integrateLineBF(a,mid,f,dx)+integrateLineBF(mid,b,f,dx);
    }
  }
  template <typename FUNC_TYPE>
  static RET_TYPE integrateTriBF(const PT3& a,const PT3& b,const PT3& c,FUNC_TYPE& f,const PT3& dx) {
    TriangleTpl<T> tri(a,b,c);
    T area=tri.area();
    if(area < dx[0]*dx[1]) {
      return area*f(tri.masscenter());
    } else {
      PT3 midAB=(a+b)/2.0f;
      PT3 midBC=(b+c)/2.0f;
      PT3 midAC=(a+c)/2.0f;
      return integrateTriBF(a,midAB,midAC,f,dx)+
             integrateTriBF(midAB,b,midBC,f,dx)+
             integrateTriBF(midAB,midBC,midAC,f,dx)+
             integrateTriBF(midAC,midBC,c,f,dx);
    }
  }
  //query
  static sizeType nrP3D(sizeType deg) {
    return (deg+1)*(deg+1)*(deg+1);
  }
  static scalarD weight3D(sizeType deg,sizeType id) {
    sizeType off=(deg+1)*(deg+1);
    return WEIGHT[deg][id/off]*
           WEIGHT[deg][(id%off)/(deg+1)]*
           WEIGHT[deg][id%(deg+1)];
  }
  static Vec3d point3D(sizeType deg,sizeType id) {
    sizeType off=(deg+1)*(deg+1);
    return Vec3d(POINT[deg][id/off],
                 POINT[deg][(id%off)/(deg+1)],
                 POINT[deg][id%(deg+1)]);
  }
};

PRJ_END

#endif

#include "DHMUtil.h"
#include "CommonFile/IO.h"
#include "CommonFile/GridBasic.h"
#include "GaussLegendre.h"

USE_PRJ_NAMESPACE

//make simple mesh
void DHMMeshMaker::makeCube(POSES& poses,CELLIDS& cells,const BBox<scalar>& bb,const Vec3i& nrCell)
{
  Vec3i nrP=nrCell+Vec3i::Ones();
  Vec3i stride(nrP[1]*nrP[2],nrP[2],1);
  Vec3 cSize=(bb.getExtent().array()/nrCell.cast<scalar>().array()).matrix();

  poses.resize(nrP.prod());
  for(sizeType x=0; x<nrP[0]; x++)
    for(sizeType y=0; y<nrP[1]; y++)
      for(sizeType z=0; z<nrP[2]; z++) {
        Vec3i id(x,y,z);
        poses[id.dot(stride)]=bb._minC+(cSize.array()*id.cast<scalar>().array()).matrix();
      }

  CELLID cid;
  Vec3i cstride(nrCell[1]*nrCell[2],nrCell[2],1);
  cells.resize(nrCell.prod());
  for(sizeType x=0; x<nrCell[0]; x++)
    for(sizeType y=0; y<nrCell[1]; y++)
      for(sizeType z=0; z<nrCell[2]; z++) {
        Vec3i id(x,y,z);
        for(char v=0; v<8; v++)
          cid[v]=stride.dot(Vec3i(id+Vec3i(v&1?1:0,v&2?1:0,v&4?1:0)));
        cells[id.dot(cstride)]=cid;
      }
}
void DHMMeshMaker::makeCube(const std::string& path,const BBox<scalar>& bb,const Vec3i& nrCell)
{
  POSES poses;
  CELLIDS cells;
  makeCube(poses,cells,bb,nrCell);
  VTKWriter<scalar> os("cube",path,false);
  os.appendPoints(poses.begin(),poses.end());
  os.appendCells(cells.begin(),cells.end(),VTKWriter<scalar>::VOXEL);
}
void DHMMeshMaker::makeSphere(const std::string& path,const scalar rad,sizeType nrPad,const Vec3i& nrCell)
{
  typedef boost::unordered_map<Vec4i,boost::shared_ptr<DHMFace>,Hash> FMAP;
  typedef FMAP::const_iterator CFITER;
  typedef FMAP::iterator FITER;

  //inner cell
  POSES poses;
  CELLIDS cells;
  scalar radI=rad*0.5f/sqrt(3.0f);
  BBox<scalar> bb(Vec3::Constant(-radI),Vec3::Constant(radI));
  makeCube(poses,cells,bb,nrCell);

  //build face
  FMAP faceMap;
  sizeType nrC=(sizeType)cells.size();
  for(sizeType i=0; i<nrC; i++) {
    for(char f=0; f<6; f++) {
      Vec4i fid=DHMCell::getFaceId(f);
      for(char v=0; v<4; v++)
        fid[v]=cells[i][fid[v]];
      sort4(fid[0],fid[1],fid[2],fid[3]);
      FITER it=faceMap.find(fid);
      if(it == faceMap.end()) {
        faceMap[fid].reset(new DHMFace);
        faceMap[fid]->_cells[0].reset(new DHMCell(i));
      } else {
        ASSERT(!faceMap[fid]->_cells[1])
        faceMap[fid]->_cells[1].reset(new DHMCell(i));
      }
    }
  }

  typedef boost::unordered_map<Vec2i,vector<sizeType>,Hash> SPMAP;
  typedef SPMAP::const_iterator CPITER;

  //build padding vertex
  CELLID cid;
  SPMAP padMap;
  for(CFITER beg=faceMap.begin(),end=faceMap.end(); beg!=end; beg++)
    if(!beg->second->_cells[1]) {
      //add vertices
      for(char v=0; v<4; v++) {
        CPITER it=padMap.find(Vec2i::Constant(beg->first[v]));
        if(it == padMap.end()) {
          Vec3 pt=poses[beg->first[v]];
          Vec3 dir=pt*rad/pt.norm()-pt;

          vector<sizeType>& pVert=padMap[Vec2i::Constant(beg->first[v])];
          pVert.push_back(beg->first[v]);
          for(sizeType p=1; p<=nrPad; p++) {
            pVert.push_back((sizeType)poses.size());
            poses.push_back(pt+dir*(scalar)p/(scalar)nrPad);
          }
        }
      }
      //add cells
      for(sizeType p=1; p<=nrPad; p++) {
        for(char v=0; v<4; v++) {
          const vector<sizeType>& pVert=padMap[Vec2i::Constant(beg->first[v])];
          cid[v*2+0]=pVert[p-1];
          cid[v*2+1]=pVert[p];
        }
        cells.push_back(cid);
      }
    }
  //write
  VTKWriter<scalar> os("cube",path,false);
  os.appendPoints(poses.begin(),poses.end());
  os.appendCells(cells.begin(),cells.end(),VTKWriter<scalar>::VOXEL);
}
void DHMMeshMaker::makeDeformedCube(const std::string& path,const Vec3i& nrCell,scalar rot)
{
  POSES poses;
  CELLIDS cells;
  makeCube(poses,cells,BBox<scalar>(-Vec3::Ones(),Vec3::Ones()),nrCell);
  for(sizeType i=0; i<(sizeType)poses.size(); i++) {
    Vec3& p=poses[i];
    p=makeRotation<scalar>(Vec3(0.0f,0.0f,p[2])*rot)*p;
  }
  VTKWriter<scalar> os("cube",path,false);
  os.appendPoints(poses.begin(),poses.end());
  os.appendCells(cells.begin(),cells.end(),VTKWriter<scalar>::VOXEL);
}
//mapping function and deformation gardient kernel
void DHMMapping::debugMapping()
{
  Eigen::Matrix<scalarD,8,1> Ms,MDs[3];
  Eigen::Matrix<scalarD,3,8> Js,Vs;
  for(sizeType i=0; i<100; i++) {
    //debug gradient
#define DELTA 1E-7f
    Vec3d c=Vec3d::Random(),cD[3]= {
      c+Vec3d::Unit(0)*DELTA,
      c+Vec3d::Unit(1)*DELTA,
      c+Vec3d::Unit(2)*DELTA,
    };
    Vs.setRandom();
    MStencil(Ms,c);
    MStencil(MDs[0],cD[0]);
    MStencil(MDs[1],cD[1]);
    MStencil(MDs[2],cD[2]);
    JStencil(Js,c);

    Vec3d M=Vec3d::Zero();
    Mat3d MD=Mat3d::Zero(),J=Mat3d::Zero();
    for(char v=0; v<8; v++) {
      M+=Vs.col(v)*Ms[v];
      MD.col(0)+=Vs.col(v)*MDs[0][v];
      MD.col(1)+=Vs.col(v)*MDs[1][v];
      MD.col(2)+=Vs.col(v)*MDs[2][v];
      J.col(0)+=Vs.col(v)*Js(0,v);
      J.col(1)+=Vs.col(v)*Js(1,v);
      J.col(2)+=Vs.col(v)*Js(2,v);
    }
    MD.col(0)=(MD.col(0)-M)/DELTA;
    MD.col(1)=(MD.col(1)-M)/DELTA;
    MD.col(2)=(MD.col(2)-M)/DELTA;
    INFOV("Err: %f %f",J.norm(),(J-MD).norm())
#undef DELTA

    //debug DFDX
    Eigen::Matrix<scalarD,3,8> ps,ps2;
    for(char v=0; v<8; v++) {
      ps.col(v).setRandom();
      ps2.col(v).setRandom();
    }

    Mat3d F,F2;
    JP(Js,ps,F);
    JP(Js,ps2,F2);
    Eigen::Matrix<scalarD,9,24> DFDX;
    calcDFDX(Js,F.inverse().eval(),DFDX);

    Mat3d FF2=F2*F.inverse().eval();
    Eigen::Matrix<scalarD,9,1> FF22=
      DFDX*Eigen::Map<Eigen::Matrix<scalarD,24,1> >(ps2.data());
    INFOV("ErrDFDX: %f %f",FF2.norm(),(FF2-Eigen::Map<Mat3d>(FF22.data())).norm());
  }
}
//uniform validity check
Eigen::Matrix<scalarD,27,1> DHMUniformValidityBound::bernstein(const Vec3d& pos)
{
  Eigen::Matrix<scalarD,27,1> ret;
  char id=0;
  scalarD f[3]= {1,2,1};
  scalarD B[3][3];
  for(char pd=0; pd<3; pd++)
    for(char d=0; d<=2; d++)
      B[pd][d]=f[d]*std::pow((1.0f-pos[pd]),2-d)*std::pow(pos[pd],d);
  for(char d1=0; d1<=2; d1++)
    for(char d2=0; d2<=2; d2++)
      for(char d3=0; d3<=2; d3++)
        ret[id++]=B[0][d1]*B[1][d2]*B[2][d3];
  return ret;
}
Eigen::Matrix<scalarD,27,1> DHMUniformValidityBound::calcValidBound(const DHMCell& cell)
{
  const scalarD* POINT=GaussLegendre::POINT[2];
  char id=0;
  Vec3d pts[27];
  for(char x=0; x<=2; x++)
    for(char y=0; y<=2; y++)
      for(char z=0; z<=2; z++)
        pts[id++]=(Vec3d(POINT[x],POINT[y],POINT[z])+Vec3d::Ones())*0.5f;

  Mat3d J;
  Eigen::Matrix<scalarD,3,8> Js;
  Eigen::Matrix<scalarD,27,1> RHS;
  Eigen::Matrix<scalarD,27,27> LHS;
  for(char v=0; v<27; v++) {
    DHMMapping::JStencil(Js,pts[v]);
    DHMMapping::J(Js,cell,J);
    RHS[v]=J.determinant();
    LHS.row(v)=bernstein(pts[v]).transpose();
  }
  return LHS.fullPivLu().solve(RHS);
}
bool DHMUniformValidityBound::isValid(const DHMCell& cell)
{
  Eigen::Matrix<scalarD,27,1> ret=calcValidBound(cell);
  for(char d=0; d<27; d++)
    if(ret[d] <= 0.0f)
      return false;
  return true;
}
void DHMUniformValidityBound::debugCorrect()
{
  //initialize a random cell
  DHMCell c(0);
  for(sizeType v=0; v<8; v++) {
    c._verts[v].reset(new DHMVertex(v));
    c._verts[v]->_pos.setRandom();
  }

  //check positive bernstein
  for(char v=0; v<100; v++) {
    Vec3d pos=Vec3d::Random();
    Eigen::Matrix<scalarD,27,1> ret=bernstein((pos+Vec3d::Ones())*0.5f);
    for(char d=0; d<27; d++)
      ASSERT(ret[d] >= 0.0f)
    }

  //check jacobian determinant regeneration
  Mat3d J;
  Eigen::Matrix<scalarD,3,8> Js;
  Eigen::Matrix<scalarD,27,1> coeff=calcValidBound(c);
  for(sizeType v=0; v<100; v++) {
    Vec3d pos=(Vec3d::Random()+Vec3d::Ones())*0.5f;
    scalarD d1=coeff.dot(bernstein(pos));

    DHMMapping::JStencil(Js,pos);
    DHMMapping::J(Js,c,J);
    scalarD d2=J.determinant();

    INFOV("Jacobian Determinant: %f %f",d1,d2)
  }
}
void DHMUniformValidityBound::debugTightness()
{
  //initialize uniform cell
  DHMCell c(0);
  Vec3 pos0[8];
  Vec3 delta[8];
  for(char v=0; v<8; v++) {
    c._verts[v].reset(new DHMVertex(v));
    pos0[v]=Vec3(v&1?1.0f:0.0f,v&2?1.0f:0.0f,v&4?1.0f:0.0f);
  }

  //search
#define ASSIGN(D)	\
for(char v=0;v<8;v++)	\
c._verts[v]->_pos=pos0[v]+delta[v]*(D);
  for(sizeType i=0; i<100; i++) {
    //search at random direction
    INFOV("Search Iter %ld",i)
    scalar a=0.0f;
    scalar b=1.0f;
    while(true) {
      for(char v=0; v<8; v++)
        delta[v]=Vec3::Random();
      for(b=1.0f; b<1E4f; b*=2.0f) {
        ASSIGN(b)
        if(!isValid(c))
          break;
      }
      if(!isValid(c))
        break;
    }

    //binary search
    while(fabs(a-b) > 0.01f) {
      scalar mid=(a+b)*0.5f;
      ASSIGN(mid)
      if(isValid(c))
        a=mid;
      else b=mid;
    }

    //write nearest invalid cell
    ASSIGN(a)
    ostringstream oss;
    oss << "./cell" << i << ".vtk";
    c.writeVTK(oss.str());
  }
}

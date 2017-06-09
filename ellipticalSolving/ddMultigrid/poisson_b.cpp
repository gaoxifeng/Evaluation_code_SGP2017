#include "DHMMesh.h"
#include "DHMLinearSolver.h"
#include "CommonFile/solvers/DACGSolver.h"
#include "CommonFile/solvers/AINVPreconditioner.h"
#include "CommonFile/solvers/PMinresQLP.h"
#include "CommonFile/solvers/BlockedSolver.h"
//#include "CommonFile/solvers/AGMGSolver.h"

USE_PRJ_NAMESPACE

template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct DivFreeKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE> {
public:
  typedef typename KERNEL_TYPE::Vec Vec;
  DivFreeKrylovMatrix(const FixedSparseMatrix<T,KERNEL_TYPE> &sB) {
    sB.toEigen(_sB);
    _sol.compute(_sB*_sB.transpose());
  }
  virtual void multiply(const Vec& b,Vec& out) const {
    Vec lambda=_sol.solve(_sB*b);
    out=b-_sB.transpose()*lambda;
  }
  virtual sizeType n() const {
    return _sB.cols();
  }
  Eigen::SparseMatrix<T,0,sizeType> _sB;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<T,0,sizeType> > _sol;
};

static Mat3 cross(const Vec3& x){
  Mat3 ret=Mat3::Zero();
  ret(1,0)=x[2];
  ret(0,1)=-x[2];

  ret(2,0)=-x[1];
  ret(0,2)=x[1];

  ret(2,1)=x[0];
  ret(1,2)=-x[0];
  return ret;
}
static void writeMatrix(const DHMLinearSystem::SMAT& lhs,const string& path)
{
  boost::filesystem::ofstream os(path);
  os << lhs.rows() << " " << lhs.cols() << endl;
  for(sizeType i=0; i<lhs.rows(); i++)
    for(ConstSMIterator<scalarD> beg=lhs.begin(i),end=lhs.end(i); beg!=end; ++beg)
      os << beg.row() << " " << beg.col() << " " << *beg << endl;
}
template <typename T,typename Vec>
static T specRad(const FixedSparseMatrix<T,Kernel<T> >& G,Vec* ev,T eps,boost::shared_ptr<KrylovMatrix<T> > kry)
{
  T delta;
  Vec tmp,tmpOut;
  tmp.resize(G.rows());
  tmpOut.resize(G.rows());
  tmp.setRandom();
  tmp.normalize();

  //power method
  for(sizeType iter=0;; iter++) {
    G.multiply(tmp,tmpOut);
    T normTmpOut=tmpOut.norm();
    if(normTmpOut < ScalarUtil<T>::scalar_eps) {
      if(ev)*ev=tmp;
      return ScalarUtil<T>::scalar_eps;
    }

    if(kry) {
      Vec tmpOut2;
      kry->multiply(tmpOut,tmpOut2);
      tmpOut=tmpOut2;
    }
    tmpOut/=normTmpOut;
    delta=(tmpOut-tmp).norm();
    INFOV("Power Iter %ld Err: %f, SpecRad: %f",iter,delta,normTmpOut)
    if(delta <= eps) {
      if(ev)*ev=tmp;
      return normTmpOut;
    }
    tmp=tmpOut;
  }
}

int mainModalAnalysis(int argc,char** argv)
{
  //read mesh
  if(argc < 2) {
    INFO("Usage: exe [mesh].vtk")
    return -1;
  }
  DHMMesh mesh;
  mesh.readVTK(argv[1]);

  //build modal analysis
  DHMEnergyPool EC;
  for(sizeType i=0; i<mesh.nrCell(); i++)
    //you can adjust the two coefficients
    EC.addVDEnergy(boost::shared_ptr<DHMElasticEnergy>(new DHMElasticEnergy(mesh.getCell(i),1000.0f,1000.0f)));
  //EC.debugEnergy();
  FixedSparseMatrix<scalarD,Kernel<scalarD> > lhs;
  lhs.resize(mesh.nrVert()*3,mesh.nrVert()*3);
  EC.buildHessian(lhs);
  {
    DHMLinearSystem::SMAT lhsT;
    lhs.transpose(lhsT);
    lhs.add(lhsT);
    lhs.mul(0.5f);
  }

  //write matrix
  writeMatrix(lhs,"./modalAnalysis.txt");

  //get minimal eigenvectors
  {
    ASSERT(lhs.isSymmetric())
    DACGSolver<scalarD,Kernel<scalarD>,SymAINVPreconSolver<scalarD> > dacg;
    dacg.setCallback(boost::shared_ptr<Callback<scalarD,Kernel<scalarD> > >(new Callback<scalarD,Kernel<scalarD> >));
    dacg.setSolverParameters(1E-20f,1E-8f,100000);

    Matd U,U0;
    Cold lambda;
    dacg.setA(lhs,true);
    U0.setZero(lhs.rows(),6);
    for(sizeType i=0; i<lhs.rows(); i+=3) {
      U0.block<3,3>(i,0).setIdentity();
      U0.block<3,3>(i,3)=cross(mesh.getVert(i/3)._pos).cast<scalarD>();
    }
    for(sizeType i=0;i<6;i++){
      Cold Au;
      Au.setZero(lhs.rows());
      lhs.multiply(U0.col(i),Au);
      cout << "Kernel Residue: " << Au.norm() << endl;
    }
    dacg.setU0(U0);
    dacg.solve(16,lambda,U);
    boost::filesystem::ofstream os("./modalAnalysisEigenMin.txt");
    os << lambda << std::endl;
    os << U << std::endl;
  }

  //get maximum eigenvectors
  {
    Cold v;
    scalarD ev=specRad<scalarD,Cold>(lhs,&v,1E-7f,boost::shared_ptr<KrylovMatrix<scalarD> >((KrylovMatrix<scalarD>*)NULL));
    boost::filesystem::ofstream os("./modalAnalysisEigenMax.txt");
    os << ev << std::endl;
    os << v << std::endl;
    cout << ev <<std::endl;
  }
  return 0;
}
int mainPoisson(int argc,char** argv)
{
  //read mesh
  if(argc < 2) {
    INFO("Usage: exe [mesh].vtk [boundary_value_file]")
    return -1;
  }
  DHMMesh mesh;
  mesh.readVTK(argv[1]);

  //build poisson
  DHMCVField vel(mesh);
  DHMPressureSolver sol(mesh,vel);
  sol.precomputeLHS();
  DHMLinearSystem::SMAT lhs=sol.getSystem()._lhs;
  {
    DHMLinearSystem::SMAT lhsT;
    lhs.transpose(lhsT);
    lhs.add(lhsT);
    lhs.mul(0.5f);
  }

  //write matrix
  writeMatrix(lhs,"./poisson.txt");

  //input boundary, solve and output
  if(argc >= 3) {
    boost::filesystem::ifstream is("./rhs.txt");
    boost::filesystem::ifstream is2("./groundTruth.txt");
    Cold rhs=Cold::Zero(mesh.nrVert()),gt=rhs;
    for(sizeType i=0; i<rhs.size(); i++) {
      is >> rhs[i];
      is2 >> gt[i];
    }
    //rhs.setOnes();

    vector<sizeType> cv;
    for(sizeType i=0; i<mesh.nrVert(); i++)
      if(mesh.getVert(i)._isSurface)
        cv.push_back(mesh.getVert(i)._index);
    sizeType nrC=(sizeType)cv.size();

    Cold RHS=Cold::Zero(mesh.nrVert()+nrC),RES;
    FixedSparseMatrix<scalarD,Kernel<scalarD> > C;
    C.resize(nrC,lhs.cols());
    for(sizeType i=0; i<nrC; i++) {
      C.addToElement(i,cv[i],1);
      RHS[mesh.nrVert()+i]=rhs[cv[i]];
    }

    boost::shared_ptr<KKTKrylovMatrix<scalarD> > kkt(new KKTKrylovMatrix<scalarD>(lhs,C));
    PMINRESSolverQLP<scalarD> sol;
    sol.setSolverParameters(1E-8f,10000);
    sol.setKrylovMatrix(kkt);
    sol.solve(RHS,RES=RHS);
    boost::filesystem::ofstream os("./out.txt");
    for(sizeType i=0; i<mesh.nrVert(); i++)
      os << RES[i] << "\n";
  }

  //get minimal eigenvectors
  {
    ASSERT(lhs.isSymmetric())
    DACGSolver<scalarD,Kernel<scalarD>,SymAINVPreconSolver<scalarD> > dacg;
    dacg.setCallback(boost::shared_ptr<Callback<scalarD,Kernel<scalarD> > >(new Callback<scalarD,Kernel<scalarD> >));
    dacg.setSolverParameters(1E-20f,1E-8f,1000);

    Matd U,U0;
    Cold lambda;
    dacg.setA(lhs,true);
    U0.setOnes(lhs.rows(),1);
    dacg.setU0(U0);
    dacg.solve(1,lambda,U);
    boost::filesystem::ofstream os("./poissonEigenMin.txt");
    os << lambda << std::endl;
    os << U << std::endl;
  }

  //get maximum eigenvectors
  {
    Cold v;
    scalarD ev=specRad<scalarD,Cold>(lhs,&v,1E-7f,boost::shared_ptr<KrylovMatrix<scalarD> >((KrylovMatrix<scalarD>*)NULL));
    boost::filesystem::ofstream os("./poissonEigenMax.txt");
    os << ev << std::endl;
    os << v << std::endl;
  }
  return 0;
}
int mainStokes(int argc,char** argv)
{
  //read mesh
  if(argc < 2) {
    INFO("Usage: exe [mesh].vtk")
    return -1;
  }
  DHMMesh mesh;
  mesh.readVTK(argv[1]);

  //build stokes
  DHMCVField vel(mesh);
  DHMPressureSolver sol(mesh,vel);
  sol.precomputeLHS();
  DHMLinearSystem::SMAT lhs=sol.getSystem()._lhs,sA,sB;
  {
    DHMLinearSystem::SMAT lhsT;
    lhs.transpose(lhsT);
    lhs.add(lhsT);
    lhs.mul(0.5f);
    lhs.kroneckerProd(Mat3d::Identity(),sA);
    sol.precomputeDiv(sB);
  }

  //write matrix
  writeMatrix(sA,"./stokesA.txt");
  writeMatrix(sB,"./stokesB.txt");

  //get minimal eigenvectors
  {
    typedef NoPreconSolver<scalarD> PRECON;//BlockedSolver<3,SymAINVPreconSolver<scalarD> > PRECON;
    ASSERT(sA.isSymmetric())
    DACGSolver<scalarD,Kernel<scalarD>,PRECON > dacg;
    dacg.setCallback(boost::shared_ptr<Callback<scalarD,Kernel<scalarD> > >(new Callback<scalarD,Kernel<scalarD> >));
    dacg.setSolverParameters(1E-20f,1E-6f,10000);
    boost::shared_ptr<KrylovMatrix<scalarD> > kry(new BlockedKrylovMatrix<3,scalarD>(lhs));

    Matd U,U0;
    Cold lambda;
    U0.resize(lhs.rows()*3,3);
    for(sizeType i=0; i<lhs.rows(); i++)
      U0.block<3,3>(i*3,0).setIdentity();
    //((PRECON*)(dacg.getPre()))->setKrylovMatrix(kry);
    dacg.setKrylovA(kry);
    dacg.setKrylovU0(boost::shared_ptr<KrylovMatrix<scalarD> >(new DivFreeKrylovMatrix<scalarD>(sB)));
    dacg.setU0(U0);
    dacg.solve(1,lambda,U);

    boost::filesystem::ofstream os("./stokesEigenMin.txt");
    os << lambda << std::endl;
    os << U << std::endl;
  }

  //get maximum eigenvectors
  {
    Cold v;
    scalarD ev=specRad<scalarD,Cold>(sA,&v,1E-7f,boost::shared_ptr<KrylovMatrix<scalarD> >(new DivFreeKrylovMatrix<scalarD>(sB)));
    boost::filesystem::ofstream os("./stokesEigenMax.txt");
    os << ev << std::endl;
    os << v << std::endl;
  }
  return 0;
}
scalarD evalV(scalar x,scalar y,scalar z)
{
  return cos(x)+y*z+exp(x+z);
}
scalarD evalV(Vec3 pos)
{
  return evalV(pos[0],pos[1],pos[2]);
}
Vec3d evalGrad(scalar x,scalar y,scalar z)
{
  return Vec3d(exp(x+z)-sin(x),z,exp(x+z)+y);
}
Vec3d evalGrad(Vec3 pos)
{
  return evalGrad(pos[0],pos[1],pos[2]);
}
int main(int argc,char** argv)//TestGroundTruth
{
  //read mesh
  if(argc < 2) {
    INFO("Usage: exe [mesh].vtk [boundary_value_file]")
    return -1;
  }
  DHMMesh mesh;
  mesh.readVTK(argv[1]);

  DHMCVField vel(mesh);
  DHMPressureSolver sol(mesh,vel);
  sol.precomputeLHS();
  for(sizeType i=0;i<mesh.nrCell();i++) {
//#define EXACT
#ifdef EXACT
    for(sizeType j=0;j<8;j++)
      vel.getCompRef(i)[j]=evalV(mesh.getCell(i)._verts[j]->_pos);
#else
    Vec3 pos=Vec3::Zero();
    for(sizeType j=0;j<8;j++)
      pos+=mesh.getCell(i)._verts[j]->_pos;
    pos/=8.0f;

    Vec3d grad=evalGrad(pos);
    for(sizeType j=0;j<8;j++)
      vel.getCompRef(i)[j]=grad.dot(mesh.getCell(i)._verts[j]->_pos.cast<scalarD>());
#endif
  }
  sol.project();

  Cold err=Cold::Zero(mesh.nrVert());
  for(sizeType i=0;i<mesh.nrVert();i++)
    err[i]=sol.getPressure().getComp(i)[0]-evalV(mesh.getVert(i)._pos);

  boost::filesystem::ofstream os("./err.txt");
  err.array()-=err.sum()/(scalarD)err.size();
  for(sizeType i=0;i<err.size();i++) {
    os << err[i] << " ";
  }

  DHMField<scalarD,false,1> ERR=sol.getPressure();
  for(sizeType i=0;i<err.size();i++) {
    ERR.getCompRef(i)[0]=err[i];
  }
  ERR.writeSVTK("./err.vtk");

#ifdef EXACT
  sol.getPressure().writeSVTK("./exact.vtk");
#else
  sol.getPressure().writeSVTK("./inexact.vtk");
#endif
  return 0;
}

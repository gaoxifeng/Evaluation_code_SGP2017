#ifndef BLOCKED_SOLVER_H
#define BLOCKED_SOLVER_H

#include "KrylovMatrix.h"
#include "LinearSolver.h"

PRJ_BEGIN

//the krylov matrix
template <int N,typename T,typename KERNEL_TYPE=Kernel<T> >
struct BlockedKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE> {
public:
  typedef typename KERNEL_TYPE::Vec Vec;
  BlockedKrylovMatrix() {}
  BlockedKrylovMatrix(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > inner):_matrix(inner) {}
  BlockedKrylovMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& smat):_matrix(new DefaultKrylovMatrix<T,KERNEL_TYPE>(smat)) {}
  virtual void multiply(const Vec& b,Vec& out) const {
    Vec bBLK,outBLK;
    bBLK.resize(_matrix->n());
    outBLK.resize(_matrix->n());
    for(int d=0; d<N; d++) {
      bBLK=Eigen::Map<const Vec,0,Eigen::InnerStride<N> >(b.data()+d,bBLK.size());
      _matrix->multiply(bBLK,outBLK);
      Eigen::Map<Vec,0,Eigen::InnerStride<N> >(out.data()+d,bBLK.size())=outBLK;
    }
  }
  virtual sizeType n() const {
    return _matrix->n()*N;
  }
  boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > _matrix;
};

//the solver
template <int N,typename SOLVER_PARENT>
struct BlockedSolver : public SOLVER_PARENT {
public:
  typedef typename SOLVER_PARENT::Vec Vec;
  typedef typename Vec::Scalar T;
  typedef Kernel<T> KERNEL_TYPE;
  using SOLVER_PARENT::SOLVER_RESULT;
  using SOLVER_PARENT::_residualOut;
  using SOLVER_PARENT::_iterationsOut;
  virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {
    ASSERT_MSG(false,"Use setKrylovMatrix for BlockedSolver!")
  }
  virtual void setKrylovMatrix(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > matrix) {
    boost::shared_ptr<KrylovMatrix<T> > kry=boost::dynamic_pointer_cast<BlockedKrylovMatrix<N,T,KERNEL_TYPE> >(matrix)->_matrix;
    boost::shared_ptr<DefaultKrylovMatrix<T,KERNEL_TYPE> > defKry=boost::dynamic_pointer_cast<DefaultKrylovMatrix<T,KERNEL_TYPE> >(kry);
    if(defKry)
      _sol.setMatrix(defKry->_fixedMatrix,true);
    else _sol.setKrylovMatrix(kry);
  }
  virtual void setSolverParameters(T toleranceFactor,sizeType maxIterations) {
    _sol.setSolverParameters(toleranceFactor,maxIterations);
  }
  virtual typename SOLVER_PARENT::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
    Vec rhsBLK,resultBLK;
    rhsBLK.resize(rhs.size()/N);
    resultBLK.resize(rhs.size()/N);
    for(int d=0; d<N; d++) {
      rhsBLK=Eigen::Map<const Vec,0,Eigen::InnerStride<N> >(rhs.data()+d,rhsBLK.size());
      typename SOLVER_PARENT::SOLVER_RESULT ret=_sol.solve(rhsBLK,resultBLK);
      _residualOut=max<T>(_residualOut,_sol.getResidual());
      _iterationsOut=max<sizeType>(_iterationsOut,_sol.getIterationsCount());
      if(ret != SOLVER_PARENT::SUCCESSFUL)
        return ret;
      Eigen::Map<Vec,0,Eigen::InnerStride<N> >(result.data()+d,resultBLK.size())=resultBLK;
    }
    return SOLVER_PARENT::SUCCESSFUL;
  }
  virtual sizeType getIterationsCount() const {
    return _iterationsOut;
  }
  virtual T getResidual() const {
    return _residualOut;
  }
  virtual void setCallback(typename boost::shared_ptr<Callback<T,KERNEL_TYPE> > cb) {
    _sol.setCallback(cb);
  }
  SOLVER_PARENT _sol;
};

PRJ_END

#endif

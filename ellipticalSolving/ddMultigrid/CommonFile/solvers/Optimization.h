#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "Objective.h"
#include "eigen3/Eigen/SparseQR"
#include <boost/property_tree/ptree.hpp>

PRJ_BEGIN

//The optimization framework tries to optimize:
//\E{argmin}_{X}&&\Sigma_i \|Term_i(X)\|^2
//\E{s.t.}&&\foreach j C_j(X)>=/==0

//Use class Term for Term_i(X)
//Use class GEQConstraint for C_j(X)>=0
//Use class EQConstraint for C_j(X)==0

//Message/Parameter passing is through boost::property_tree
//To choose optimizer, use: pt.put<string>("optimizer",TYPE)
//Here we support several types of optimizer:
//use TYPE=LBFGS to use buildin CommonFile LBFGS solver (just a port of commerical software)
//use TYPE=BLEIC to use alglib BLEIC solver (support linear/box constraints)
//use TYPE=AUGLAG to use augmented-lagrangian solver with alglib BLEIC as subproblem solver (support all constraints)
//use TYPE=LM to use alglib levenberg-marquardt solver (support box constraints)
//use TYPE=ELM to use eigen levenberg-marquardt solver (no constraints supported)
//Checkout Optimizer.cpp to find more tunable parameters
//Checkout Optimizer::runExamples() to find usage

class Optimizer;
class Term
{
public:
  typedef scalarD T;
  typedef Objective<scalarD>::Vec Vec;
  //internal functions
  Term(const Optimizer& opt);
  virtual Vec val(const Vec& x);
  virtual void grad(const Vec& x,Matd& grad);
  virtual void addGrad(const Vec& x,Vec& grad,Vec alpha=Vec());
  virtual Vec valGrad(const Vec& x,Matd* grad=NULL);
  void setReuse(bool reuse);
  //user implemented functions
  virtual sizeType values() const=0;
  //default to true
  virtual bool isLinear() const;
  //for constraint, default to exit(-1)
  virtual char sign() const;
protected:
  //just implement this: return the value of Term_i(X) or C_j(X)
  //and set grad=grad+\FPP{Term_i/C_j}{X}*alpha (not needed for gradient free optimization)
  virtual Vec valGradInner(const Vec& x,Matd& grad)=0;
  const boost::property_tree::ptree& _pt;
  bool _reuse;
  Vec _valCache;
private:
  Matd _gradCache;
};
class SparseTerm : public Term
{
public:
  using Term::T;
  using Term::Vec;
  typedef vector<Eigen::Triplet<scalarD,sizeType> > SMatd;
  SparseTerm(const Optimizer& opt);
  virtual void addGrad(const Vec& x,Vec& grad,Vec alpha=Vec());
  virtual Vec valGrad(const Vec& x,Matd* grad=NULL);
  virtual Vec svalGrad(const Vec& x,SMatd* grad=NULL);
protected:
  virtual Vec valGradInner(const Vec& x,Matd& grad);
  virtual Vec valGradInner(const Vec& x,SMatd& grad)=0;
private:
  SMatd _sgradCache;
};
class GEQConstraint : public Term
{
public:
  GEQConstraint(const Optimizer& opt):Term(opt) {}
  virtual char sign() const;
};
class EQConstraint : public Term
{
public:
  EQConstraint(const Optimizer& opt):Term(opt) {}
  virtual char sign() const;
};
class SparseGEQConstraint : public SparseTerm
{
public:
  SparseGEQConstraint(const Optimizer& opt):SparseTerm(opt) {}
  virtual char sign() const;
};
class SparseEQConstraint : public SparseTerm
{
public:
  SparseEQConstraint(const Optimizer& opt):SparseTerm(opt) {}
  virtual char sign() const;
};
class Optimizer : public Objective<scalarD>
{
public:
  typedef scalarD T;
  typedef scalarD Scalar;
  typedef sizeType Index;
  typedef Objective<scalarD>::Vec InputType;
  typedef Objective<scalarD>::Vec ValueType;
  typedef Eigen::SparseMatrix<Scalar,Eigen::ColMajor,int> JacobianType;
  typedef Eigen::SparseQR<JacobianType,Eigen::COLAMDOrdering<int> > QRSolver;
  //methods
  Optimizer();
  Optimizer(const boost::property_tree::ptree& pt);
  scalarD minimize(Vec &x);
  void addTerm(boost::shared_ptr<Term> term);
  void addConstraint(boost::shared_ptr<Term> term);
  const boost::property_tree::ptree& getProps() const;
  boost::property_tree::ptree& getProps();
  //Overwritten of Objective<scalarD>
  //For LBFGS Optimization
  virtual int operator()(const Vec& x,T& FX,Vec& DFDX,const T& step,bool wantGradient);
  //For LM Optimization
  virtual int operator()(const Vec& x,Vec& fvec);
  virtual int df(const Vec& x,Eigen::MatrixXd& fjac);
  virtual int df(const Vec& x,JacobianType& jacobian);
  //For Finite Difference Debug
  virtual void debugFD(scalarD delta=1E-8f,bool debugJac=false,bool debugGrad=false);
  //Dimension Info
  virtual int inputs() const;
  virtual int values() const;
  virtual int constraints() const;
  virtual void profileLineSearch(const sizeType& k,const Vec& x,const Vec& d,const T& step);
private:
  void extractConstraint(Vec& lb,Vec& ub,Matd& c,Coli& ct,vector<boost::shared_ptr<Term> >& nl) const;
  boost::property_tree::ptree _pt;
  vector<boost::shared_ptr<Term> > _terms,_cons;
  Vec _lastX;
};
class OptimizerTest
{
public:
  typedef Optimizer::Vec Vec;
  typedef Optimizer::T T;
  //example1: 1D diffusion with constraint
  static Vec runExample1(sizeType N,bool cons,bool consTotal,bool consNorm,const string& type);
  static void runExample1();
  //example2: 2D cloth with constraint
  static void runExample2(sizeType N,bool cons,bool consTwoPlane,bool consBall,const string& type);
  static void runExample2();
};

PRJ_END

#endif

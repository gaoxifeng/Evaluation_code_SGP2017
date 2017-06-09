#include "Optimization.h"
//LBFGS
#include "Minimizer.h"
//EM Levenberg-Marquardt
#include "eigen3/Eigen/IterativeLinearSolvers"
//alglib
#include "alglib/optimization.h"

USE_PRJ_NAMESPACE

//Term interface
Term::Term(const Optimizer& opt):_pt(opt.getProps()) {}
Term::Vec Term::val(const Vec& x)
{
  return valGrad(x);
}
void Term::grad(const Vec& x,Matd& grad)
{
  valGrad(x,&grad);
}
void Term::addGrad(const Vec& x,Vec& grad,Vec alpha)
{
  valGrad(x,NULL);
  if(alpha.size() == 0)
    alpha.setOnes(values());
  grad+=_gradCache.transpose()*alpha;
}
Term::Vec Term::valGrad(const Vec& x,Matd* grad)
{
  if(!_reuse) {
    _gradCache.setZero(values(),x.size());
    _valCache=valGradInner(x,_gradCache);
    ASSERT(values() == _valCache.size())
  }
  if(grad)
    *grad=_gradCache;
  return _valCache;
}
void Term::setReuse(bool reuse)
{
  _reuse=reuse;
}
bool Term::isLinear() const
{
  return true;
}
char Term::sign() const
{
  ASSERT_MSG(false,"Inherit (Sparse)(G)EQConstraint for constraints!")
  return 0;
}

//SparseTerm interface
SparseTerm::SparseTerm(const Optimizer& opt):Term(opt) {}
void SparseTerm::addGrad(const Vec& x,Vec& grad,Vec alpha)
{
  valGrad(x,NULL);
  if(alpha.size() == 0)
    alpha.setOnes(values());
  for(sizeType i=0; i<(sizeType)_sgradCache.size(); i++)
    grad[_sgradCache[i].col()]+=_sgradCache[i].value()*alpha[_sgradCache[i].row()];
}
Term::Vec SparseTerm::valGrad(const Vec& x,Matd* grad)
{
  Vec ret=svalGrad(x);
  if(grad) {
    grad->setZero(values(),x.size());
    for(sizeType i=0; i<(sizeType)_sgradCache.size(); i++)
      (*grad)(_sgradCache[i].row(),_sgradCache[i].col())+=_sgradCache[i].value();
  }
  return ret;
}
Term::Vec SparseTerm::svalGrad(const Vec& x,SMatd* grad)
{
  if(!_reuse) {
    _sgradCache.clear();
    _valCache=valGradInner(x,_sgradCache);
    ASSERT(values() == _valCache.size())
  }
  if(grad)
    *grad=_sgradCache;
  return _valCache;
}
Term::Vec SparseTerm::valGradInner(const Vec& x,Matd& grad)
{
  ASSERT_MSG(false,"Use the sparse version for SparseTerm!")
  return Vec();
}

//constraint interface
char EQConstraint::sign() const
{
  return 0;
}
char GEQConstraint::sign() const
{
  return 1;
}
char SparseEQConstraint::sign() const
{
  return 0;
}
char SparseGEQConstraint::sign() const
{
  return 1;
}

//term utility
sizeType COUNT_TERM(const vector<boost::shared_ptr<Term> >& terms)
{
  sizeType nrTerm=0;
  for(sizeType i=0; i<(sizeType)terms.size(); i++)
    nrTerm+=terms[i]->values();
  return nrTerm;
}
#define LOOP_TERM(v) for(sizeType i=0,off=0;i<(sizeType)v.size();off+=v[i]->values(),i++)
bool IS_BOX_CONSTRAINT(Term::Vec& lb,Term::Vec& ub,SparseTerm& term)
{
  typedef Eigen::SparseMatrix<scalarD,Eigen::RowMajor,sizeType> SMat;
  SparseTerm::SMatd grad;
  Term::Vec val=term.svalGrad(Term::Vec::Zero(lb.size()),&grad);

  //the sparse matrix
  SMat gradM;
  gradM.resize(val.size(),lb.size());
  gradM.setFromTriplets(grad.begin(),grad.end());

  //detect
  Term::Vec coef=Term::Vec::Zero(val.size());
  Coli id=Coli::Constant(val.size(),-1);
  for(sizeType r=0; r<val.size(); r++) {
    for(SMat::InnerIterator it(gradM,r); it; ++it) {
      //simply ignore zero terms
      //user should call prune explicitly
      if(it.value() == 0)
        continue;
      else if(id[r] == -1) {
        coef[r]=it.value();
        id[r]=it.col();
      } else {
        id[r]=-1;
        break;
      }
    }
    if(id[r] == -1)
      return false;
  }

  //setup
  val.array()/=-coef.array();
  for(sizeType r=0; r<val.size(); r++) {
    if(term.sign() == 0) {
      lb[id[r]]=max(lb[id[r]],val[r]);
      ub[id[r]]=min(ub[id[r]],val[r]);
    } else if(coef[r] > 0)
      lb[id[r]]=max(lb[id[r]],val[r]);
    else ub[id[r]]=min(ub[id[r]],val[r]);
  }
  return true;
}

//alglib interface
Cold fromAlglib(const alglib::real_1d_array& xAlglib)
{
  Eigen::Map<const Cold> xMap(xAlglib.getcontent(),xAlglib.length());
  return xMap;
}
void toAlglib(alglib::real_1d_array& ret,const Cold& x)
{
  if(ret.length() != x.size())
    ret.setlength(x.size());
  memcpy(ret.getcontent(),x.data(),x.size()*sizeof(Cold::Scalar));
}
alglib::real_1d_array toAlglib(const Cold& x)
{
  alglib::real_1d_array ret;
  toAlglib(ret,x);
  return ret;
}
void toAlglib(alglib::real_2d_array& ret,const Matd& x)
{
  if(ret.rows() != x.rows() || ret.cols() != x.cols())
    ret.setlength(x.rows(),x.cols());
  Matd xt=x.transpose();
  for(sizeType i=0; i<ret.rows(); i++)
    memcpy(ret[i],xt.data()+i*ret.cols(),ret.cols()*sizeof(Matd::Scalar));
}
alglib::real_2d_array toAlglib(const Matd& x)
{
  alglib::real_2d_array ret;
  toAlglib(ret,x);
  return ret;
}
alglib::integer_1d_array toAlglib(const Coli& x)
{
  alglib::integer_1d_array ret;
  ret.setcontent(x.size(),x.cast<alglib::ae_int_t>().eval().data());
  return ret;
}
void gradAlglib(const alglib::real_1d_array& xAlglib,scalarD& val,alglib::real_1d_array& grad,void* ptr)
{
  Optimizer* opt=(Optimizer*)ptr;
  Cold x=fromAlglib(xAlglib),gradEigen;
  (*opt)(x,val,gradEigen,1,true);
  toAlglib(grad,gradEigen);
}
void fvecAlglib(const alglib::real_1d_array& xAlglib,alglib::real_1d_array& fvec,void* ptr)
{
  Optimizer* opt=(Optimizer*)ptr;
  Cold x=fromAlglib(xAlglib);
  Cold fvecEigen;
  (*opt)(x,fvecEigen);
  toAlglib(fvec,fvecEigen);
}
void fjacAlglib(const alglib::real_1d_array& xAlglib,alglib::real_1d_array& fvec,alglib::real_2d_array& fjac,void* ptr)
{
  Optimizer* opt=(Optimizer*)ptr;
  Cold x=fromAlglib(xAlglib);
  Cold fvecEigen;
  Matd fjacEigen;
  (*opt)(x,fvecEigen);
  (*opt).df(x,fjacEigen);
  toAlglib(fvec,fvecEigen);
  toAlglib(fjac,fjacEigen);
}

//auglag interface
class AugLag
{
public:
  typedef Optimizer::Vec Vec;
  AugLag(vector<boost::shared_ptr<Term> >& terms,vector<boost::shared_ptr<Term> >& nl,Optimizer& inner)
    :_terms(terms),_nl(nl),_inner(inner) {}
  void minimize(alglib::minbleicstate& state,alglib::real_1d_array& x,scalarD epsg,scalarD epsf,scalarD epsx,sizeType nrIter) {
    _lambda=Vec::Zero(COUNT_TERM(_nl));
    _mu=1E3f;

    sizeType iter=0;
    for(; iter<nrIter;) {
      //inner loop
      alglib::minbleicsetcond(state,epsg,epsf,epsx,nrIter-iter);
      if(iter > 0)
        alglib::minbleicrestartfrom(state,x);
      alglib::minbleicoptimize(state,valGrad,NULL,this);

      //update lambda
      alglib::minbleicreport rep;
      alglib::minbleicresults(state,x,rep);
      const Cold xEigen=fromAlglib(x);
      LOOP_TERM(_nl)
      _lambda.segment(off,_nl[i]->values())=clampCons(_lambda.segment(off,_nl[i]->values())+_mu*_nl[i]->val(xEigen),_nl[i]->sign());

      //count number of iterations consumed
      if(rep.iterationscount == 0)break;
      iter+=rep.iterationscount;

      //fallback to BLEIC
      if(_nl.empty())
        break;
    }
  }
private:
  static void valGrad(const alglib::real_1d_array& xAlglib,scalarD& val,alglib::real_1d_array& grad,void* ptr) {
    //init
    AugLag* opt=(AugLag*)(ptr);
    const Vec& lambda=opt->_lambda;
    const scalarD mu=opt->_mu;

    //toEigen
    const Cold xEigen=fromAlglib(xAlglib);
    Vec gradEigen=Vec::Zero(grad.length());

    //compute gradient
    opt->_inner(xEigen,val,gradEigen,1,true);
    LOOP_TERM(opt->_nl) {
      Term& nl=*(opt->_nl[i]);
      nl.setReuse(false);
      Vec cValClamp=clampCons(nl.val(xEigen)+lambda.segment(off,opt->_nl[i]->values())/mu,nl.sign());
      val+=mu*0.5f*cValClamp.squaredNorm();
      nl.setReuse(true);
      nl.addGrad(xEigen,gradEigen,mu*cValClamp);
    }
    toAlglib(grad,gradEigen);
  }
  static Vec clampCons(Vec x,char sign) {
    return (sign<0) ? x.cwiseMax(Vec::Zero(x.size())) :
           (sign>0) ? x.cwiseMin(Vec::Zero(x.size())) : x;
  }
  vector<boost::shared_ptr<Term> >& _terms;
  vector<boost::shared_ptr<Term> >& _nl;
  Optimizer& _inner;
  Vec _lambda;
  scalarD _mu;
};

//Optimizer interface
Optimizer::Optimizer() {}
Optimizer::Optimizer(const boost::property_tree::ptree& pt):_pt(pt) {}
scalarD Optimizer::minimize(Vec& x)
{
  scalarD gtol=_pt.get<scalarD>("gtol",0);
  scalarD ftol=_pt.get<scalarD>("ftol",sqrt(numeric_limits<scalarD>::epsilon()));
  scalarD xtol=_pt.get<scalarD>("xtol",sqrt(numeric_limits<scalarD>::epsilon()));
  int nrIter=_pt.get<int>("nrIter",numeric_limits<int>::max());

  //constraints
  Vec lb,ub;
  Matd c;
  Coli ct;

  string type=_pt.get<string>("type");
  if(type == "LBFGS") {
    ASSERT_MSG(_cons.empty(),"LBFGS doesn't support any constraints!")
    scalarD fx=0;
    NoCallback<scalarD,Kernel<scalarD> > cb;
    LBFGSMinimizer<scalarD> sol;
    sol.maxIterations()=nrIter;
    sol.nrCorrect()=_pt.get<int>("nrCorrect",10);
    sol.epsilon()=_pt.get<scalarD>("gtol_over_xtol",1E-5f);
    sol.minimize(x,fx,*this,cb);
  } else if(type == "BLEIC" || type == "AUGLAG") {
    alglib::real_1d_array xAlglib=toAlglib(x);
    alglib::minbleicstate state;
    alglib::minbleicreport rep;
    alglib::minbleiccreate(xAlglib,state);
    vector<boost::shared_ptr<Term> > nl;
    {
      extractConstraint(lb,ub,c,ct,nl);
      alglib::minbleicsetbc(state,toAlglib(lb),toAlglib(ub));
      if(ct.size() > 0)
        alglib::minbleicsetlc(state,toAlglib(c),toAlglib(ct));
      if(type != "AUGLAG")
        ASSERT_MSG(nl.empty(),"alglib BLEIC doesn't support nonlinear constraints!")
      }
    AugLag opt(_terms,nl,*this);
    opt.minimize(state,xAlglib,gtol,ftol,xtol,nrIter);
    alglib::minbleicresults(state,xAlglib,rep);
    x=fromAlglib(xAlglib);
  } else if(type == "LM") {
    alglib::real_1d_array xAlglib=toAlglib(x);
    alglib::minlmstate state;
    alglib::minlmreport rep;
    alglib::minlmcreatevj(values(),xAlglib,state);
    {
      vector<boost::shared_ptr<Term> > nl;
      extractConstraint(lb,ub,c,ct,nl);
      alglib::minlmsetbc(state,toAlglib(lb),toAlglib(ub));
      ASSERT_MSG(ct.size() == 0,"alglib LM doesn't support general linear constraints!")
      ASSERT_MSG(nl.empty(),"alglib LM doesn't support nonlinear constraints!")
    }
    alglib::minlmsetcond(state,gtol,ftol,xtol,nrIter);
    alglib::minlmoptimize(state,fvecAlglib,fjacAlglib,NULL,this);
    alglib::minlmresults(state,xAlglib,rep);
    x=fromAlglib(xAlglib);
  } else if(type == "ELM") {
    ASSERT_MSG(_cons.empty(),"Eigen LM doesn't support any constraints!")
    //Eigen::LevenbergMarquardt<Optimizer> lm(*this);
    //lm.setMaxfev(nrIter);
    //lm.setFtol(ftol);
    //lm.setXtol(xtol);
    //lm.setGtol(gtol);
    //lm.minimize(x);
    Vec fvec,RHS,dx;
    JacobianType fjac;
    operator()(x,fvec);
    scalarD F=fvec.squaredNorm();
    scalarD offDiag=200.0f;
    for(int it=0; it<nrIter; it++) {
      //update
      operator()(x,fvec);
      df(x,fjac);
      JacobianType LHS=fjac.transpose()*fjac;
      OMP_PARALLEL_FOR_
      for(sizeType i=0; i<LHS.rows(); i++)
        LHS.coeffRef(i,i)+=offDiag;
      RHS=fjac.transpose()*fvec;
      if(RHS.norm() < gtol)
        break;
      dx=-Eigen::ConjugateGradient<JacobianType>(LHS).solve(RHS);
      //dx=-Eigen::SimplicialCholesky<JacobianType>(LHS).solve(RHS);
      //INFOV("%f",(LHS*dx+RHS).norm())
      if(dx.norm() < xtol)
        break;
      x+=dx;

      //check
      operator()(x,fvec);
      scalarD FNew=fvec.squaredNorm(),pred=2*RHS.dot(dx);
      INFOV("SparseLM Achieved F=%f!",FNew)
      if(abs(F-FNew) < ftol)
        break;
      else if(F-FNew > 0.75f*pred) {
        offDiag=max(offDiag*0.75f,(scalarD)0.001f);
      } else if(F-FNew < 0.25f*pred) {
        offDiag/=0.75f;
        if(F-FNew <= 0) {
          x-=dx;
          continue;
        }
      }
      F=FNew;
    }
  } else ASSERT_MSG(false,"Unknown Optimizer Type!")

    Vec dfdx=x;
  scalarD fx;
  operator()(x,fx,dfdx,1,false);
  return fx;
}
void Optimizer::addTerm(boost::shared_ptr<Term> term)
{
  _terms.push_back(term);
}
void Optimizer::addConstraint(boost::shared_ptr<Term> term)
{
  _cons.push_back(term);
}
const boost::property_tree::ptree& Optimizer::getProps() const
{
  return _pt;
}
boost::property_tree::ptree& Optimizer::getProps()
{
  return _pt;
}
//Overwritten of Objective<scalarD>
int Optimizer::operator()(const Vec& x,T& FX,Vec& DFDX,const T& step,bool wantGradient)
{
  bool reuse=_lastX.size() == x.size() && _lastX == x;
  FX=0;
  DFDX.setZero(inputs());
  LOOP_TERM(_terms) {
    _terms[i]->setReuse(reuse);
    Vec val=_terms[i]->val(x);
    FX+=val.squaredNorm();

    if(wantGradient) {
      _terms[i]->setReuse(true);
      _terms[i]->addGrad(x,DFDX,2*val);
    }
  }
  //INFOV("F=%f",FX)
  _lastX=x;
  return 0;
}
//For LM Optimization
int Optimizer::operator()(const Vec& x,Vec& fvec)
{
  bool reuse=_lastX.size() == x.size() && _lastX == x;
  fvec.setZero(values());
  LOOP_TERM(_terms) {
    _terms[i]->setReuse(reuse);
    fvec.segment(off,_terms[i]->values())=_terms[i]->val(x);
  }
  _lastX=x;
  return 0;
}
int Optimizer::df(const Vec& x,Eigen::MatrixXd& fjac)
{
  //JacobianType jac;
  //df(x,jac);
  //fjac=jac.toDense();
  //return 0;

  bool reuse=_lastX.size() == x.size() && _lastX == x;
  fjac.setZero(values(),inputs());
  LOOP_TERM(_terms) {
    Matd df=Matd::Zero(_terms[i]->values(),inputs());
    _terms[i]->setReuse(reuse);
    _terms[i]->grad(x,df);
    fjac.block(off,0,df.rows(),df.cols())=df;
  }
  _lastX=x;
  return 0;
}
int Optimizer::df(const Vec& x,JacobianType& fjac)
{
  bool reuse=_lastX.size() == x.size() && _lastX == x;
  SparseTerm::SMatd trips,currTrips;
  LOOP_TERM(_terms) {
    _terms[i]->setReuse(reuse);
    if(boost::dynamic_pointer_cast<SparseTerm>(_terms[i])) {
      boost::dynamic_pointer_cast<SparseTerm>(_terms[i])->svalGrad(x,&currTrips);

      sizeType offT=(sizeType)trips.size();
      trips.resize(trips.size()+currTrips.size());
      OMP_PARALLEL_FOR_
      for(sizeType i=0; i<(sizeType)currTrips.size(); i++)
        trips[i+offT]=Eigen::Triplet<scalarD,sizeType>(currTrips[i].row()+off,currTrips[i].col(),currTrips[i].value());
    } else {
      Matd df=Matd::Zero(_terms[i]->values(),inputs());
      _terms[i]->grad(x,df);

      for(sizeType r=0; r<df.rows(); r++)
        for(sizeType c=0; c<df.cols(); c++)
          trips.push_back(Eigen::Triplet<scalarD,sizeType>(r+off,c,df(r,c)));
    }
  }
  fjac.resize(values(),inputs());
  fjac.setFromTriplets(trips.begin(),trips.end());
  _lastX=x;
  return 0;
}
//For Finite Difference Debug
void Optimizer::debugFD(scalarD delta,bool debugJac,bool debugGrad)
{
  Matd fjac;
  Vec x=Vec::Random(inputs()),fvec,fvec2;

  //add all constraints to terms
  vector<boost::shared_ptr<Term> > terms=_terms;
  _terms.insert(_terms.end(),_cons.begin(),_cons.end());

  //pass 1
  scalarD maxErr=0;
  if(debugJac) {
    operator()(x,fvec);
    df(x,fjac);
    for(sizeType i=0; i<inputs(); i++) {
      Vec tmp=x;
      tmp[i]+=delta;
      operator()(tmp,fvec2);

      scalarD jacAbs=fjac.col(i).norm();
      scalarD jacErr=(fjac.col(i)-(fvec2-fvec)/delta).norm();
      if(_pt.get<bool>("printFDError",false))
        INFOV("fjac-Norm: %f Err: %f!",jacAbs,jacErr);
      if(jacAbs > ScalarUtil<scalarD>::scalar_eps)
        maxErr=max(maxErr,jacErr/jacAbs);
    }
    INFOV("Max Jacobian Error: %f!",maxErr)
  }

  //pass 2
  maxErr=0;
  if(debugGrad) {
    scalarD fx,fx2;
    Vec grad,grad2;
    operator()(x,fx,grad,1,true);
    for(sizeType i=0; i<inputs(); i++) {
      Vec tmp=x;
      tmp[i]+=delta;
      operator()(tmp,fx2,grad2,1,false);

      scalarD gradErr=(grad[i]-(fx2-fx)/delta);
      if(_pt.get<bool>("printFDError",false))
        INFOV("fjac-Norm: %f Err: %f!",grad[i],gradErr);
      if(abs(grad[i]) > ScalarUtil<scalarD>::scalar_eps)
        maxErr=max(maxErr,gradErr/grad[i]);
    }
    INFOV("Max Grad Error: %f!",maxErr)
  }

  //restore
  _terms=terms;
}
//Dimension Info
int Optimizer::inputs() const
{
  return _pt.get<int>("inputs",-1);
}
int Optimizer::values() const
{
  return COUNT_TERM(_terms);
}
int Optimizer::constraints() const
{
  return COUNT_TERM(_cons);
}
void Optimizer::profileLineSearch(const sizeType& k,const Vec& x,const Vec& d,const T& step)
{
  if(!_pt.get<bool>("debug",false))
    return;
  INFOV("Line Search Step Size=%f!",step)
}
void Optimizer::extractConstraint(Vec& lb,Vec& ub,Matd& c,Coli& ct,vector<boost::shared_ptr<Term> >& nl) const
{
  Vec x=Vec::Zero(inputs());
  lb.setConstant(inputs(),-ScalarUtil<scalarD>::scalar_max);
  ub.setConstant(inputs(),ScalarUtil<scalarD>::scalar_max);
  c.setZero(0,0);
  ct.setZero(0);

  //pass 1: build box constraints/collect number of general linear constraints
  sizeType nrBox=0;
  vector<boost::shared_ptr<Term> > gc;
  for(sizeType i=0; i<(sizeType)_cons.size(); i++) {
    if(!_cons[i]->isLinear()) {
      nl.push_back(_cons[i]);
      continue;
    }
    //you must use SparseTerm for box constraint
    boost::shared_ptr<SparseTerm> term=boost::dynamic_pointer_cast<SparseTerm>(_cons[i]);
    if(!term) {
      gc.push_back(_cons[i]);
      continue;
    }
    if(!IS_BOX_CONSTRAINT(lb,ub,*term))
      gc.push_back(term);
    else nrBox+=term->values();
  }

  //pass 2: build general linear constraints
  ct.resize(COUNT_TERM(gc));
  c.resize(COUNT_TERM(gc),x.size()+1);
  LOOP_TERM(gc) {
    Matd grad;
    c.block(off,inputs(),gc[i]->values(),1)=-gc[i]->valGrad(x,&grad);
    c.block(off,0,gc[i]->values(),inputs())=grad;
    ct.segment(off,gc[i]->values()).setConstant(gc[i]->sign());
  }

  //profile
  if(!_pt.get<bool>("quiet",false))
    INFOV("Found %ld box constraints, %ld linear constraints, %ld nonlinear constraints!",nrBox,c.rows(),COUNT_TERM(nl))
  }

//example1: 1D Diffusion with constraint
class DirichletTerm : public SparseTerm
{
public:
  DirichletTerm(const Optimizer& opt,sizeType i,sizeType j):SparseTerm(opt),_i(i),_j(j) {}
  virtual sizeType values() const {
    return 1;
  }
protected:
  virtual Vec valGradInner(const Vec& x,SMatd& grad) {
    grad.clear();
    grad.push_back(Eigen::Triplet<scalarD,sizeType>(0,_i,1));
    grad.push_back(Eigen::Triplet<scalarD,sizeType>(0,_j,-1));

    Vec ret=Vec::Zero(1);
    ret[0]=x[_i]-x[_j];
    return ret;
  }
  sizeType _i,_j;
};
class FixTerm : public SparseEQConstraint
{
public:
  FixTerm(const Optimizer& opt,sizeType i,scalarD goal,scalarD strength)
    :SparseEQConstraint(opt),_i(i),_goal(goal),_strength(strength) {}
  virtual sizeType values() const {
    return 1;
  }
protected:
  virtual Vec valGradInner(const Vec& x,SMatd& grad) {
    grad.clear();
    grad.push_back(Eigen::Triplet<scalarD,sizeType>(0,_i,_strength));

    Vec ret=Vec::Zero(1);
    ret[0]=(x[_i]-_goal)*_strength;
    return ret;
  }
  sizeType _i;
  scalarD _goal,_strength;
};
class TotalTerm : public EQConstraint
{
public:
  TotalTerm(const Optimizer& opt,scalarD goal):EQConstraint(opt),_goal(goal) {}
  virtual sizeType values() const {
    return 1;
  }
protected:
  virtual Vec valGradInner(const Vec& x,Matd& grad) {
    grad.setConstant(1,x.size(),-1);
    Vec ret=Vec::Zero(1);
    ret[0]=_goal-x.sum();
    return ret;
  }
  scalarD _goal;
};
class NormTerm : public GEQConstraint
{
public:
  NormTerm(const Optimizer& opt,scalarD goal):GEQConstraint(opt),_goal(goal) {}
  virtual sizeType values() const {
    return 1;
  }
  bool isLinear() const {
    return false;
  }
protected:
  virtual Vec valGradInner(const Vec& x,Matd& grad) {
    grad.resize(1,x.size());
    grad.row(0)=-2*x;
    Vec ret=Vec::Zero(1);
    ret[0]=_goal-x.squaredNorm();
    return ret;
  }
  scalarD _goal;
};
OptimizerTest::Vec OptimizerTest::runExample1(sizeType N,bool cons,bool consTotal,bool consNorm,const string& type)
{
  Optimizer opt;
  opt.getProps().put<int>("inputs",N);
  opt.getProps().put<string>("type",type);
  for(sizeType i=1; i<N; i++)
    opt.addTerm(boost::shared_ptr<DirichletTerm>(new DirichletTerm(opt,i-1,i)));
  if(cons) { //add soft constraint
    opt.addConstraint(boost::shared_ptr<FixTerm>(new FixTerm(opt,0,0,10.0f)));
    opt.addConstraint(boost::shared_ptr<FixTerm>(new FixTerm(opt,N-1,1,10.0f)));
  } else { //add box hard constraint
    opt.addTerm(boost::shared_ptr<FixTerm>(new FixTerm(opt,0,0,10.0f)));
    opt.addTerm(boost::shared_ptr<FixTerm>(new FixTerm(opt,N-1,1,10.0f)));
  }
  if(consTotal) //add dense hard constraint
    opt.addConstraint(boost::shared_ptr<TotalTerm>(new TotalTerm(opt,1)));
  if(consNorm)  //add nonlinear constraint
    opt.addConstraint(boost::shared_ptr<NormTerm>(new NormTerm(opt,1.5)));
  Vec x=Vec::Random(N);
  scalarD fx=opt.minimize(x);
  opt.debugFD();
  INFOV("Energy=%f, Total=%f, Norm=%f, using %s!",fx,x.sum(),x.squaredNorm(),type.c_str())
  return x;
}
void OptimizerTest::runExample1()
{
  INFO("Unconstrained")
  runExample1(128,false,false,false,"LBFGS");
  runExample1(128,false,false,false,"BLEIC");
  runExample1(128,false,false,false,"AUGLAG");
  runExample1(128,false,false,false,"LM");
  runExample1(128,false,false,false,"ELM");

  INFO("Box Constraint")
  runExample1(128,true,false,false,"BLEIC");
  runExample1(128,true,false,false,"AUGLAG");
  runExample1(128,true,false,false,"LM");

  INFO("Linear Constraint")
  runExample1(128,true,true,false,"BLEIC");
  runExample1(128,true,true,false,"AUGLAG");

  INFO("Nonlinear Constraint")
  runExample1(128,true,false,true,"AUGLAG");
}
//example2: 2D cloth with constraint
#define ID(X,Y) (Y)*_N+(X)
#define IDV(X,Y) ((Y)*_N+(X))*3
template <typename MAT>
void addBLK(SparseTerm::SMatd& smat,sizeType r,sizeType c,const MAT& m)
{
  for(sizeType R=0; R<m.rows(); R++)
    for(sizeType C=0; C<m.cols(); C++)
      smat.push_back(Eigen::Triplet<scalarD,sizeType>(r+R,c+C,m(R,C)));
}
class SpringTerm : public SparseTerm
{
public:
  SpringTerm(const Optimizer& opt,sizeType N,scalarD K):SparseTerm(opt),_N(N),_K(K) {}
  virtual sizeType values() const {
    return _N*(_N-1)*2;
  }
protected:
  virtual Vec valGradInner(const Vec& x,SMatd& grad) {
    grad.clear();
    Vec ret=Vec::Zero(values());
    sizeType off=0;
    for(sizeType i=0; i<_N; i++)
      for(sizeType j=0; j<_N; j++) {
        if(i>0) {
          ret[off]=(x.segment<3>(IDV(i,j))-x.segment<3>(IDV(i-1,j))).squaredNorm()-1;
          addBLK(grad,off,IDV(i,j),2*(x.segment<3>(IDV(i,j))-x.segment<3>(IDV(i-1,j))).transpose()*_K);
          addBLK(grad,off,IDV(i-1,j),2*(x.segment<3>(IDV(i-1,j))-x.segment<3>(IDV(i,j))).transpose()*_K);
          off++;
        }
        if(j>0) {
          ret[off]=(x.segment<3>(IDV(i,j))-x.segment<3>(IDV(i,j-1))).squaredNorm()-1;
          addBLK(grad,off,IDV(i,j),2*(x.segment<3>(IDV(i,j))-x.segment<3>(IDV(i,j-1))).transpose()*_K);
          addBLK(grad,off,IDV(i,j-1),2*(x.segment<3>(IDV(i,j-1))-x.segment<3>(IDV(i,j))).transpose()*_K);
          off++;
        }
      }
    return ret*_K;
  }
  sizeType _N;
  scalarD _K;
};
class GravityTerm : public SparseTerm
{
public:
  GravityTerm(const Optimizer& opt,sizeType N,scalarD g):SparseTerm(opt),_N(N),_g(g) {}
  virtual sizeType values() const {
    return _N*_N;
  }
protected:
  virtual Vec valGradInner(const Vec& x,SMatd& grad) {
    grad.clear();
    scalarD BIG_NUMBER=1E3f;
    Vec ret=Vec::Zero(values());
    for(sizeType i=0; i<_N; i++)
      for(sizeType j=0; j<_N; j++) {
        ret[ID(i,j)]=sqrt(BIG_NUMBER-x.segment<3>(IDV(i,j))[2]*_g);
        grad.push_back(Eigen::Triplet<scalarD,sizeType>(ID(i,j),IDV(i,j)+2,-_g/ret[ID(i,j)]/2));
      }
    return ret;
  }
  sizeType _N;
  scalarD _g;
};
class FixPointTerm : public SparseEQConstraint
{
public:
  FixPointTerm(const Optimizer& opt,sizeType N,scalarD strength)
    :SparseEQConstraint(opt),_N(N),_strength(strength) {}
  virtual sizeType values() const {
    return 6;
  }
protected:
  virtual Vec valGradInner(const Vec& x,SMatd& grad) {
    grad.clear();
    addBLK(grad,0,IDV(0,0),Mat3::Identity()*_strength);
    addBLK(grad,3,IDV(_N-1,0),Mat3d::Identity()*_strength);

    Vec ret=Vec::Zero(6);
    ret.segment<3>(0)=(x.segment<3>(IDV(0,0))-Vec3d(0,0,0))*_strength;
    ret.segment<3>(3)=(x.segment<3>(IDV(_N-1,0))-Vec3d(_N-1,0,0))*_strength;
    return ret;
  }
  sizeType _N;
  scalarD _strength;
};
class TwoPlaneTerm : public SparseGEQConstraint
{
public:
  TwoPlaneTerm(const Optimizer& opt,sizeType N)
    :SparseGEQConstraint(opt),_N(N) {}
  virtual sizeType values() const {
    return 2;
  }
protected:
  virtual Vec valGradInner(const Vec& x,SMatd& grad) {
    grad.clear();
    addBLK(grad,0,IDV(_N-1,_N-1),Vec3d(1,0,1).transpose());
    addBLK(grad,1,IDV(0,_N-1),Vec3d(-1,0,1).transpose());

    Vec3d x0=Vec3d(7.5f,0,-2.0f);
    Vec ret=Vec::Zero(2);
    ret[0]=(x.segment<3>(IDV(_N-1,_N-1))-x0).dot(Vec3d(1,0,1));
    ret[1]=(x.segment<3>(IDV(0,_N-1))-x0).dot(Vec3d(-1,0,1));
    return ret;
  }
  sizeType _N;
};
class InBallTerm : public SparseGEQConstraint
{
public:
  InBallTerm(const Optimizer& opt,sizeType N):SparseGEQConstraint(opt),_N(N) {}
  virtual sizeType values() const {
    return 2;
  }
  bool isLinear() const {
    return false;
  }
protected:
  virtual Vec valGradInner(const Vec& x,SMatd& grad) {
    grad.clear();
    addBLK(grad,0,IDV(_N-1,_N-1),-2*(x.segment<3>(IDV(_N-1,_N-1))-Vec3d(_N-1,_N-1,0)).transpose());
    addBLK(grad,1,IDV(0,_N-1),-2*(x.segment<3>(IDV(0,_N-1))-Vec3d(0,_N-1,0)).transpose());

    Vec ret=Vec::Zero(2);
    ret[0]=0.1f-(x.segment<3>(IDV(_N-1,_N-1))-Vec3d(_N-1,_N-1,0)).squaredNorm();
    ret[1]=0.1f-(x.segment<3>(IDV(0,_N-1))-Vec3d(0,_N-1,0)).squaredNorm();
    return ret;
  }
  sizeType _N;
};
void OptimizerTest::runExample2(sizeType N,bool cons,bool consTwoPlane,bool consBall,const string& type)
{
  Optimizer opt;
  opt.getProps().put<int>("inputs",N*N*3);
  opt.getProps().put<string>("type",type);
  opt.addTerm(boost::shared_ptr<SpringTerm>(new SpringTerm(opt,N,100.0f)));
  opt.addTerm(boost::shared_ptr<GravityTerm>(new GravityTerm(opt,N,-9.81f)));
  if(cons)  //add soft constraint
    opt.addConstraint(boost::shared_ptr<FixPointTerm>(new FixPointTerm(opt,N,10.0f)));
  else  //add box hard constraint
    opt.addTerm(boost::shared_ptr<FixPointTerm>(new FixPointTerm(opt,N,10.0f)));
  if(consTwoPlane) //constraint above two plane points
    opt.addConstraint(boost::shared_ptr<TwoPlaneTerm>(new TwoPlaneTerm(opt,N)));
  if(consBall)  //constraint for points in a ball
    opt.addConstraint(boost::shared_ptr<InBallTerm>(new InBallTerm(opt,N)));
#undef ID
#undef IDV
#define ID(X,Y) (Y)*N+(X)
#define IDV(X,Y) ((Y)*N+(X))*3
  Vec x=Vec::Zero(N*N*3);
  for(sizeType i=0; i<N; i++)
    for(sizeType j=0; j<N; j++)
      x.segment<3>(IDV(i,j))=Vec3d(i,j,0);
  opt.debugFD();
  opt.minimize(x);

  //output
  ostringstream ss;
  ss << (cons ? "cons_" : "") << (consTwoPlane ? "twoPlane_" : "") << (consBall ? "ball_" : "") << type << ".vtk";
  VTKWriter<scalarD> wrt("cloth",ss.str(),true);
  vector<Vec3d,Eigen::aligned_allocator<Vec3d> > vss;
  vector<Vec2i,Eigen::aligned_allocator<Vec2i> > iss;
  for(sizeType i=0; i<N; i++)
    for(sizeType j=0; j<N; j++) {
      if(i>0)iss.push_back(Vec2i(ID(i,j),ID(i-1,j)));
      if(j>0)iss.push_back(Vec2i(ID(i,j),ID(i,j-1)));
      vss.push_back(x.segment<3>(IDV(i,j)));
    }
  wrt.appendPoints(vss.begin(),vss.end());
  wrt.appendCells(iss.begin(),iss.end(),VTKWriter<scalarD>::LINE);
}
#undef ID
#undef IDV
void OptimizerTest::runExample2()
{
  INFO("Unconstrained")
  runExample2(16,false,false,false,"LBFGS");
  runExample2(16,false,false,false,"BLEIC");
  runExample2(16,false,false,false,"AUGLAG");
  //These two guys doesn't work
  //runExample2(16,false,false,false,"LM");
  //runExample2(16,false,false,false,"ELM");

  INFO("Box Constraint")
  runExample2(16,true,false,false,"BLEIC");
  runExample2(16,true,false,false,"AUGLAG");

  INFO("Linear Constraint")
  runExample2(16,true,true,false,"BLEIC");
  runExample2(16,true,true,false,"AUGLAG");

  INFO("Nonlinear Constraint")
  runExample2(16,true,false,true,"AUGLAG");
}

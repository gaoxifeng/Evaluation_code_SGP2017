#include "ImplicitFunc.h"
#include "geom/StaticGeomCell.h"

USE_PRJ_NAMESPACE

//ImplicitFuncPlane
ImplicitFuncPlane::ImplicitFuncPlane() {}
ImplicitFuncPlane::ImplicitFuncPlane(const Vec3& x0,const Vec3& n):_p(x0,n) {}
ImplicitFuncPlane::ImplicitFuncPlane(const Vec3& a,const Vec3& b,const Vec3& c):_p(a,b,c) {}
scalar ImplicitFuncPlane::operator()(const Vec3& pos) const
{
  return _p.side(pos);
}
//ImplicitFuncCSG
ImplicitFuncCSG::ImplicitFuncCSG(OP_TYPE op):_a((ImplicitFuncCSG*)NULL),_b((ImplicitFuncCSG*)NULL),_alpha(0.8f),_op(op) {}
ImplicitFuncCSG::ImplicitFuncCSG(const BBox<scalar,3>& bb,OP_TYPE op):_alpha(0.8f),_op(op)
{
  BBox<scalar,2> bb2(bb._minC.block(0,0,2,1),bb._maxC.block(0,0,2,1));
  boost::shared_ptr<ImplicitFuncCSG> axis01(new ImplicitFuncCSG(bb2,op));
  boost::shared_ptr<ImplicitFuncCSG> axis2(new ImplicitFuncCSG(op));
  if(op == INTERSECT) {
    axis2->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._minC.z()),Vec3(0.0f,0.0f,-1.0f)));
    axis2->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._maxC.z()),Vec3(0.0f,0.0f, 1.0f)));
  } else {
    axis2->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._minC.z()),Vec3(0.0f,0.0f, 1.0f)));
    axis2->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._maxC.z()),Vec3(0.0f,0.0f,-1.0f)));
  }
  _a=axis01;
  _b=axis2;
}
ImplicitFuncCSG::ImplicitFuncCSG(const BBox<scalar,2>& bb,OP_TYPE op):_alpha(0.8f)
{
  boost::shared_ptr<ImplicitFuncCSG> axis0(new ImplicitFuncCSG(op));
  boost::shared_ptr<ImplicitFuncCSG> axis1(new ImplicitFuncCSG(op));
  if(op == INTERSECT) {
    axis0->_a.reset(new ImplicitFuncPlane(Vec3(bb._minC.x(),0.0f,0.0f),Vec3(-1.0f,0.0f,0.0f)));
    axis0->_b.reset(new ImplicitFuncPlane(Vec3(bb._maxC.x(),0.0f,0.0f),Vec3( 1.0f,0.0f,0.0f)));
    axis1->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._minC.y(),0.0f),Vec3(0.0f,-1.0f,0.0f)));
    axis1->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._maxC.y(),0.0f),Vec3(0.0f, 1.0f,0.0f)));
  } else {
    axis0->_a.reset(new ImplicitFuncPlane(Vec3(bb._minC.x(),0.0f,0.0f),Vec3( 1.0f,0.0f,0.0f)));
    axis0->_b.reset(new ImplicitFuncPlane(Vec3(bb._maxC.x(),0.0f,0.0f),Vec3(-1.0f,0.0f,0.0f)));
    axis1->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._minC.y(),0.0f),Vec3(0.0f, 1.0f,0.0f)));
    axis1->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._maxC.y(),0.0f),Vec3(0.0f,-1.0f,0.0f)));
  }
  _a=axis0;
  _b=axis1;
}
scalar ImplicitFuncCSG::operator()(const Vec3& pos) const
{
  if(!_a && !_b)
    return 1.0f;
  else if(!_b)
    return (*_a)(pos);
  else if(!_a) {
    scalar vb=(*_b)(pos);
    if(_op == SUBTRACT)
      vb*=-1.0f;
    return vb;
  } else {
    scalar va=(*_a)(pos);
    scalar vb=(*_b)(pos);
    if(_op == SUBTRACT)
      vb*=-1.0f;

    scalar sgn=_op == UNION ? -1.0f : 1.0f;
    return (va+vb+sgn*sqrt(abs(va*va+vb*vb-2.0f*va*vb*_alpha)))/(1.0f+_alpha);
  }
}
void ImplicitFuncCSG::setAlpha(const scalar& alpha)
{
  _alpha=alpha;
  ImplicitFuncCSG* a=dynamic_cast<ImplicitFuncCSG*>(_a.get());
  ImplicitFuncCSG* b=dynamic_cast<ImplicitFuncCSG*>(_b.get());
  if(a)a->setAlpha(alpha);
  if(b)b->setAlpha(alpha);
}
//ImplicitFuncGridRef
scalar ImplicitFuncGridRef::operator()(const Vec3& pos) const
{
  return _lsRef.sampleSafe(pos);
}
//ImplicitFuncReinit
ImplicitFuncReinit::ImplicitFuncReinit(scalar cellSz,const ImplicitFunc<scalar>& inner)
{
  BBox<scalar> bb=inner.getBB();
  bb.enlarged(cellSz*3,bb.getExtent()[2] == 0 ? 2 : 3);
  Vec3i nrCell=ceilV(Vec3(bb.getExtent()/cellSz));
  _ls.reset(nrCell,bb,0.0f);

  GridOp<scalar,scalar>::copyFromImplictFunc(_ls,inner);
  GridOp<scalar,scalar>::reinitialize(_ls);
}
ImplicitFuncReinit::ImplicitFuncReinit(const Grid<scalar,scalar> &tpl,const ImplicitFunc<scalar>& inner)
{
  _ls.makeSameGeometry(tpl);
  GridOp<scalar,scalar>::copyFromImplictFunc(_ls,inner);
  GridOp<scalar,scalar>::reinitialize(_ls);
}
scalar ImplicitFuncReinit::operator()(const Vec3& pos) const
{
  return _ls.sampleSafe(pos);
}
//ImplicitFuncMesh
ImplicitFuncMeshRef::ImplicitFuncMeshRef(const ObjMeshGeomCell& mesh)
  :_mesh(mesh),_ext(mesh.getBB().getExtent()*10) {}
scalar ImplicitFuncMeshRef::operator()(const Vec3& pos) const
{
  Vec3 n;
  bool inside=_mesh.closest(pos,n);
  scalar c=n.norm();
  if(_mesh.rayQuery(pos,-Vec3::Unit(0)*_ext[0]) == 1)
    return c;
  if(_mesh.rayQuery(pos, Vec3::Unit(0)*_ext[0]) == 1)
    return c;
  if(_mesh.rayQuery(pos,-Vec3::Unit(1)*_ext[1]) == 1)
    return c;
  if(_mesh.rayQuery(pos, Vec3::Unit(1)*_ext[1]) == 1)
    return c;
  if(_mesh.dim() == 3) {
    if(_mesh.rayQuery(pos,-Vec3::Unit(2)*_ext[2]) == 1)
      return c;
    if(_mesh.rayQuery(pos, Vec3::Unit(2)*_ext[2]) == 1)
      return c;
  }
  return c*(inside ? -1 : 1);
}
BBox<scalar> ImplicitFuncMeshRef::getBB() const
{
  return _mesh.getBB();
}
//ImplicitFuncRosy
ImplicitFuncRosy::ImplicitFuncRosy(const Vec3& origin,const Vec3& X,const scalar& step,const scalar& coef)
  :_step(step),_coef(coef),_origin(origin),_X(X) {}
scalar ImplicitFuncRosy::operator()(const Vec3& pos) const
{
  Vec3 rel=pos-_origin;
  return _coef*dist(Vec2(rel.dot(_X),(rel-rel.dot(_X)*_X).norm()));
}
scalar ImplicitFuncRosy::dist(const Vec2& p) const
{
  //scalar minX=p.x();
  scalar dist=ScalarUtil<scalar>::scalar_max;
  for(scalar currX=p.x();; currX-=_step) {
    scalar newDist=(Vec2(currX,y(currX))-p).norm();
    if(newDist < dist) {
      dist=newDist;
      //minX=currX;
    } else break;
  }
  for(scalar currX=p.x();; currX+=_step) {
    scalar newDist=(Vec2(currX,y(currX))-p).norm();
    if(newDist < dist) {
      dist=newDist;
      //minX=currX;
    } else break;
  }
  return (y(p.x()) < p.y()) ? dist : -dist;
}

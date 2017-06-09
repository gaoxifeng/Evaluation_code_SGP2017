#ifndef DHM_TRAITS_H
#define DHM_TRAITS_H

#include "CommonFile/MathBasic.h"
#include "CommonFile/solvers/MatVec.h"
#include <Eigen/Sparse>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

PRJ_BEGIN

//hash for eigen vector
struct Hash {
  size_t operator()(const sizeType& key) const {
    return _h(key);
  }
  size_t operator()(const Vec2i& key) const {
    return _h(key[0])+_h(key[1]);
  }
  size_t operator()(const Vec3i& key) const {
    return _h(key[0])+_h(key[1])+_h(key[2]);
  }
  size_t operator()(const Vec4i& key) const {
    return _h(key[0])+_h(key[1])+_h(key[2])+_h(key[3]);
  }
  size_t operator()(const Eigen::Matrix<sizeType,8,1>& key) const {
    return _h(key[0])+_h(key[1])+_h(key[2])+_h(key[3])+
           _h(key[4])+_h(key[5])+_h(key[6])+_h(key[7]);
  }
private:
  boost::hash<sizeType> _h;
};
//an implementation of multi_map
template <typename K,typename V,bool unique>
class multi_unordered_map
{
public:
  //keys
  typedef boost::fast_pool_allocator<std::pair<const K,sizeType> > ALLOC_KEYS;
  typedef typename boost::unordered_map<K,sizeType,Hash,std::equal_to<K>,ALLOC_KEYS> KEYS;
  typedef typename KEYS::const_iterator const_iterator;
  typedef typename KEYS::iterator iterator;
  //values
  typedef std::pair<V,sizeType> PAIR;
  typedef vector<V> SET;
  //methods
  const_iterator begin() const {
    return _keys.begin();
  }
  const_iterator end() const {
    return _keys.end();
  }
  const_iterator find(const K& key) const {
    return _keys.find(key);
  }
  void reserve(sizeType nr) {
    _keys.reserve(nr);
    _lists.reserve(nr);
  }
  //if unique=false, get all element IN ORDER
  void get(const K& key,SET& vss) const {
    vss.clear();
    const_iterator it=_keys.find(key);
    if(it == _keys.end())
      return;
    sizeType head=it->second;
    while(head != -1) {
      vss.push_back(_lists[head].first);
      head=_lists[head].second;
    }
    if(unique) {
      std::sort(vss.begin(),vss.end());
      typename SET::iterator last=std::unique(vss.begin(),vss.end());
      vss.erase(last,vss.end());
    }
  }
  //this method insert from beg to end IN ORDER
  template <typename ITER>
  void set(const K& key,ITER beg,ITER end) {
    if(beg==end) {
      _keys.erase(key);
    } else {
      _keys[key]=(sizeType)_lists.size();
      for(; beg!=end; ++beg)
        _lists.push_back(make_pair(*beg,(sizeType)_lists.size()+1));
      _lists.back().second=-1;
    }
  }
  //insert val to the back of list[key]
  void insert(const K& key,const V& val) {
    iterator it=_keys.find(key);
    if(it == _keys.end()) {
      _keys[key]=(sizeType)_lists.size();
      _lists.push_back(make_pair(val,-1));
    } else {
      sizeType head=it->second;
      while(_lists[head].second != -1)
        head=_lists[head].second;
      _lists[head].second=(sizeType)_lists.size();
      _lists.push_back(make_pair(val,-1));
    }
  }
protected:
  KEYS _keys;
  vector<PAIR,boost::fast_pool_allocator<PAIR> > _lists;
};
//basic types
struct DHMVertex;
struct DHMCell;
struct DHMFace;
class DHMTraits
{
public:
  //basic data structure
  typedef vector<boost::shared_ptr<DHMVertex>,boost::fast_pool_allocator<boost::shared_ptr<DHMVertex> > > VERTPTRS;
  typedef vector<boost::shared_ptr<DHMCell>,boost::fast_pool_allocator<boost::shared_ptr<DHMCell> > > CELLPTRS;
  typedef vector<boost::shared_ptr<DHMFace>,boost::fast_pool_allocator<boost::shared_ptr<DHMFace> > > FACEPTRS;
  //poses
  typedef vector<Vec3,boost::fast_pool_allocator<Vec3> > POSES;
  typedef boost::fast_pool_allocator<sizeType> ALLOC_POSSET;
  typedef boost::unordered_set<sizeType,Hash,std::equal_to<sizeType>,ALLOC_POSSET> POSSET;
  typedef POSSET::const_iterator CPOSITER;
  typedef POSSET::iterator POSITER;
  //cells
  typedef Eigen::Matrix<sizeType,8,1> CELLID;
  typedef vector<CELLID,boost::fast_pool_allocator<CELLID> > CELLIDS;
  //cell map
  typedef boost::fast_pool_allocator<std::pair<const CELLID,sizeType> > ALLOC_CMAP ;
  typedef boost::unordered_map<CELLID,sizeType,Hash,std::equal_to<CELLID>,ALLOC_CMAP> CMAP;
  typedef CMAP::const_iterator CCITER;
  typedef CMAP::iterator CITER;
  //texture coords
  typedef vector<Vec3i,boost::fast_pool_allocator<Vec3i> > TCOORDS;
  //face neighbor map
  typedef multi_unordered_map<sizeType,boost::shared_ptr<DHMFace>,true> FNMAP;
  //face id set
  typedef boost::fast_pool_allocator<Vec4i> ALLOC_FIDSET;
  typedef boost::unordered_set<Vec4i,Hash,std::equal_to<Vec4i>,ALLOC_FIDSET> FIDSET;
  //padding vertex
  typedef multi_unordered_map<Vec2i,sizeType,false> PMAP;
  typedef PMAP::const_iterator CPITER;
};
//basic types for energy
template <typename T,int rows=0>
class DHMEnergyTraits
{
public:
  static const int ROWS=rows;
  typedef Eigen::Matrix<T,ROWS,ROWS> HESS;
  typedef Eigen::Matrix<T,ROWS,1> GRAD;
  typedef Eigen::Matrix<T,ROWS,ROWS+1> BUNDLE;
  //Eigen column
  typedef Eigen::Matrix<T,-1,1> COL;
  typedef FixedSparseMatrix<T,Kernel<T> > SMAT;
  typedef vector<Eigen::Triplet<T,sizeType> > TRIPS;
};

PRJ_END

#endif

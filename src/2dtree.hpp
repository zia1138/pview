#ifndef __2dtree_hpp__
#define __2dtree_hpp__

#include <vector>
#include <limits>
#include <algorithm>

// This header file contains the main data structure, a kd-tree that
// works for k=2 or a 2d-tree.  

/* Implements in-place kd-tree idea from the following paper and master's
thesis:
Towards in-place geometric algorithms and data structures, with
1.) T. Chan and E. Chen. Proc. 20th Annu. ACM Symp. Comput. Geom.,
Brooklyn, 2004, pp. 239-246.
2.) Space-efficient geometric algorithms and data structures, NYU-Poly
Masters Thesis by Ilya Katz, December 2005.
*/

/*
  < class T > must implement the methods
  double x() const { return xval; }
  double y() const { return yval; }
  Implements a 2d-tree that stores values at the leaf nodes (some call this a 2d-trie).
  This 2d-tree correctly handles duplicate points. 
*/

using namespace std;

enum dim {x_axis, y_axis}; // Axis on which to conduct the split.
// Specifies the side of the spliting region left is less than right is greater than.
enum region_t { right_region, left_region };

// 2-d orthogonal region
template <class T> struct region {
  double x_min, x_max, y_min, y_max;

  // Simplest constructor. It just assigns max and min elements as specified.
  region(double x_min_, double x_max_, double y_min_, double y_max_){ 
    x_min = x_min_; x_max = x_max_;  y_min = y_min_; y_max = y_max_; 
  }
  // If no ranges are specified, the region spans the entire (x, y) space.
  region() { 
    x_min = -(numeric_limits<double>::max()); x_max = numeric_limits<double>::max();
    y_min = -(numeric_limits<double>::max()); y_max = numeric_limits<double>::max();
  }

  // Create a half plane region left, right, top or bottom.
  region(double splitline, dim splitdim, region_t r) {
    x_min = -(numeric_limits<double>::max()); x_max = numeric_limits<double>::max();
    y_min = -(numeric_limits<double>::max()); y_max = numeric_limits<double>::max();
    // Figure out which axis the node splits on.
    if(splitdim == x_axis){
      // Right if a "right" region region (greater than or equal to).
      if(r == right_region) x_min = splitline;
      else if( r == left_region )  x_max = splitline;
      // Left if a "left" region (strictly less than). 
    }
    else if(splitdim == y_axis){
      // Top if a "right" region region (greater than or equal to).
      if(r == right_region) y_min = splitline;
      else if( r == left_region )  y_max = splitline;
      // Bottom if a "left" region (strictly less than). 
    }
  }

  // Returns true if the current region encloses the passed region R. 
  bool encloses( region R ){ 
    return ( (R.x_max <= x_max) && (R.x_min >= x_min) ) && 
           ( (R.y_max <= y_max) && (R.y_min >= y_min) );
  }

  // Returns true if the current region contains a point (or pointer to a point).  
  bool encloses( T p ) { return ( (p.x()  <= x_max) && (p.x()  >= x_min) ) && ( (p.y()  <= y_max) && (p.y()  >= y_min) ); }
  bool encloses( T *p ){ return ( (p->x() <= x_max) && (p->x() >= x_min) ) && ( (p->y() <= y_max) && (p->y() >= y_min) ); }

  // Returns true of the current region and R intersect. 
  bool intersects( region &R ){    
    // Negation of (not intersection). 
    return !( R.x_max < x_min || R.x_min > x_max) && !(R.y_max < y_min || R.y_min > y_max ); 
  }
  // Constructor creates a new region that is the intersection of the two given regions.
  region ( region a, region b) {
    x_max = min(a.x_max, b.x_max); x_min = max(a.x_min, b.x_min);
    y_max = min(a.y_max, b.y_max); y_min = max(a.y_min, b.y_min);
  }
};

template <class T> struct _2dtree_inplace {
  struct x_cmp { bool operator () (T a, T b) const { return a.x() < b.x(); } };
  struct y_cmp { bool operator () (T a, T b) const { return a.y() < b.y(); } };
  typedef typename vector<T>::iterator it_t;

  void build(it_t first, it_t last) {
    if(last - first <= 3) return; // <= 3, since need mid, Lmid, and Rmid

    it_t mid  = first + ( last - first + 1 ) / 2;
    it_t Lmid = first + ( mid - first + 1 ) / 2;
    it_t Rmid = mid   + ( last - mid + 1 ) / 2;
    
    // Note that this code relies heavily on the orthogonal range
    // queries being inclusive of the borders [x1,x2] x [y1,y2]
    // instead of (x1,x2) x (y1, y2). This assures that intersections
    // will handle repeats correctly.
    nth_element(first, mid, last,  x_cmp());
    nth_element(first, Lmid, mid,  y_cmp()); // split left half by y-coordniate
    nth_element(mid+1, Rmid, last, y_cmp()); // split right half by y-coordinate

    /// Use Lmid + 1, mid+1, and Rmid+1 because they become spliting
    /// plane elements.
    build(first, Lmid);  // recursively partition Lbot
    build(Lmid+1, mid);  // recursively partition Ltop
    build(mid+1, Rmid);  // recursively partition Rbot
    build(Rmid+1, last); // recursively partition Rtop
  }

  vector<T> &P;
  region<T> region0;

  _2dtree_inplace(vector<T> &P_, bool skip_build = false) : P(P_) {   
    // reverse region0
    region0.x_min = numeric_limits<double>::max(); region0.x_max = -(numeric_limits<double>::max());
    region0.y_min = numeric_limits<double>::max(); region0.y_max = -(numeric_limits<double>::max());
    for(size_t i = 0; i < P.size(); i++) {
      double x = P[i].x(), y = P[i].y();
      if(x < region0.x_min) region0.x_min = x; if(y < region0.y_min) region0.y_min = y;
      if(x > region0.x_max) region0.x_max = x; if(y > region0.y_max) region0.y_max = y;
    }
    if(skip_build == false) build(P.begin(), P.end());   
  }

  void query(it_t first, it_t last, region<T> vr, region<T> R, vector<T*> &res)  {
    if(R.encloses(vr)) { // Region is entirely enclosed, report everyone.
      for(it_t it = first; it != last; it++) res.push_back(&*it);
      return;
    }
    if( last - first <= 3 ) { // <= 3, since need mid, Lmid, and Rmid
      for(it_t it = first; it != last; it++) if(R.encloses(*it)) res.push_back(&*it);
      return;
    }
    it_t mid  = first + ( last - first + 1 ) / 2;
    it_t Lmid = first + ( mid - first + 1 ) / 2;
    it_t Rmid = mid   + ( last - mid + 1 ) / 2;

    // Create regions based on specified splitting planes.
    region<T> left( region<T>((*mid).x(), x_axis, left_region),  vr);
    region<T> right(region<T>((*mid).x(), x_axis, right_region), vr);

    region<T> Lbot(region<T>((*Lmid).y(), y_axis, left_region), left);
    region<T> Ltop(region<T>((*Lmid).y(), y_axis, right_region), left);

    region<T> Rbot(region<T>((*Rmid).y(), y_axis, left_region), right);
    region<T> Rtop(region<T>((*Rmid).y(), y_axis, right_region), right);

    if(R.intersects(Lbot)) query(first, Lmid, Lbot, R, res);
    if(R.intersects(Ltop)) query(Lmid+1, mid, Ltop, R, res);
    if(R.intersects(Rbot)) query(mid+1, Rmid, Rbot, R, res);
    if(R.intersects(Rtop)) query(Rmid+1, last, Rtop, R, res);
    
    if(R.encloses(*mid))  res.push_back(&*mid);
    if(R.encloses(*Rmid)) res.push_back(&*Rmid);
    if(R.encloses(*Lmid)) res.push_back(&*Lmid);
  }

  void Query ( region< T > R, vector<T*> &res) { query(P.begin(), P.end(), region0, R, res);  }

  // Find point nearest in retention time axis within a given query rectangle.
  T* NearestX(double x, double y, double dx, double dy) {
    T *nearest = NULL;
    vector<T *> in_region;
    region<T> query(x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2);
    Query(query, in_region);
    if(in_region.size() > 0) {
      double xmin = numeric_limits<double>::max();
      for(size_t r = 0; r < in_region.size(); r++) {
	if(fabs(in_region[r]->x() - x) < xmin) {
	  xmin = fabs(in_region[r]->x() - x);
	  nearest = in_region[r];
	}
      }
    }
    return nearest;
  }

};

// Pointer interface for 2d-tree.  TODO: Use this pointer interface to
// eliminate repeated code.
template <class T> struct _2dtree_ptr {
  T *ptr;
  double x() { return ptr->x(); }  double y() { return ptr->y(); }
};

// Implements a 2d-tree that uses pointers to the data. This avoids
// repeating all of the code in 2dtree_inplace!
template <class T> struct _2dtree {
  _2dtree_inplace< _2dtree_ptr<T> > *root;
  vector< _2dtree_ptr<T> > P;

  _2dtree(vector<T> &P_) {   
    P.resize(P_.size());
    for(size_t i = 0; i < P_.size(); i++) P[i].ptr = &P_[i];
    root = new _2dtree_inplace< _2dtree_ptr<T> >(P);
  }
  ~_2dtree() { delete root; }

  void Query ( region< T > R, vector<T*> &res) { 
    vector< _2dtree_ptr<T> *> res2;
    region< _2dtree_ptr<T> > R2(R.x_min, R.x_max, R.y_min, R.y_max);
    root->Query(R2, res2);
    for(size_t i = 0; i < res2.size(); i++) res.push_back(res2[i]->ptr);
  }
  T* NearestX(double x, double y, double dx, double dy) { 
    _2dtree_ptr<T> *nearest = root->NearestX(x,y,dx,dy); 
    return nearest == NULL ? NULL : nearest->ptr;
  }
};

#endif //  __2dtree_hpp__

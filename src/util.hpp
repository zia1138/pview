#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <vector>
#include <stack>
#include <iomanip>
#include <sstream> 

// Utility code. Generally useful functions, that don't belong anywhere else.
namespace util {
  using namespace std;

  // Generates all multi-combinations of n elements taken t elements
  // at a time *WITH REPETITIONS*.  NOTE: There is probably a faster
  // algorithm involving Gray codes.
  struct multi_comb {
    vector<int> m; // current muli-combination
    vector< vector<int> > &res;
    multi_comb(int n, int t, vector< vector<int> > &res_) : m(t+1), res(res_) {
      m[t] = n-1;
      run(n, t, t-1);
    }
    // Recursive procedure generates multi-combinations in lexicographical order.
    // m[0] <= m[1] <= m[2] <= ... <= m[t-1] <= n-1
    void run(int n, int t, int j) {
      if(j < 0) {
	vector<int> r(t);
	for(int i = 0; i < t; i++) r[i] = m[i];	
	res.push_back(r);
	return;
      }
      for(int i = 0; i <= m[j+1]; i++) { m[j] = i; run(n, t, j-1); }
    }
  };

  // Simple adjacency list undirected graph, useful for computing
  // connected components for designating protein groups.
  struct ugraph_t {
    vector< vector<int> > A;
    ugraph_t(size_t N) : A(N) {}
    // Does the edge exists?
    bool edge_exists(int i, int j) {
      bool ij = false, ji = false;
      for(size_t k = 0; k < A[i].size(); k++) if(A[i][k] == j) ij = true;
      for(size_t k = 0; k < A[j].size(); k++) if(A[j][k] == i) ji = true;
      if(ij != ji) cerr << "ij ji problem" << endl;
      return ij;
    }
    // No self-edges and no parallel edges.
    void add_edge(int i, int j) {
      if(i == j || edge_exists(i,j)) return;
      A[i].push_back(j); A[j].push_back(i);
    }
    // Get adjacent nodes.
    vector<int> &adjacent(int c) { return A[c]; }
    // Depth first search for finding connected components.
    void dfs_label(vector<int> &labels, int cur_label, int start) {
      stack<int> L; L.push(start);
      while(L.empty() == false) {
	// Get current node from stack.
	int cur = L.top(); L.pop();
	labels[cur] = cur_label; // Label node.
	// Get adjacent nodes.
	vector<int> &adj = adjacent(cur); 
	for(size_t j = 0; j < adj.size(); j++) {
	  // Nodes that haven't been labeled, add them to the stack.
	  if(labels[adj[j]] == 0) { 
	    labels[adj[j]] = 1; L.push(adj[j]); 
	  }
	}
      }
    }
    // Simple DFS based algorithm for computing connected components.
    // Assigns each node i in label[i] a label according component
    // membership. Returns total number of components.
    int connected_components(vector<int> &labels) {
      int num_components = 0;
      labels.resize(A.size()); for(size_t i = 0;  i < labels.size(); i++) labels[i] = 0;
      for(size_t i = 0;  i < labels.size(); i++) {
	if(labels[i] != 0) continue; // Already visited vertex, skip.
	// Apply DFS to label component.
	dfs_label(labels, num_components + 2, i);
	num_components++; // Advance to next component label.
      }
      // Adjust for 0 = visit 1=traverse labelling.
      for(size_t i = 0; i < labels.size(); i++) labels[i] -= 2;
      return num_components; 
    }
  };

  // Enumerates cross product of given ranges.  0..L[0] x 0..L[1] x
  // 0..L[2] x ... x 0..L[n-1] with some entries constrained as
  // combinations.
  struct cross {
    vector<int> &L, c;
    vector< vector<int> > &res;
    vector< bool> backref;
    cross(vector<int> &L, vector<bool> &backref, vector< vector<int> > &res) : 
      L(L), c(L.size(), 0), res(res), backref(backref) {
      run(0); 
    }
    cross(vector<int> &L, vector< vector<int> > &res) : 
      L(L), c(L.size(), 0), res(res), backref(L.size(), false) {
      run(0);
    }
    // Compute cross product recursively.
    void run(int p) {
      if(p == (int)L.size()) {
	vector<int> r(L.size());
	for(size_t i = 0; i < L.size(); i++) r[i] = c[i];	
	res.push_back(r);
 	return;
      }
      // If backref is true start from the value at (c[p-1]+1).  This
      // constrains some entries as comabinations.
      if(backref[p] && p > 0)
	for(int i = c[p-1]+1; i < L[p]; i++) { c[p] = i; run(p+1); }
      else
	for(int i = 0; i < L[p]; i++) { c[p] = i; run(p+1); }
    }
  };

  template<typename T> string toString(const T& t, bool *ok = NULL) {
    ostringstream stream;
    stream << t;
    if(ok != NULL) *ok = stream.fail() == false;
    return stream.str();
  }

  inline string toStringDouble(double t, int prec, bool *ok = NULL) {
    ostringstream stream;
    stream << fixed << setprecision(prec)  << t;
    if(ok != NULL) *ok = stream.fail() == false;
    return stream.str();
  }

  template<typename T> T fromString(const string& s, bool *ok = NULL) {
    istringstream stream (s);
    T t;
    stream >> t;
    if(ok != NULL) *ok = stream.fail() == false;
    return t;
  }  

  template <typename T> T fromString(const string& s, T &defval) {
    bool ok;
    T val = fromString<T>(s, &ok);
    if(!ok) return defval; else return val;
  }

  // Note median actually modifies the order of the elements in the
  // given vector!! This happens because nth_element's O(N) selection
  // algorithm will swap elements around.
  inline double _median(vector<double>::iterator begin, vector<double>::iterator end) {
    size_t N = end - begin;
    if(N == 1) return *begin;
    if(N % 2 == 0) { // even number of elements
      vector<double>::iterator mid = begin + ( end - begin ) / 2 - 1;
      nth_element(begin, mid, end);
      vector<double>::iterator mid1 = begin + ( end - begin ) / 2;	
      nth_element(begin, mid1, end);
      return *mid + (*mid1 - *mid) / 2.0;
    }
    else { // odd number of elements
      vector<double>::iterator mid = begin + ( end - begin ) / 2;
      nth_element(begin, mid, end);
      return *mid;
    }
  }

  // This version of median will change the order of elements.  
  inline double median_unsafe(vector<double> &data) {
    if(data.size() == 0) { std::cerr << "no data for median()" << std::endl; exit(1); }
    return _median(data.begin(), data.end());
  }

  // Computes the median in O(N) time. By making a copy this version assures no re-ordering.
  inline double median(vector<double> &data) {
    vector<double> data2(data.size());
    for(size_t i = 0; i < data.size(); i++) data2[i] = data[i];
    return median_unsafe(data2);
  }

  // O(KN) running median computation (can be improved to O(log K * N)
  // with a min-heap and max-heap.
  inline void runmed(vector<double> &x, vector<double> &m) {
    double n = x.size();
    long k = (long)min( ( n - 1 ) / 2, ceil(0.1*n) ); // Automatically estimate K for runs.
    m.resize(x.size());
    long k2 = k / 2;
    vector<double> x2(k);
    // Advance along vector.
    for(long i = 0; i < (long)x.size(); i++) {
      long start = max(0L, i - k2);
      long end = min((long)x.size(), i + k2 + 1);
      x2.resize(end - start);
      // Copy run into another vector (becaused _median modifies
      // original order).
      int j = 0;
      while(start < end) { x2[j] = x[start]; j++; start++; }
      // Compute median using that vector. 
      m[i] = _median(x2.begin(), x2.end());
    }
  }

  // There's no log2() in Windows. 
  inline double Log2(double x) { return log(x)/log(double(2)); }
};



#endif // __UTIL_HPP__

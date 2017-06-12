#ifndef FINDER_H_
#define FINDER_H_


#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>

#include "patricia.h"

using std::set;
using std::map;
using std::pair;
using std::string;
using std::vector;
using std::greater;
using std::unordered_map;

struct IntPair {
public:
  int g;
  int s;
  int rank;
};

inline bool operator< (const IntPair& left, const IntPair& right){
  return (left.rank < right.rank);
}

struct Counterpart {
public:
  int target;
  int rank;
};

inline bool operator< (const Counterpart& left, const Counterpart& right){
  return (left.rank < right.rank);
}

struct Triplet {
public:
  int x, y, z;
  Triplet reverse();
  explicit Triplet(){}
  explicit Triplet(int _x,int _y,int _z): x(_x),y(_y),z(_z){}
};

struct Pair { public: int a, b; void set(int, int); };
inline void Pair::set(int _a, int _b){ a = _a; b = _b; }

struct VertexPair: public Pair { public: int id; };
struct Edge   { public: int to, id; Triplet labels; };
struct DFSCode { public: Triplet labels; Pair time; };

typedef vector< vector<Edge> > AdjacentList;

class Graph: public AdjacentList {
 public:
  float value;
  string name;
  vector<int> label;
  int num_of_edges;
  void resize(unsigned int s){ AdjacentList::resize(s); label.resize(s); }
};

class EdgeTracer {
 public:
  VertexPair  vpair;
  EdgeTracer* predec;
  explicit EdgeTracer(){};
  void set(int,int,int,EdgeTracer*);
};

inline void EdgeTracer::set(int a,int b,int id, EdgeTracer* pr){
  vpair.a = a;
  vpair.b = b;
  vpair.id = id;
  predec = pr;
}

typedef vector<EdgeTracer> Tracers;
typedef map<int,Tracers> GraphToTracers;
typedef map<Pair,GraphToTracers> PairSorter;

inline Triplet Triplet::reverse(){ return Triplet(z,y,x); }

inline bool operator< (const Triplet& left, const Triplet& right){
  if (left.x!=-1 && right.x!=-1 && left.x != right.x) return (left.x < right.x);
  if (left.y!=-1 && right.y!=-1 && left.y != right.y) return (left.y < right.y);
  return (left.z < right.z);
}

inline bool operator<= (const Triplet& left, const Triplet& right){
  return !(right < left);
}

inline bool operator< (const Pair& left, const Pair& right){
  if (left.a != right.a) return (left.a < right.a);
  return (left.b < right.b);
}

inline bool operator== (const DFSCode& left, const DFSCode& right){
  if(left.time.a != right.time.a) return false;
  if(left.time.b != right.time.b) return false;	
  if(left.labels.x != right.labels.x) return false;
  if(left.labels.y != right.labels.y) return false;
  return (left.labels.z == right.labels.z);	
}

inline bool operator!= (const DFSCode& x, const DFSCode& y){
  return !(x==y);
}

inline std::ostream& operator<< (std::ostream& os, const vector<DFSCode> pattern){
  if(pattern.empty()) return os;
  os << "(" << pattern[0].labels.x << ") " << pattern[0].labels.y << " (0f" << pattern[0].labels.z << ")";
  for(unsigned i=1; i<pattern.size(); ++i){
    if(pattern[i].time.a < pattern[i].time.b){
      os << " " << pattern[i].labels.y << " (" << pattern[i].time.a << "f" << pattern[i].labels.z << ")";
    }else{
      os << " " << pattern[i].labels.y << " (b" << pattern[i].time.b << ")";
    }
  }
  return os;
}

typedef map<int,vector<int> > occurence;

struct datapack {
  int cmp;
  int seq;
  int label;
};
/*
bool operator<(const Pair &a, const Pair &b){
  if(a.a < b.a){
    return true;
  }else if(a.a==b.a && a.b < b.b){
    return true;
  }
  return false;
  }*/
class Finder {
 public:
  unordered_map<int,string> g_tname;
  unordered_map<string,int> g_tid;
  map<int,string> feature;
  unordered_map<string,int> fdic;
  unordered_map<int,float> coeff;
  float alphasum, bias;
  vector<Graph>  gdata;
  vector<datapack>  trans;
  map<Triplet,GraphToTracers> gheap;
  patricia gpat_tree;
  patricia spat_tree;
  map<Pair,float> cooceff;
  vector<vector<string> > fvec;
  // member functions
  template<typename T> void read_graphs(T& ins);
  template<typename T> void read_features(T& ins);
};

template<typename T> void Finder::read_features(T& ins){
  alphasum = 0.0;
  bias = 0.0;

  string line;
  vector<string> vg;
  int gcount = 1;
  std::string tmp, dfscode;
  float val;
  while( getline(ins,line) ){
    std::list<std::string> vec;
    boost::split(vec, line, boost::is_any_of("\t"));
    tmp = vec.front();
    vec.pop_front();
    val = atof(tmp.c_str());
    bias -= val;
    alphasum += std::abs(val);
    //std::cout<<vec.size()<<std::endl;
    if(vec.size()==1){
      dfscode = boost::algorithm::join(vec, " ");
      fvec.resize(fvec.size()+1);
      fvec[fvec.size()-1].push_back(dfscode);
      vg.clear();
      if(fdic.find(dfscode)==fdic.end()){
	feature[gcount] = dfscode;
	fdic[dfscode]   = gcount;
	coeff[gcount]   = val;
	gtokenize(vg,dfscode);
	gpat_tree.insert(vg,patricia::root,gcount);
	gcount += 1;
      }else{
	int fid = fdic[dfscode];
	coeff[fid] += val;
      }
    }else{
      vector<int> cooc;
      cooc.resize(vec.size());
      fvec.resize(fvec.size()+1);
      string pat[1];
      int i = 0;
      while(vec.size()>0){
	pat[0] = vec.front();
	dfscode = pat[0];//boost::algorithm::join(pat, " ");
	//std::cout<<dfscode<<std::endl;
	fvec[fvec.size()-1].push_back(dfscode);
	vg.clear();
	if(fdic.find(dfscode)==fdic.end()){
	  feature[gcount] = dfscode;
	  fdic[dfscode]   = gcount;
	  coeff[gcount]   = 0;
	  gtokenize(vg,dfscode);
	  gpat_tree.insert(vg,patricia::root,gcount);
	  cooc[i]=gcount;
	  gcount += 1;
	}else{
	int fid = fdic[dfscode];
	cooc[i]=fid;
	}
	++i;
	vec.pop_front();
      }
      Pair p;
      if(cooc[0] >= cooc[1]){
	p.set(cooc[1],cooc[0]);
      }else{
	p.set(cooc[0],cooc[1]);
      }
      if(cooceff.find(p)==cooceff.end()){
	cooceff[p]=val;
      }else{
	std::cout << "!!!!!!!"<<val<<std::endl;
	cooceff[p]+=val;
      }
    }
  }
  bias /= alphasum;
  std::cout << "#features: " << gcount-1 << std::endl;
  std::cout << "#features: " << gpat_tree.node.size() << std::endl;
}

template<typename T>
void Finder::read_graphs(T& ins){
  Graph g;
  Triplet labels;
  Triplet wlabels;
  Edge edge;

  float val;
  char c; unsigned int i, j;
  std::string line, gname;
  int eid=0;
  int gid=0;

  EdgeTracer cursor;

  bool reading = false;
  while (getline(ins, line)) {
    if (line.empty()) {
      g.num_of_edges = eid;
      gdata.push_back(g);
      gid++;
      reading = false;
    }
    std::stringstream stream(line);
    if (line[0] == 't') {
      reading = true;
      g.clear();
      eid = 0;
      // t # 0 -1 name
      stream >> c >> c >> c >> val >> gname;
      //std::cout << gid << ") " << gname << " " << val << std::endl;
      g.name = gname;
      g.value = val;
    } else if (line[0] == 'v') {
      stream >> c >> i >> j;
      g.resize(i+1);
      g.label[i] = j;
    } else if (line[0] == 'e') {
      stream >> c >> i >> j >> labels.y;
      labels.x = g.label[i];
      labels.z = g.label[j];
      int etmp = eid;
      edge.to = j; edge.labels = labels; edge.id = eid++;
      g[i].push_back(edge);
      edge.to = i; edge.labels = labels.reverse(); edge.id = eid++;
      g[j].push_back(edge);
      // edge i->j
      cursor.vpair.a = i;
      cursor.vpair.b = j;
      cursor.vpair.id = etmp;
      cursor.predec = 0;
      gheap[labels][gid].push_back(cursor);
      wlabels = labels;
      wlabels.z = 999;
      gheap[wlabels][gid].push_back(cursor);
      // edge j->i
      cursor.vpair.a = j;
      cursor.vpair.b = i;
      cursor.vpair.id = etmp+1;
      cursor.predec = 0;
      gheap[labels.reverse()][gid].push_back(cursor);
      wlabels = labels.reverse();
      wlabels.z = 999;
      gheap[wlabels][gid].push_back(cursor);
    }
  }
  if(reading==true){
    g.num_of_edges = eid;
    gdata.push_back(g);
    gid++;
    reading = false;
  }
}

int buildDFSCode(const string& str, vector<DFSCode>& query, int _cur, map<int,int>& label, vector<int>& addkeys){

  int cur = _cur;
  DFSCode tmp;
  std::istringstream in(str);
  char c;
  int i, j, elab=-200;
  while(in.get(c)){
    string token(1,c);
    unsigned int n = 1;
    while(in.get(c) && c != ' '){
      token += c;
      n += 1;
    }
    if (token[0] == '(' && token[n-1] == ')'){
      string sstr = token.substr(1,n-2);
      string::size_type r;
      if ((r = sstr.find("f")) != string::npos){
	i = atoi(sstr.substr(0,r).c_str());
	j = atoi(sstr.substr(r+1,string::npos).c_str());
	//std::cout << "call label[" << cur << "]=" << j << std::endl;
	addkeys.push_back(cur);
	label[cur] = j;
	tmp.time.a = i;
	tmp.time.b = cur;
	tmp.labels.x = label[i];
	tmp.labels.y = elab;
	tmp.labels.z = j;
	query.push_back(tmp);
	cur += 1;
      }else if ((r = sstr.find("b")) != string::npos){
	i = atoi(sstr.substr(r+1,string::npos).c_str());
	tmp.time.a = cur-1;
	tmp.time.b = i;
	tmp.labels.x = label[cur-1];
	tmp.labels.y = elab;
	tmp.labels.z = label[i];
	query.push_back(tmp);
      }else{
	i = atoi(sstr.c_str());
	label[cur] = i;
	addkeys.push_back(cur);
	cur += 1;
      }
    }else{
      elab = atoi(token.c_str());
    }
  }
  //std::cout << "done" << std::endl;
  
  return cur;
}

#endif /*FINDER_H_*/

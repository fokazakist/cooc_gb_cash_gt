#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <unordered_map>
#include <boost/xpressive/xpressive.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "finder.h"

using std::deque;
using std::vector;
using std::string;

using std::istringstream;

using std::unordered_map;

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#define USAGE " features graphs"

void gtraverse(int node_id, vector<DFSCode>& query, unsigned pnt, int cur,
	       GraphToTracers& g2tracers, vector<Graph>& gdata,
	       const patricia& P, map<int,set<int> >& saver, int treelevel,
	       map<int,int>& dfslabel){

  //if (g2tracers.size()==0) return;

  if (query.size() > pnt){
    DFSCode& dfscur = query[pnt];
    int s = dfscur.time.a;
    int t = dfscur.time.b;
    EdgeTracer* tracer;
    EdgeTracer cursor;
    GraphToTracers newtracer;

    for(GraphToTracers::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x){
      int gid = x->first;
      Graph& g = gdata[gid];
      for(vector<EdgeTracer>::iterator it = x->second.begin(); it != x->second.end(); ++it){
	
	tracer = &(*it);

	vector<unsigned> discovered(g.size(),0);
	int growpoint = -999;
	int backpoint = -999;
	for(int j = pnt-1; j>=0; --j, tracer = tracer->predec){
	  discovered[tracer->vpair.a] = 1;
	  discovered[tracer->vpair.b] = 1;
	  if(query[j].time.a == s) growpoint = tracer->vpair.a;
	  if(query[j].time.b == s) growpoint = tracer->vpair.b;
	  if(s > t){
	    if(query[j].time.a == t) backpoint = tracer->vpair.a;
	    if(query[j].time.b == t) backpoint = tracer->vpair.b;
	  }
	}
	
	for(unsigned int i=0; i<g[growpoint].size(); ++i){
	  const Edge& e = g[growpoint][i];
	  if(dfscur.labels.z!=999){
	    if (e.labels.y == dfscur.labels.y && e.labels.z == dfscur.labels.z){
	      if (s < t){
		if (discovered[e.to] != 1){
		  cursor.vpair.a = growpoint;
		  cursor.vpair.b = e.to;
		  cursor.vpair.id = e.id;
		  cursor.predec  = &(*it);
		  newtracer[gid].push_back(cursor);
		}
	      }else{
		if (backpoint == e.to){
		  cursor.vpair.a = growpoint;
		  cursor.vpair.b = e.to;
		  cursor.vpair.id = e.id;
		  cursor.predec  = &(*it);
		  newtracer[gid].push_back(cursor);
		}
	      }
	    }
	  }else{
	    if (e.labels.y == dfscur.labels.y){
	      if (s < t){
		if (discovered[e.to] != 1){
		  cursor.vpair.a = growpoint;
		  cursor.vpair.b = e.to;
		  cursor.vpair.id = e.id;
		  cursor.predec  = &(*it);
		  newtracer[gid].push_back(cursor);
		}
	      }else{
		if (backpoint == e.to){
		  cursor.vpair.a = growpoint;
		  cursor.vpair.b = e.to;
		  cursor.vpair.id = e.id;
		  cursor.predec  = &(*it);
		  newtracer[gid].push_back(cursor);
		}
	      }
	    }
	  }
	}
      }
    }
    gtraverse(node_id,query,pnt+1,cur,newtracer,gdata,P,saver,treelevel,dfslabel);
  }else{

    if(P.outQ[node_id]>0){
      //std::cout << "==>" << g2tracers.size() << " " << query << std::endl;

      for(GraphToTracers::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x){
	saver[x->first].insert(P.outQ[node_id]);
	//saver[P.outQ[node_id]].insert(x->first);
      }
    }
    
    const list<int>& child = P.children[node_id];

    for(list<int>::const_iterator it=child.begin(); it != child.end(); ++it){
      const vector<string>& vec = P.node[*it];
      int addcount = 0;
      int c2 = cur;
      vector<int> addkeys;
      for(vector<string>::const_iterator it2=vec.begin(); it2 != vec.end(); ++it2){
	c2 = buildDFSCode(*it2, query, c2, dfslabel, addkeys);
	addcount += 1;
      }
      gtraverse(*it,query,pnt,c2,g2tracers,gdata,P,saver,treelevel+1,dfslabel);
      for(vector<int>::iterator lit=addkeys.begin(); lit!=addkeys.end(); ++lit){
	//std::cout << "erase" << *lit << std::endl;
	dfslabel.erase(*lit);
      }
      for(int j=0; j<addcount; ++j){
	query.pop_back();
      }
    }
  }
}


int main(int argc, char **argv) {

  int opt;
  string outdir=".";
  string filename = "temp.txt";
  while ((opt = getopt(argc, argv, "d:f:o:")) != -1) {
    switch (opt) {
    case 'f':
      filename = string(optarg);
      break;
    case 'o':
      outdir = string(optarg);
      break;
    default:
      std::cerr << "Usage: "<< argv[0] << USAGE << std::endl;
      return -1;
    }
  }
  
  // *********** i/o check ***********

  if(argc-optind != 2){
    std::cerr << "Usage: "<< argv[0] << USAGE << std::endl;
    return -1;
  }

  // *********** mode check ***********

  Finder f;

  string filefeat   (argv[optind]);
  string filegraph  (argv[optind+1]);

  // *********** file check ***********

  std::ifstream feat_file(filefeat.c_str());
  if(feat_file.fail()){
    std::cerr << "File not found: " << filefeat << std::endl;
    return -1;
  }
  std::ifstream graph_file(filegraph.c_str());
  if(graph_file.fail()){
    std::cerr << "File not found: " << filegraph << std::endl;
    return -1;
  }

  std::cout << "Finder version 0.4->1.0"         << std::endl;
  std::cout << "      Okazaki" << std::endl;
  std::cout << "  Feature => " << filefeat   << std::endl;
  std::cout << "  Graph   => " << filegraph  << std::endl;

  // *********** main ***********

   f.read_features( feat_file );
   f.read_graphs(graph_file);

   std::cout << "# graphs: " << f.gdata.size() << std::endl;

   map<int,set<int> > gsave;

   vector<DFSCode> gquery;
   map<int,int> dfslabel;

   list<int>& childg = f.gpat_tree.children[patricia::root];
   for(list<int>::iterator it=childg.begin(); it != childg.end(); ++it){
     vector<string>& vec = f.gpat_tree.node[*it];
     int cur = 0;
     int addcount = 0;
     vector<int> addkeys;
     for(vector<string>::iterator it2=vec.begin(); it2 != vec.end(); ++it2){
       cur = buildDFSCode(*it2, gquery, cur, dfslabel, addkeys);
       addcount += 1;
     }
     GraphToTracers& newtracer = f.gheap[gquery[0].labels];
     gtraverse(*it,gquery,1,cur,newtracer,f.gdata,f.gpat_tree,gsave,1,dfslabel);
     for(vector<int>::iterator lit=addkeys.begin(); lit!=addkeys.end(); ++lit){
       //std::cout << "erase" << *lit << std::endl;
       dfslabel.erase(*lit);
     }
     for(int j=0; j<addcount; ++j){
       gquery.pop_back();
     }
   }
   std::ofstream ofs(filename);
   ////std::cout << "gsave size: " << gsave.size() << std::endl;
   /*
   std::ofstream ofs1("m1");

   for(vector<vector<string> >::iterator p = f.fvec.begin();p !=f.fvec.end();++p){
     if((*p).size()==1){
       ofs1 <<(*p)[0]<<"\t";
       int fid = f.fdic[(*p)[0]];
       //std::cout << fid<<std::endl;
       for(unsigned int i=0;i<f.gdata.size();++i){
	 if(gsave[i].count(fid)==1){
	   ofs1 << i<<",";
	 }
       }
       ofs1 << std::endl;
     }else{
       ofs1 <<(*p)[0]<<"\t"<<(*p)[1]<<"\t";
       int fid1,fid2;
       fid1 = f.fdic[(*p)[0]];
       fid2 = f.fdic[(*p)[1]];
       
       for(unsigned int i=0;i<f.gdata.size();++i){
	 if(gsave[i].count(fid1)==1){
	   if(gsave[i].count(fid2)==1){
	     ofs1 << i<<",";
	   }
	 }
       }
       ofs1 << std::endl;
       
     }
   }
   */
   float cor = 0, err = 0;
   
   for(unsigned i=0; i<f.gdata.size(); ++i){
     std::cout << "GRAPH #" << i+1 << " (" << gsave[i].size() << ")" << std::endl;
     float sum = f.bias;
     //std::cout << "   bias=" << sum << std::endl; 
     
     for(set<int>::iterator itr = gsave[i].begin(); itr != gsave[i].end(); ++itr){
       std::cout << "   alpha[i]=" << f.coeff[*itr] << " i=" << *itr << std::endl; 
       sum += 2*f.coeff[*itr]/f.alphasum;
     }
     for(map<Pair,float>::iterator itr = f.cooceff.begin();itr!=f.cooceff.end();++itr){
       //std::cout << "cooc val" << itr->second << " i=" << itr->first.a<<","<<itr->first.b << std::endl; 
       //std::cout << "occurence" << gsave[i].count(itr->first.a)<<gsave[i].count(itr->first.b) << std::endl; 
       if(gsave[i].count(itr->first.a)==0) {continue;}
       if(gsave[i].count(itr->first.b)==0) {continue;}
       
       std::cout << "   alpha[i]=" << itr->second << " i=" << itr->first.a<<","<<itr->first.b << std::endl; 
       sum += 2*itr->second/f.alphasum;
     }
     
     //std::cout << "   cooc_sum=" << f.cooceff.size() << std::endl; 
     
     ofs << f.gdata[i].value<< "\t" << sum  << std::endl;
     if (sum >= 0 && f.gdata[i].value > 0 ){
       cor += 1;
       std::cout << "res + g=" << i+1 << " p=" << sum << " y=" << f.gdata[i].value << std::endl;
     }else if (sum < 0 && f.gdata[i].value < 0 ){
       cor += 1;
       std::cout << "res + g=" << i+1 << " p=" << sum << " y=" << f.gdata[i].value << std::endl;
     }else if (sum >= 0 && f.gdata[i].value < 0 ){
       err += 1;
       std::cout << "res - g=" << i+1 << " p=" << sum << " y=" << f.gdata[i].value << std::endl;
     }else{
       err += 1;
       std::cout << "res - g=" << i+1 << " p=" << sum << " y=" << f.gdata[i].value << std::endl;
     }
   }
   
   std::cout << "---" << std::endl;
   std::cout << "correct:" << cor << " wrong:" << err << "  " << cor/(cor+err) << std::endl;
   
  return 0;
}

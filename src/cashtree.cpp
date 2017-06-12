#include "gspan.h"
#include <sstream>
#include <set>
#include <strstream>

void Ctree::print(){
  std::cout << g2tracers.size() << ":" << pat.rebuild() <<std::endl;
  if(children.size()!=0){
    for(list<Ctree*>::iterator it=children.begin();it!=children.end();++it){
      std::cout << "-";
      (*it)->print();
    }
  }
}
void CRoot::print(){
    for(list<Ctree*>::iterator it=one_edge_graphs.begin();it!=one_edge_graphs.end();++it){
      (*it)->print();
    }
}
void Tdelete(Ctree* tree){
  if(tree->children.size()!=0){
    for(list<Ctree*>::iterator it=tree->children.begin();it!=tree->children.end();++it){
      Tdelete(*it);
    }
  }
  delete tree;
}
void Rdelete(CRoot* croot){
    if(croot->one_edge_graphs.size()!=0){
    for(list<Ctree*>::iterator it=croot->one_edge_graphs.begin();it!=croot->one_edge_graphs.end();++it){
      Tdelete(*it);
    }
  }
  delete croot;
}

void Gspan::Crun(){
  if(first_flag==true){
    first_tree_make();
    first_flag=false;
    std::cout <<"CashTree Node size : " << TNnum <<std::endl;
    return;
  }
  can_grow.resize(0);
  CashTree_search();
  can_grow_search();
  std::cout <<"CashTree Node size : " << TNnum <<std::endl;

  //cooc search...
  if(!cooc) return;
  std::cout <<"Cooc searching ..."  <<std::endl;
  cooc_is_opt = false;
  opt_pat_cooc.gain = 0.0;
  opt_pat_cooc.sumofsize = 0;
  CoocSearch();
}

void Gspan::first_tree_make(){
  /***   init CRoot         ***/
  croot = new CRoot;
  TNnum = 1;
  /****  construct CRoot   ****/
  map<Triplet,GraphToTracers> heap;
  for(unsigned int gid = 0; gid < gdata.size(); ++gid){
    EdgeTracer cursor;
    Triplet wild_edge;
    Graph& g = gdata[gid];

    for(unsigned int v=0; v<g.size(); ++v){
      for(vector<Edge>::iterator e = g[v].begin(); e != g[v].end(); ++e){
	if (e->labels.x <= e->labels.z){
	  cursor.set(v,e->to,e->id,0);
	  heap[e->labels][gid].push_back(cursor);
	  if(wildcard_r>0){
	    wild_edge = e->labels;
	    wild_edge.z =999;
	    heap[wild_edge][gid].push_back(cursor);
	    wild_edge = e->labels.reverse();
	    wild_edge.z =999;
	    cursor.set(e->to,v,e->id,0);
	    heap[wild_edge][gid].push_back(cursor);
	  }
	}
      }
    }
  }
  pattern.resize(1);
  int lwild = wildcard_r;
  for(map<Triplet,GraphToTracers>::iterator it = heap.begin(); it != heap.end(); ++it){
    if(support(it->second) < minsup){ continue;}
    pattern[0].labels = it->first;
    pattern[0].time.set(0,1);
    
    Ctree *node = new Ctree;
    TNnum++;
    node->wildcard_res = lwild;
    nDFSCode d;
    d.set(pattern[0],NULL);
    node->pat = d;
    if(pattern[0].labels.z == 999){
      node->wildcard_res = lwild - 1;
    }
    GraphToTracers& g2 = node->g2tracers;
    g2 = it->second;
    node->children.resize(0);
    croot->one_edge_graphs.push_back(node);
    edge_grow(*node);
    pattern.resize(1);
  }
}

void Gspan::edge_grow(Ctree& nnode){
  
  if(can_prune(nnode)) { return;}

  if(pattern.size() == maxpat){return;}
  
  PairSorter b_heap;
  map<int,PairSorter,greater<int> > f_heap;
  wildcard_r = nnode.wildcard_res;
  int maxtoc = scan_gspan(nnode.g2tracers,b_heap,f_heap);
  
  // projecting...
  DFSCode  dcode;
  for(PairSorter::iterator it = b_heap.begin(); it != b_heap.end(); ++it){	
    dcode.labels = Triplet(-1,it->first.b,-1);
    dcode.time.set(maxtoc, it->first.a);
    pattern.push_back(dcode);

    if(support(it->second) < minsup){pattern.pop_back(); continue;}
    if(!is_min()){pattern.pop_back(); continue;}

    //new Node projecting ...
    Ctree *node = new Ctree;
    TNnum++;
    node->wildcard_res = nnode.wildcard_res;
    nDFSCode d;
    d.set(dcode,&nnode.pat);
    node->pat = d;
    GraphToTracers& g2 = node->g2tracers;
    g2 = it->second;
    nnode.children.push_back(node);
    
    edge_grow(*node);
    pattern.pop_back();
  }
	
  for(map<int,PairSorter,greater<int> >::iterator it = f_heap.begin(); it != f_heap.end(); ++it){	
    for(PairSorter::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){		
      dcode.labels = Triplet(-1,it2->first.a,it2->first.b);
      dcode.time.set(it->first,maxtoc+1);
      pattern.push_back(dcode);

      if(support(it2->second) < minsup){pattern.pop_back(); continue;}
      if(!is_min()){pattern.pop_back(); continue;}
      
      //new Node projecting ...
      Ctree *node = new Ctree;
      TNnum++;
      node->wildcard_res = nnode.wildcard_res;
      nDFSCode d;
      d.set(dcode,&nnode.pat);
      node->pat = d;
      
      if(pattern[pattern.size()-1].labels.z == 999){
	node->wildcard_res = nnode.wildcard_res - 1;
      }
      GraphToTracers& g2 = node->g2tracers;
      g2 = it2->second;
      nnode.children.push_back(node);
      
      edge_grow(*node);
      pattern.pop_back();
    }
  }
}

bool Gspan::can_prune(Ctree& node){

  double gain=0.0;
  double upos=0.0;
  double uneg=0.0;

  gain=-wbias;
  upos=-wbias;
  uneg=wbias;

  for(GraphToTracers::iterator it=node.g2tracers.begin();it!=node.g2tracers.end();++it){
    int gid = it->first;
    gain += 2 * corlab[gid] * weight[gid];
    if(corlab[gid]>0){
      upos += 2 * weight[gid];
	}else{
      uneg += 2 * weight[gid];
    }
  }
  node.max_gain = std::max(upos,uneg);
  if(fabs(opt_pat.gain) > node.max_gain ){ return true;}

  double gain_abs = fabs(gain);
  if(gain_abs > fabs(opt_pat.gain) || (fabs(gain_abs - fabs(opt_pat.gain))<1e-10 && pattern.size() < opt_pat.size)){
    opt_pat.gain = gain;
    opt_pat.optimalplace = &node;
    opt_pat.size = pattern.size();
  }
  return false;
}

void Gspan::CashTree_search(){
  //std::cout << "cash tree search..." << std::endl;
  pattern.resize(0);
  for(list<Ctree*>::iterator it = croot->one_edge_graphs.begin();it != croot->one_edge_graphs.end();++it){
    pattern.push_back((*it)->pat.dcode);
    node_search(*(*it));
    pattern.pop_back();
  }
}

void Gspan::node_search(Ctree& node){
  
  if(can_prune(node)) { return;}
  if(pattern.size() >= maxpat){ return;}
  if(node.children.size() == 0){
    can_grow.push_back(&node);
  }
  
  for(list<Ctree*>::iterator it = node.children.begin();it != node.children.end();++it){
    pattern.push_back((*it)->pat.dcode);
    node_search(**it);
    pattern.pop_back();
  }
}

void Gspan::can_grow_search(){
  for(list<Ctree*>::iterator cit = can_grow.begin();cit!=can_grow.end();++cit){
    pattern = (*cit)->pat.rebuild();
    //std::cout << pattern << std::endl;
    //edge_grow(**it);

    if(fabs(opt_pat.gain) > (*cit)->max_gain ){ continue;}
    
    PairSorter b_heap;
    map<int,PairSorter,greater<int> > f_heap;
    wildcard_r = (*cit)->wildcard_res;
    int maxtoc = scan_gspan((*cit)->g2tracers,b_heap,f_heap);
  
    // projecting...
    DFSCode  dcode;
    for(PairSorter::iterator it = b_heap.begin(); it != b_heap.end(); ++it){	
      dcode.labels = Triplet(-1,it->first.b,-1);
      dcode.time.set(maxtoc, it->first.a);
      pattern.push_back(dcode);

      if(support(it->second) < minsup){pattern.pop_back(); continue;}
      if(!is_min()){pattern.pop_back(); continue;}

      //new Node projecting ...
      Ctree *node = new Ctree;
      TNnum++;
      node->wildcard_res = (*cit)->wildcard_res;
      nDFSCode d;
      d.set(dcode,&(*cit)->pat);
      node->pat = d;
      GraphToTracers& g2 = node->g2tracers;
      g2 = it->second;
      (*cit)->children.push_back(node);
      
      edge_grow(*node);
      pattern.pop_back();
  }
	
  for(map<int,PairSorter,greater<int> >::iterator it = f_heap.begin(); it != f_heap.end(); ++it){	
    for(PairSorter::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){		
      dcode.labels = Triplet(-1,it2->first.a,it2->first.b);
      dcode.time.set(it->first,maxtoc+1);
      pattern.push_back(dcode);

      if(support(it2->second) < minsup){pattern.pop_back(); continue;}
      if(!is_min()){pattern.pop_back(); continue;}
      
      //new Node projecting ...
      Ctree *node = new Ctree;
      TNnum++;
      node->wildcard_res = (*cit)->wildcard_res;
      nDFSCode d;
      d.set(dcode,&(*cit)->pat);
      node->pat = d;
      
      if(pattern[pattern.size()-1].labels.z == 999){
	node->wildcard_res = (*cit)->wildcard_res - 1;
      }
      GraphToTracers& g2 = node->g2tracers;
      g2 = it2->second;
      (*cit)->children.push_back(node);
      
      edge_grow(*node);
      pattern.pop_back();
    }
  }


    
  }
}

void update(DPat& op){
  op.locsup.clear();
  for(GraphToTracers::iterator it=op.optimalplace->g2tracers.begin();it!=op.optimalplace->g2tracers.end();++it){
    op.locsup.push_back(it->first);
  }
    std::ostrstream ostrs;
    ostrs << op.optimalplace->pat.rebuild();
    ostrs << std::ends;
    op.dfscode = ostrs.str();

}
void update(CDPat& op){
  
  std::ostrstream ostrs1;
  ostrs1 << op.optimalplaces[0]->pat.rebuild();
  ostrs1 << std::ends;
  op.dfscodes[0] = ostrs1.str();
  std::ostrstream ostrs2;
  ostrs2 << op.optimalplaces[1]->pat.rebuild();
  ostrs2 << std::ends;
  op.dfscodes[1] = ostrs2.str();

}


void Gspan::CoocSearch(){
  
  //edge_grow_before_cooc();
  count = 0;
  pattern.resize(0);
  cpattern.resize(0);
  opt_gain = opt_pat.gain;
  for(list<Ctree*>::iterator it = croot->one_edge_graphs.begin();it != croot->one_edge_graphs.end();++it){
    pattern.push_back((*it)->pat.dcode);
    edge_grow_before_cooc(*(*it));
    pattern.pop_back();
  }
  std::cout << "cooc check count : "<<count << std::endl;
}

void Gspan::edge_grow_before_cooc(Ctree& left){
  if(left.max_gain < opt_gain) return;
  //std::cout << left.pat.rebuild() << std::endl;

  for(list<Ctree*>::iterator it = croot->one_edge_graphs.begin();it != croot->one_edge_graphs.end();++it){
    if((*it)->pat.dcode == pattern[0]){
      if(pattern.size() == cpattern.size() + 1) break;
      cpattern.push_back((*it)->pat.dcode);
      edge_grow_cooc(left,*(*it));
      cpattern.pop_back();
      break;
    }
    cpattern.push_back((*it)->pat.dcode);
    edge_grow_cooc_(left,*(*it));
    cpattern.pop_back();
  }
  
  for(list<Ctree*>::iterator it = left.children.begin();it != left.children.end();++it){
    pattern.push_back((*it)->pat.dcode);
    edge_grow_before_cooc(**it);
    pattern.pop_back();
  }
}

void Gspan::edge_grow_cooc(Ctree& left,Ctree& right){
  if(right.max_gain < opt_gain) return;
  if(pattern.size() + cpattern.size() > cmaxpat) return;
  //count++;
  for(list<Ctree*>::iterator it = right.children.begin();it != right.children.end();++it){
    if((*it)->pat.dcode == pattern[cpattern.size()]){
      if(pattern.size() == cpattern.size() + 1) break;
      cpattern.push_back((*it)->pat.dcode);
      edge_grow_cooc(left,**it);
      cpattern.pop_back();
      break;
    }
    cpattern.push_back((*it)->pat.dcode);
    edge_grow_cooc_(left,**it);
    cpattern.pop_back();
  }
}
void Gspan::edge_grow_cooc_(Ctree& left,Ctree& right){
  if(right.max_gain < opt_gain) return;
  if(pattern.size() + cpattern.size() > cmaxpat) return;
  count++;//cooc check count
  if(can_prune(left,right)) return;
  for(list<Ctree*>::iterator it = right.children.begin();it != right.children.end();++it){
    cpattern.push_back((*it)->pat.dcode);
    edge_grow_cooc(left,**it);
    cpattern.pop_back();
  }
}

bool Gspan::can_prune(Ctree& left,Ctree& right){

  double gain=0.0;
  double upos=0.0;
  double uneg=0.0;

  gain=-wbias;
  upos=-wbias;
  uneg=wbias;
  
  vector<int> lsup;
  GraphToTracers::iterator rit = right.g2tracers.begin();
  
  for(GraphToTracers::iterator it=left.g2tracers.begin();it!=left.g2tracers.end();++it){
    int gid = it->first;
    while(rit != right.g2tracers.end() && gid >rit->first){++rit;}
    if(gid != rit->first){continue;}
    lsup.push_back(gid);
    gain += 2 * corlab[gid] * weight[gid];
    if(corlab[gid]>0){
      upos += 2 * weight[gid];
	}else{
      uneg += 2 * weight[gid];
    }
  }
  double mgain = std::max(upos,uneg);
  if(fabs(opt_pat.gain) > mgain ){ return true;}
  if(fabs(opt_pat_cooc.gain) > mgain ){ return true;}

  double gain_abs = fabs(gain);
  if(!cooc_is_opt){//max pattern is nomal
    double d = gain_abs - fabs(opt_pat.gain);
    if(d < -1e-10) return false;
    if((fabs(d)<1e-10 && pattern.size()+cpattern.size() > opt_pat.size)) return false;
    opt_pat_cooc.gain = gain;
    opt_pat_cooc.optimalplaces[0]=&left;
    opt_pat_cooc.optimalplaces[1]=&right;
    opt_pat_cooc.sumofsize = pattern.size() + cpattern.size();
    cooc_is_opt = true;
    opt_pat_cooc.locsup = lsup;
    
  }else{//max pattern is cooc    
    double d = gain_abs - fabs(opt_pat_cooc.gain);
    if(d < -1e-10) return false;
    if((fabs(d)<1e-10 && pattern.size()+cpattern.size() > opt_pat_cooc.sumofsize)) return false;
    opt_pat_cooc.gain = gain;
    opt_pat_cooc.optimalplaces[0]=&left;
    opt_pat_cooc.optimalplaces[1]=&right;
    opt_pat_cooc.sumofsize = pattern.size() + cpattern.size();
    opt_pat_cooc.locsup = lsup;
  }

  return false;
}

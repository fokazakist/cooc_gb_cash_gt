#include <iostream>
#include <sstream>
#include <fstream>
#include "getopt.h"
#include "gspan.h"
#include <strstream>
extern "C" {
#include "../glpk-4.8/include/glpk.h"
}
#define OPT " [-m minsup] [-x maxpat] [-w wildcard] [-n v] [-e conv_epsilon] [-o] graph-file"
#define ROW(i) ((i)+1)
#define COL(j) ((j)+1)

int main(int argc, char **argv) {
  unsigned int minsup = 1;
  unsigned int minp = 0;
  unsigned int maxpat = 10;
  bool out_instances = false;
  bool percent=false;
  int opt;
  int wildcard_num = 0;
  unsigned int maxitr = 500000;
  double nu = 0.4; 
  double conv_epsilon = 1e-2;
  bool cooc = false;
  unsigned int cmaxpat = 11;

  //clock_t allstart, allend;

  while ((opt = getopt(argc, argv, "m:p:w:e:n:x:s:io")) != -1) {
    switch (opt) {
    case 'm':
      minsup = atoi (optarg);
      break;
    case 'p':
      minp = atoi (optarg);
      percent = true;
      break;
    case 'n':
      nu = atof(optarg);
      break;
    case 'e':
      conv_epsilon = atof(optarg);
      break;
    case 'x':
      maxpat = atoi (optarg);
      break;
    case 's':
      cmaxpat = atoi (optarg);
      break;
    case 'w':
      wildcard_num = atoi (optarg);
      break;
    case 'i':
      out_instances = true;
      break;
    case 'o':
      cooc = true;
      break;
    default:
      std::cerr << "Usage: "<< argv[0] << OPT<< std::endl;
      return -1;
    }
  }
	
  if(argc-optind != 1){
    std::cerr << "Usage: "<< argv[0] << OPT<< std::endl;
    return -1;
  }
	
  std::ifstream graph_file(argv[optind++]);
  if(graph_file.fail()){
    std::cerr << "File not found: " << argv[optind-1] << std::endl;
    return -1;
  }
  Gspan gspan;
  gspan.wildcard_r = wildcard_num;
  gspan.maxpat = maxpat;
  gspan.cmaxpat = cmaxpat;
  gspan.out_instances = out_instances;
  gspan.max_itr = maxitr;
  
  gspan.set_data(graph_file);
  gspan.minsup = minsup;
  gspan.nu = nu;
  gspan.conv_epsilon = conv_epsilon;
  gspan.cooc = cooc;
  if(percent==true){
    gspan.minsup = gspan.gdata.size() * minp / 100;
  }
  gspan.lpboost();
  
  std::cout << "Given options::" << "maxpat: " << maxpat << " minsup: " << gspan.minsup << " nu: " << nu << " conv_epsilon: " << conv_epsilon <<" maxitr: " << maxitr << std::endl;
  Rdelete(gspan.croot);
  return 0;
}

struct Hypothesis {
public:
  vector<double> weight;
  vector<std::string> dfs_vector;
  vector<int> flag;
  map<int,vector<int> > tmp;
  explicit Hypothesis(){
    weight.resize(0);
    dfs_vector.resize(0);
    flag.resize(0);
  }
};

void Gspan::lpboost(){
  // "out" is name of output model file
  const char *out = "model";
  
  //initialize
  unsigned int gnum = gdata.size(); 
  weight.resize(gnum);
  std::fill(weight.begin(),weight.end(),1.0);
  corlab.resize(gnum);
  for(unsigned int gid=0;gid<gnum;++gid){
    corlab[gid]=gdata[gid].class_label;
  }
  wbias=0.0;
  Hypothesis model;
  first_flag=true;

  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  std::cout.precision(8);
  
  //Initialize GLPK
  int* index = new int[gnum+2]; double* value = new double[gnum+2];
  LPX* lp = lpx_create_prob();
		       
  lpx_add_cols(lp, gnum+1); // set u_1,...u_l, beta
  for (unsigned int i = 0; i < gnum; ++i){
    lpx_set_col_bnds(lp, COL(i), LPX_DB, 0.0, 1/(nu*gnum));
    lpx_set_obj_coef(lp, COL(i), 0); // u ... lambda
  }
  lpx_set_col_bnds(lp, COL(gnum), LPX_FR, 0.0, 0.0);
  lpx_set_obj_coef(lp, COL(gnum), 1); // beta ... gamma
  lpx_set_obj_dir(lp, LPX_MIN); //optimization direction: min objective
		       
  lpx_add_rows(lp,1); // Add one row constraint s.t. sum_u == 1
  for (unsigned int i = 0; i < gnum; ++i){
    index[i+1] = COL(i);
    value[i+1] = 1;
  }
  lpx_set_mat_row(lp, ROW(0), gnum, index, value);
  lpx_set_row_bnds(lp, ROW(0), LPX_FX, 1, 1);
		       
  double beta = 0.0;
  double margin = 0.0;
  unsigned int litr = 0; // the number of iterations
  
  //main loop
  for(unsigned int itr=0;itr < max_itr;++itr){
    std::cout <<"itrator : "<<itr+1<<std::endl;
    
    opt_pat.gain=0.0;//gain init
    opt_pat.size=0;
    //std::cout <<"cooc gain : "<<opt_pat_cooc.gain<<std::endl;
    pattern.resize(0);
    
    Crun(); //search optimal pattern
    
    std::vector <int>     result (gnum);
    int _y;
    vector<int> locvec;  //= opt_pat.locsup;
    std::string dfscode; //= opt_pat.dfscode;
    // _y is omega
    if(cooc_is_opt){
      _y = opt_pat_cooc.gain > 0 ? +1 :-1;
      update(opt_pat_cooc); //rewrite opt_pat_cooc
      locvec = opt_pat_cooc.locsup;
      dfscode = opt_pat_cooc.dfscodes[0]+"\t"+opt_pat_cooc.dfscodes[1];
    }else{
      _y = opt_pat.gain > 0 ? +1 :-1;
      update(opt_pat); //rewrite opt_pat
      locvec = opt_pat.locsup;
      dfscode = opt_pat.dfscode;
    }
    
    //print new hypo
    if(!cooc_is_opt){
      std::cout << "not cooc :::";
      std::cout << opt_pat.gain << " : " << opt_pat.locsup.size();
      std::cout << " *  " << opt_pat.dfscode << std::endl;
    }else{
      std::cout << "    cooc :::";
      std::cout << opt_pat_cooc.gain << " : " << locvec.size();
      std::cout << " *  " << dfscode << std::endl;
    }
   
    model.flag.resize(itr+1);
    model.flag[itr] = _y;
    model.tmp[itr]  = locvec;

    std::fill(result.begin (), result.end(), -_y);
    for (unsigned int i = 0; i < locvec.size(); ++i){ result[locvec[i]] = _y;}
    double uyh = 0;
    for (unsigned int i = 0; i < gnum;  ++i) { // summarizing hypotheses
      uyh += weight[i]*corlab[i]*result[i];
    }
      
    std::cout << "Stopping criterion: " << uyh << "<=?";
    std::cout << beta << " + " << conv_epsilon << std::endl;

    if( (uyh <= beta + conv_epsilon ) ){
      std::cout << "*********************************" << std::endl;
      std::cout << "Convergence ! at iteration: " << itr+1 << std::endl;
      std::cout << "*********************************" << std::endl;
      if(cooc){
	litr = itr;
	break;
      }
      cooc = true;
    }
      
    lpx_add_rows(lp,1); // Add one row constraint s.t. sum( uyh - beta ) <= 0
    for (unsigned int i = 0; i < gnum; ++i){
      index[i+1] = COL(i);
      value[i+1] = result[i] * corlab[i];
    }
    index[gnum+1] = COL(gnum);
    value[gnum+1] = -1;
    lpx_set_mat_row(lp, ROW(itr+1), gnum+1, index, value);
    lpx_set_row_bnds(lp, ROW(itr+1), LPX_UP, 0.0, 0.0);

    model.weight.push_back(0);
    model.dfs_vector.push_back(dfscode);
      
    lpx_simplex(lp); 
    beta = lpx_get_obj_val(lp);
    for (unsigned int i = 0; i < gnum; ++i){
      double new_weight;
      new_weight = lpx_get_col_prim(lp, COL(i));
      if(new_weight < 0) new_weight = 0; // weight > 0
      weight[i] = new_weight;
    }
    margin = lpx_get_row_dual(lp, ROW(0));
    double margin_error = 0.0;
    for (unsigned int i = 0; i < gnum;  ++i) { // summarizing hypotheses
      if (corlab[i]*result[i] < margin){
	++margin_error;
      }
    }
    margin_error /= gnum;

    //next rule is estimated
    wbias = 0.0;
    for (unsigned int i = 0; i < gnum; ++i){
      wbias += corlab[i] * weight[i];
    }

    std::cout << "After iteration " << itr+1 << std::endl;
    std::cout << "Margin: " << margin << std::endl;
    std::cout << "Margin Error: " << margin_error << std::endl;
  }

  // output model
  std::ofstream os (out);
  if (! os) {
    std::cerr << "FATAL: Cannot open output file: " << out << std::endl;
    return;
  }
  os.setf(std::ios::fixed,std::ios::floatfield);
  os.precision(12);
  
  std::vector<float> pred(gnum);
  std::fill (pred.begin (), pred.end(), 0.0);
  for (unsigned int r = 0; r < litr; ++r){
    model.weight[r] = - lpx_get_row_dual(lp, ROW(r+1));
    if(model.weight[r] < 0) model.weight[r] = 0; // alpha > 0
    os << model.flag[r] * model.weight[r] << "\t" << model.dfs_vector[r] << std::endl;
  }
  
  delete [] index; delete [] value;
  lpx_delete_prob(lp);
}

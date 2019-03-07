#include "single_module.h"
#include "setup_config.h"
#include "TBReader.h"
#include "fitter.h"
#include <fstream>
#include <iostream>
#include "TCanvas.h"

// General Setup
string main_outpath = "./";
string main_default_rootname = "TPro.root";
string main_datainput = "./data_input.txt"; 
string Module_configfile = "./configs/all_config.csv";
int main_config = 3; // Run 179 ~ 722(config1), Run 751 ~ 1064(config2),
                // Run 1079~1167(config3)
// Using the first configuration so all boards are included
// 6-Nov-2018 copy from
// https://docs.google.com/spreadsheets/d/1KcvFr3JG69plQeVy4AA8nI63ZuTMc6hGDWGwo-F5xPE/edit#gid=337823033


// Usage
void main_Print_help();
void main_make_TProfile(string TProfile_name);
void main_make_Inj_TProfile(string TProfile_name);
void main_make_module_ntuple();
void main_fitter(string TProfile_name);
bool main_check_fname(string fname);
string arg_string;
vector<string> all_args;
bool DBG = false;

int main(int argc, char* argv[]){
  TApplication *app = new TApplication("app",0,0);
  app->Init();
  string TProfile_name = main_default_rootname;
  for(int i = 0 ; i < argc ; ++i){
    arg_string = argv[i];
    all_args.push_back(arg_string); }
  if(argc == 1){
    main_Print_help();
    return 1; }
  else if(argc == 2){
    if(all_args[1] == "-h" || all_args[1] == "-H" ){ main_Print_help(); }
    else if(all_args[1] == "-t" || all_args[1] == "-T" ){
      main_make_TProfile(TProfile_name); }
    else if(all_args[1] == "-i" || all_args[1] == "-I" ){
      main_make_Inj_TProfile(TProfile_name); }
    else if(all_args[1] == "-n" || all_args[1] == "-N" ){
      main_make_module_ntuple(); }
    else if(all_args[1] == "-f" || all_args[1] == "-F" ){
      main_fitter(TProfile_name); }
    else{
      std::cout << "Unknown option... print usage" << std::endl;
      main_Print_help(); }
    return 1;
  }
  else if(argc == 3){
    bool check = main_check_fname(TProfile_name);
    if(!check){ TProfile_name = main_default_rootname; }
    
    if(all_args[1] == "-t" || all_args[1] == "-T" ){
      TProfile_name = all_args[2]; 
      main_make_TProfile(TProfile_name); }
    else if(all_args[1] == "-i" || all_args[1] == "-I" ){
      TProfile_name = all_args[2];
      main_make_Inj_TProfile(TProfile_name); }
    else if(all_args[1] == "-f" || all_args[1] == "-F" ){
      TProfile_name = all_args[2];
      main_fitter(TProfile_name);}
    else{
      std::cout << "Unknown option... print usage" << std::endl;
      main_Print_help();
    }
  }
  else{ cout << "unexpected number of option! QUIT!" << endl; }
}

void main_Print_help(){
  std::cout
  << "Usage: \n" << "(1) " << all_args[0] << " -h  : Print this message.\n"
  << "(2) "<< all_args[0] << " -t filename.root  : "
  << "Create root with TProfile in HGLG/LGTOT for each Board, "
  << "default filename = outpath/TPro.root.\n"
  << "(3) " << all_args[0] << " -n  : Create root ntuple with different"
  << "module.root default path is " << main_outpath
  << "Module_Ntuple/TB(or Inj)/module***.root \n"
  << "(4) " << all_args[0] << " -f filename.root  : Use filename.root as"
  << " input file for fitting( Should be an output file of "
  << all_args[0] << " -t)\n"
  << "(5) " << all_args[0] << " -i : Create Injection run TProfile\n\n"
  << "For Usage (1) and (2), One should write the testbeam ntuple "
  << "into " << main_datainput << std::endl;
};
bool main_check_fname(string fname){
  if( fname.find(".root") ){ return true; }
  else{
    std::cout << "You better choose a ***.root filename."
	      << "Set input/output name to default.(" << main_default_rootname
	      << ")" << std::endl;
    return false; }
}

void main_make_TProfile(string TProfile_name){
    
  // Initialize output directory
  setup_config *SC = new setup_config;
  SC->dirpath = main_outpath;
  SC->Make_dir();
  SC->Read_Module_List(Module_configfile,main_config); // Set ModuleID List && map
  
  int TB_member     = 0;
  
  string inputfile = main_datainput;
  ifstream infile(inputfile.c_str());

  TProfile_name = string( main_outpath + string("Module_TProfile/") + TProfile_name );
  cout << "Output file with be " << TProfile_name << endl;
  if(DBG){
    cout << "Press any key to continue...\n\n" << endl;
    getchar();}
  MakePlots *M = new MakePlots(SC);
  bool turefile = M->Init_TFile(TProfile_name);
  
  if(!turefile){ return; }

  M->root_logon();
  string filename;
  while(true){
    
    infile >> filename;
    if(infile.eof()) {
      M-> Write_TProfile();
      break;}
    if( filename.length() > 2){
      cout << "input file: " << filename << endl;

      TFile f( filename.c_str() );
      //check if root directories exist...
      bool TB_ntuple,trackimpactntupler;
      TB_ntuple   = f.GetListOfKeys()->Contains("rechitntupler");
      trackimpactntupler = f.GetListOfKeys()->Contains("trackimpactntupler");
      if(TB_ntuple){
	TB_member++;
	// if(TB_member != 0){
	//   cout << "DO NOT merge both TB and Injection runs in "
	//        << main_datainput << "!!!!\nBREAK!\n" << endl;
	//   break; }

	TChain *chain  = new TChain("rechitntupler/hits");
	chain ->Add(filename.c_str());
	TChain *chain2 = new TChain("trackimpactntupler/impactPoints");
	if(trackimpactntupler)
	  chain2->Add(filename.c_str());
	// ^^^^^^^^ Will simply not initialize DWC tree if it doesn't exist
	// Since we are not doing any event selection
	TBReader TBReader(chain,chain2,filename);
	if(!TBReader.Check_Config(main_config)){ continue ;}
	TBReader.dirpath = main_outpath;
	TBReader.TProfile_Maker(SC,M);
	if(TB_member %20 == 0){
	  M-> Write_TProfile();}
	delete chain;
	delete chain2;
      }
      else{
	cout << filename.c_str() << " contains unknown tree to me ..." << endl;
      }
      f.Close();
    }
    else{
      cout << "file " << filename << " is not available, please check "
	   << inputfile << endl;}
  }
  infile.close();
}

void main_make_Inj_TProfile(string TProfile_name){
  string inputfile = main_datainput;
  ifstream infile(inputfile.c_str());
  
  TProfile_name = string( main_outpath + string("Module_TProfile/") + TProfile_name );
  cout << "Output file with be " << TProfile_name << endl;
  if(DBG){
    cout << "Press any key to continue...\n\n" << endl;
    getchar();}
  
  string filename;

  while(true){    
    infile >> filename;
    if(infile.eof()) {
      break;}
    if( filename.length() > 2){
      cout << "input file: " << filename << endl;
    
      TFile f( filename.c_str() );
      //check if root directories exist...
      bool pulseshape_tree;
      pulseshape_tree = f.GetListOfKeys()->Contains("pulseshapeplotter");
      if(pulseshape_tree){
	TChain *chain_single  = new TChain("pulseshapeplotter/tree");
	chain_single->Add(filename.c_str());
	single_module S(chain_single,filename,TProfile_name);
	S.Loop();
	delete chain_single;
      }
      else{
	cout << filename.c_str() << " contains unknown tree to me ..." << endl;}
      f.Close();
    }
    else{
      cout << "file " << filename << " is not available, please check "
	   << inputfile << endl;
    }
  }
  infile.close();
}
void main_make_module_ntuple(){
  std::cout << "Under construction now QQ" << std::endl;
}
void main_fitter(string TProfile_name){


  cout << "Fitting file with be " << TProfile_name << endl;
  if(DBG){
    cout << "Press any key to continue...\n\n" << endl;
    getchar();}
  
  // Initialize output directory
  setup_config *SC  = new setup_config;
  SC->Read_Module_List(Module_configfile,main_config); // Set ModuleID List && map, config == 1 for fitter
  
  MakePlots *M = new MakePlots(SC);
  bool turefile = M->Init_TFile(TProfile_name);
  if(!turefile){ return; }

  M->root_logon();
  fitter f(SC,M,TProfile_name);  
  //f.fit_LGTOT();
  f.fit_output();
};

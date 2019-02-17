#ifndef MakePlots_h
#define MakePlots_h

#include "setup_config.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
class MakePlots{
 public:
  MakePlots(setup_config *SC);
  ~MakePlots();

  //Histogram members
  TH2D *h_HG_LG [MAXBOARDS][MAXSKI][MAXCH];
  //TH2D *h_LG_TOT[MAXBOARDS][MAXSKI][MAXCH];
  
  //TProfile members
  TProfile *HG_LG [MAXBOARDS][MAXSKI][MAXCH];
  TProfile *LG_TOT[MAXBOARDS][MAXSKI][MAXCH];
 
  // Function
  bool Init_TFile(string TPro_outputname);
  void Write_TProfile();  // Takes time, fill every few runs
  bool Check_Run(int RunN); //For TB runs
  bool Check_Name(string fname); // For injection runs
  bool Check_Fit_Type(); // false = injection run, true = TB  

  // Drawing Function
  void root_logon();
  
  // member
  bool file_exist;
  int  history_Run;
  string TPro_output;
  string m_filename;

  
 private:
  // Function
  void Init_Pointers();

  TFile *TPro_root;
  TTree *TPro_history;
  TTree *TPro_fname;
  
  setup_config *mysetup;

};

#endif

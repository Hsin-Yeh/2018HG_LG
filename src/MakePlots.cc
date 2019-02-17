#include "MakePlots.h"
#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"

MakePlots::MakePlots(setup_config *SC){
  cout << "Constructor of MakePlots ..." << endl;
  mysetup = SC;
}
MakePlots::~MakePlots(){
  cout << "Destructor of MakePlots ..." << endl;
}

bool MakePlots::Init_TFile(string TPro_outputname){
  Init_Pointers();
  ifstream f_check(TPro_outputname.c_str());
  file_exist = f_check.good();
  TPro_output = TPro_outputname;
  char p_name[200],HGLG_name[200],LGTOT_name[200];
  if(file_exist){
    TPro_root = new TFile(TPro_outputname.c_str(),"update");
    // Check TTree exist
    TPro_history = (TTree*)TPro_root->Get("history");
    if(TPro_history == NULL){ 
      cout << "\nInput file has no object named " << "history! " << endl;
      cout << "Please choose another root file or create one.\n\n" << endl;
      return false;}
    TPro_history -> SetBranchAddress("history_Run",&history_Run);

    // TPro_fname = (TTree*)TPro_root->Get("Filename");
    // if(TPro_fname == NULL){ 
    //   cout << "\nInput file has no object named " << "Filename! " << endl;
    //   cout << "Please choose another root file or create one.\n\n" << endl;
    //   return false;}
    // TPro_fname -> SetBranchAddress("Filename_Inj",&m_filename);

    
    // Get the TProfiles
    for(int BD = 0 ; BD < MAXBOARDS ; ++BD){
      int moduleID = mysetup->Module_List[BD];
      if( moduleID == 0 ) { continue; }
      for(int chip = 0 ; chip < MAXSKI ; ++chip){
	for(int ch = 0 ; ch < MAXCH ; ++ch){
	  sprintf(HGLG_name,"Module%d/HGLG_chip%d_ch%d",moduleID,chip,ch*2);
	  sprintf(LGTOT_name,"Module%d/LGTOT_chip%d_ch%d",moduleID,chip,ch*2);
	  HG_LG[BD][chip][ch] = (TProfile*)TPro_root->Get(HGLG_name);
	  LG_TOT[BD][chip][ch] = (TProfile*)TPro_root->Get(LGTOT_name);
	  // Check TProfile exist
	  if(TPro_history ->GetEntries() != 0){
	    if(HG_LG[BD][chip][ch] == NULL || LG_TOT[BD][chip][ch] == NULL){
	      cout << "\nFile " << TPro_outputname << " has no object named "
		   << HGLG_name << " or " << LGTOT_name << "\n"
		   << "Please choose another root file or create one.\n"
		   << endl;
	      return false; }
	  }
	}
      }
    }
  }
  else{
    TPro_root = new TFile(TPro_outputname.c_str(),"recreate");
    TPro_history = new TTree("history","history");
    TPro_history-> Branch("history_Run",&history_Run);
    // TPro_fname   = new TTree("Filename","Filename");
    // TPro_history-> Branch("Filename_Inj",&m_filename);

    
    int HGLGBIN  = 400;
    int LGTOTBIN = 200;
  
    for(int BD = 0 ; BD < MAXBOARDS ; ++BD){
      int moduleID = mysetup->Module_List[BD];
      if(moduleID == 0){ continue;}
      for(int chip = 0 ; chip < MAXSKI ; ++chip){
	for(int ch = 0 ; ch < MAXCH ; ++ch){
	  sprintf(p_name,"HG_LG_Module%d_chip%d_ch%d",moduleID,chip,ch*2);
	  HG_LG[BD][chip][ch] = new TProfile(p_name,"",HGLGBIN,0,800,0,4000);
	  sprintf(p_name,"Histo_HG_LG_Module%d_chip%d_ch%d",moduleID,chip,ch*2);
	  h_HG_LG[BD][chip][ch] = new TH2D(p_name,"",HGLGBIN,0,800,HGLGBIN,0,4000);
	  sprintf(p_name,"LG_TOT_Module%d_chip%d_ch%d",moduleID,chip,ch*2);
	  LG_TOT[BD][chip][ch] = new TProfile(p_name,"",LGTOTBIN,0,800,0,2000);
	  sprintf(p_name,"Histo_LG_TOT_Module%d_chip%d_ch%d",moduleID,chip,ch*2);
	  h_LG_TOT[BD][chip][ch] = new TH2D(p_name,"",LGTOTBIN,0,800,LGTOTBIN,0,2000);
	  HG_LG[BD][chip][ch]->SetMarkerStyle(22);
	  HG_LG[BD][chip][ch]->SetMarkerColor(chip+1);
	  HG_LG[BD][chip][ch]->SetMarkerSize(1);
	  HG_LG[BD][chip][ch]->SetYTitle("HG(ADC Counts)");
	  HG_LG[BD][chip][ch]->SetXTitle("LG(ADC Counts)");
	  
	  LG_TOT[BD][chip][ch]->SetMarkerStyle(21);
	  LG_TOT[BD][chip][ch]->SetMarkerColor(chip+1);
	  LG_TOT[BD][chip][ch]->SetMarkerSize(1);
	  LG_TOT[BD][chip][ch]->SetXTitle("TOT(ADC Counts)");
	  LG_TOT[BD][chip][ch]->SetYTitle("LG(ADC Counts)");
	}
      }
    }
  }
  return true;
}
void MakePlots::Init_Pointers(){
  for(int BD = 0 ; BD < MAXBOARDS ; ++BD){
    for(int chip = 0 ; chip < MAXSKI ; ++chip){
      for(int ch = 0 ; ch < MAXCH ; ++ch){
	HG_LG [BD][chip][ch] = 0;
        LG_TOT[BD][chip][ch] = 0;
      }
    }
  }
}
bool MakePlots::Check_Fit_Type(){
  if( TPro_history->GetEntries() == 0){ return false; }
  else{ return true; }
}
bool MakePlots::Check_Run(int RunN){

  bool doublefill = false;
  int history_run = TPro_history->GetEntries(); 
  for(int i = 0 ; i < history_run ; ++i){
    TPro_history->GetEntry(i);
    if(history_Run == RunN){ doublefill = true; }
  }
  if(!doublefill){
    history_Run = RunN;
    TPro_history->Fill();}
  
  return doublefill;
}
bool MakePlots::Check_Name(string fname){

  bool doublefill = false;
  
  int history_runs = TPro_fname->GetEntries(); 
  for(int i = 0 ; i < history_runs ; ++i){
    TPro_fname->GetEntry(i);
    if(fname == m_filename){ doublefill = true; }
  }
  if(!doublefill){
    m_filename = fname;
    TPro_fname->Fill();}
  
  return doublefill;
}

void MakePlots::Write_TProfile(){
  TPro_root->cd();
  cout << "Filling TProfiles ..." << endl;
  TCanvas* C1 = new TCanvas();
  TDirectory *dir;
  char title[50];
  
  for(int BD = 0 ; BD < MAXBOARDS ; BD++){
    sprintf(title,"Module%i",mysetup->Module_List[BD]);
    if(mysetup->Module_List[BD] == 0){ continue; } //If module doesn't exist
    if(!file_exist){
      dir = new TDirectory();
      dir = TPro_root->mkdir(title); }
    else{ dir = (TDirectory*)TPro_root->Get(title); }
    dir->cd();
    
    for(int chip = 0 ; chip < MAXSKI ; ++chip){
      for(int ch = 0 ; ch < MAXCH ; ++ch){
	sprintf(title,"HGLG_chip%i_ch%i",chip,ch*2);
	HG_LG[BD][chip][ch]->SetTitle(title);
	HG_LG[BD][chip][ch]->Write(title,TObject::kOverwrite);
	sprintf(title,"Histo_HGLG_chip%i_ch%i",chip,ch*2);
	h_HG_LG[BD][chip][ch]->SetTitle(title);
	h_HG_LG[BD][chip][ch]->Write(title,TObject::kOverwrite);
	//h_HG_LG[BD][chip][ch]->Draw("CANDLEX(30)");
	//h_HG_LG[BD][chip][ch]->Draw();
	//C1->Update();
	//gPad->WaitPrimitive();
	sprintf(title,"LGTOT_chip%i_ch%i",chip,ch*2);
	LG_TOT[BD][chip][ch]->SetTitle(title);
	LG_TOT[BD][chip][ch]->Write(title,TObject::kOverwrite);
	sprintf(title,"Histo_LGTOT_chip%i_ch%i",chip,ch*2);
	h_LG_TOT[BD][chip][ch]->SetTitle(title);
	h_LG_TOT[BD][chip][ch]->Write(title,TObject::kOverwrite);
      }
    }
  }
  
  TPro_root->cd();
  TPro_history->Write("history",TObject::kOverwrite);
  
  file_exist = true;
}
void MakePlots::root_logon(){

cout << endl << "Welcome to the rootlogon.C" << endl;
//
// based on a style file from BaBar
//

//..BABAR style from RooLogon.C in workdir
TStyle *atlasStyle= new TStyle("ATLAS","Atlas style");

// use plain black on white colors
 Int_t icol=0;
atlasStyle->SetFrameBorderMode(icol);
atlasStyle->SetCanvasBorderMode(icol);
atlasStyle->SetPadBorderMode(icol);
atlasStyle->SetPadColor(icol);
atlasStyle->SetCanvasColor(icol);
atlasStyle->SetStatColor(icol);
//atlasStyle->SetFillColor(icol);

// set the paper & margin sizes
atlasStyle->SetPaperSize(20,26);
atlasStyle->SetPadTopMargin(0.1);
//atlasStyle->SetPadRightMargin(0.05);
atlasStyle->SetPadRightMargin(0.12);
atlasStyle->SetPadBottomMargin(0.16);
atlasStyle->SetPadLeftMargin(0.12);

// use large fonts
//Int_t font=72;
Int_t font=32;
Double_t tsize=0.05;
atlasStyle->SetTextFont(font);


atlasStyle->SetTextSize(tsize);
atlasStyle->SetLabelFont(font,"x");
atlasStyle->SetTitleFont(font,"x");
atlasStyle->SetLabelFont(font,"y");
atlasStyle->SetTitleFont(font,"y");
atlasStyle->SetLabelFont(font,"z");
atlasStyle->SetTitleFont(font,"z");

atlasStyle->SetLabelSize(tsize,"x");
atlasStyle->SetTitleSize(tsize,"x");
atlasStyle->SetLabelSize(tsize,"y");
atlasStyle->SetTitleSize(tsize,"y");
atlasStyle->SetLabelSize(tsize,"z");
atlasStyle->SetTitleSize(tsize,"z");
//atlasStyle->SetTitleOffset(1.2,"y");

//use bold lines and markers
atlasStyle->SetMarkerStyle(20);
atlasStyle->SetMarkerSize(1.2);
atlasStyle->SetHistLineWidth(2.);
atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

//get rid of X error bars and y error bar caps
//atlasStyle->SetErrorX(0.001);

//do not display any of the standard histogram decorations
//atlasStyle->SetOptTitle(0);
//atlasStyle->SetOptStat(1111);
atlasStyle->SetOptStat(0);
//atlasStyle->SetOptFit(1111);
atlasStyle->SetOptFit(0);

// put tick marks on top and RHS of plots
atlasStyle->SetPadTickX(1);
atlasStyle->SetPadTickY(1);
 

gROOT->SetStyle("Plain");

//gStyle->SetPadTickX(1);
//gStyle->SetPadTickY(1);

}


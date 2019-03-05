#include "TBReader.h"
#include "MakePlots.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"

TBReader::TBReader( TChain *c1, TChain *c2, string filename):T_Rechit(c1),T_DWC(c2){
  cout << "Constructor of TBReader, Test Beam Run ... \n\n" << endl;
  fname = filename;
}

TBReader::~TBReader(){
  cout << "\n\n";
  cout << "Destructor of TBReader ... " << endl;
}

void TBReader::Init_Pointers(){
  rechit_detid = 0;
  rechit_module = 0;
  rechit_layer = 0;
  rechit_chip = 0;
  rechit_channel = 0;
  rechit_type = 0;
  rechit_x = 0;
  rechit_y = 0;
  rechit_z = 0;
  rechit_iu = 0;
  rechit_iv = 0;
  rechit_energy = 0;
  rechit_energy_noHG = 0;
  rechit_amplitudeHigh = 0;
  rechit_amplitudeLow = 0;
  rechit_hg_goodFit = 0;
  rechit_lg_goodFit = 0;
  rechit_hg_saturated = 0;
  rechit_lg_saturated = 0;
  rechit_fully_calibrated = 0;
  rechit_TS2High = 0;
  rechit_TS2Low = 0;
  rechit_TS3High = 0;
  rechit_TS3Low = 0;
  rechit_Tot = 0;
  rechit_time = 0;
  rechit_timeMaxHG = 0;
  rechit_timeMaxLG = 0;
  rechit_toaRise = 0;
  rechit_toaFall = 0;
   
}

void TBReader::SetRootBranch(){
  nevents = T_Rechit->GetEntries();

  T_Rechit->SetBranchAddress("event", &event);
  T_Rechit->SetBranchAddress("run", &run);
  T_Rechit->SetBranchAddress("pdgID", &pdgID);
  T_Rechit->SetBranchAddress("beamEnergy", &beamEnergy);
  T_Rechit->SetBranchAddress("trueBeamEnergy", &trueBeamEnergy);
  //T_Rechit->SetBranchAddress("NRechits", &NRechits);
  T_Rechit->SetBranchAddress("rechit_detid", &rechit_detid);
  T_Rechit->SetBranchAddress("rechit_module", &rechit_module);
  T_Rechit->SetBranchAddress("rechit_layer", &rechit_layer);
  T_Rechit->SetBranchAddress("rechit_chip", &rechit_chip);
  T_Rechit->SetBranchAddress("rechit_channel", &rechit_channel);
  T_Rechit->SetBranchAddress("rechit_type", &rechit_type);

  T_Rechit->SetBranchAddress("rechit_x", &rechit_x);
  T_Rechit->SetBranchAddress("rechit_y", &rechit_y);
  T_Rechit->SetBranchAddress("rechit_z", &rechit_z);
  T_Rechit->SetBranchAddress("rechit_iu", &rechit_iu);
  T_Rechit->SetBranchAddress("rechit_iv", &rechit_iv);
  T_Rechit->SetBranchAddress("rechit_energy", &rechit_energy);
  T_Rechit->SetBranchAddress("rechit_energy_noHG", &rechit_energy_noHG);

  T_Rechit->SetBranchAddress("rechit_amplitudeHigh", &rechit_amplitudeHigh);
  T_Rechit->SetBranchAddress("rechit_amplitudeLow", &rechit_amplitudeLow);
  T_Rechit->SetBranchAddress("rechit_hg_goodFit", &rechit_hg_goodFit);
  T_Rechit->SetBranchAddress("rechit_lg_goodFit", &rechit_lg_goodFit);
  T_Rechit->SetBranchAddress("rechit_hg_saturated", &rechit_hg_saturated);
  T_Rechit->SetBranchAddress("rechit_lg_saturated", &rechit_lg_saturated);
  T_Rechit->SetBranchAddress("rechit_fully_calibrated", &rechit_fully_calibrated);
  T_Rechit->SetBranchAddress("rechit_TS2High", &rechit_TS2High);
  T_Rechit->SetBranchAddress("rechit_TS2Low", &rechit_TS2Low);
  T_Rechit->SetBranchAddress("rechit_TS3High", &rechit_TS3High);
  T_Rechit->SetBranchAddress("rechit_TS3Low", &rechit_TS3Low);
    
  T_Rechit->SetBranchAddress("rechit_Tot", &rechit_Tot);
  T_Rechit->SetBranchAddress("rechit_time", &rechit_time);
  T_Rechit->SetBranchAddress("rechit_timeMaxHG", &rechit_timeMaxHG);
  T_Rechit->SetBranchAddress("rechit_timeMaxLG", &rechit_timeMaxLG);
  T_Rechit->SetBranchAddress("rechit_toaRise", &rechit_toaRise);
  T_Rechit->SetBranchAddress("rechit_toaFall", &rechit_toaFall);

  
  T_DWC->SetBranchAddress("ntracks", &ntracks);
  T_DWC->SetBranchAddress("trackChi2_X", &trackChi2_X);
  T_DWC->SetBranchAddress("trackChi2_Y", &trackChi2_Y);
  T_DWC->SetBranchAddress("dwcReferenceType", &dwcReferenceType);
  T_DWC->SetBranchAddress("m_x", &m_x);
  T_DWC->SetBranchAddress("m_y", &m_y);
  T_DWC->SetBranchAddress("b_x", &b_x);
  T_DWC->SetBranchAddress("b_y", &b_y);

}

void TBReader::Init_Beaminfo(bool quite){
  T_Rechit->GetEntry(0);
  beamE = beamEnergy;
  if( pdgID == 11 ){
    beam_str = "Electron";
    PID = 0;}
  else if( pdgID == 13){
    beam_str = "Muon";
    PID = 2;}
  else if( pdgID == 211){
    beam_str = "Pion";
    PID = 1;}
  else{
    cout << "unknown PDGID QQ" << endl;
    beam_str = "??";
    PID = -1;}
  RunN = run;
  if(!quite){
    cout << beam_str.c_str() << " Run, "<< beamE << "GeV , with " << nevents
	 << " events." << endl;}
}

bool TBReader::Check_Config(int given_config){
  Init(true);
  int setup_config;
  if(RunN <= 722 ) {
    setup_config = 1;
    if(RunN <= 257){
      cout << "TOA threshold will be changed after Run 257." << endl;}
  }
  else if(RunN > 722 && RunN <= 1057){
    setup_config = 2;}
  else if(RunN > 1057 && RunN <= 1078){
    setup_config = 5;
    return false; //Not going to handle muon run
  }
  else{
    setup_config = 3;
  }

  if(setup_config == given_config){ return true; }
  else{
    cout << "conflict config with main function, please modify the main config."
	 << endl;
    return false;}
  
}

void TBReader::Init(bool quite){
  Init_Pointers();
  SetRootBranch();
  Init_Beaminfo(quite);
}

void TBReader::Ntuple_Maker(setup_config *SC){
  
  string outpath = string( dirpath + string("Module_Ntuple") );
  TFile *outNtuple[MAXBOARDS];
  TTree *outTree[MAXBOARDS];
  TTree *outTree_history[MAXBOARDS];
  for(int i = 0 ; i < MAXBOARDS ; ++i){
    outNtuple[i] = NULL;
    outTree  [i] = NULL;
    outTree_history[i] = NULL; }

  // Variable for output trees
  // vector<unsigned int> *TB_chip;
  // vector<unsigned int> *TB_channel;
  // vector<unsigned int> *TB_type;
  
  for(int ifile = 0 ; ifile < 1 ; ifile++){
    char fpath[200];
    
    // Check root file exist
    int moduleID = SC->Module_List[ifile];
    if(moduleID == 0) { continue; }
    sprintf(fpath,"%s/TB/Module%d_Oct18.root",outpath.c_str(),moduleID);
    
    ifstream f_check(fpath);
    if( !f_check.good() ){
      outNtuple[ifile] = new TFile(fpath,"recreate");
      outTree[ifile]   = new TTree("TB","TB");

      //TODO: Create outtree branch

      //Fill Run history
      outTree_history[ifile] = new TTree("Run_history","Run_history");
      int Run_hist;
      outTree_history[ifile]->Branch("RunNumber",&Run_hist);
      Run_hist = RunN; 
      outTree_history[ifile]->Fill();
    }
    else{
      outNtuple[ifile] = new TFile(fpath,"update");
      outTree[ifile]   = (TTree*)outNtuple[ifile]->Get("TB");
      //TODO: Set outtree branch Address
      
      // Check and update run_history      
      outTree_history[ifile] = (TTree*)outNtuple[ifile]->Get("Run_history");
      bool already_filled = Check_run_filled(outTree_history[ifile]);
      if(already_filled){
	cout << "Run " << RunN << " has already filled in "
	     << fpath << ", skip it!"<< endl;
	outTree[ifile] = NULL;}
    }
  }
  
  for(int ev = 0 ; ev < nevents ; ++ev){
    T_Rechit->GetEntry(ev);
    //T_DWC->GetEntry(ev);
    
    //Event selection
    //if(dwcReferenceType != 15) continue;
      
    for(int ifile = 0 ; ifile < 1 ; ifile++){
      if(outTree[ifile] != NULL){	
	// InitoutTreeBranch();
	// SetoutTreeBranchAddress();
	
	outTree[ifile]->Fill();}
    }
  }
  
  for(int ifile = 0 ; ifile < 1 ; ifile++){
    if(outTree[ifile] != NULL){
      outNtuple[ifile]->cd();
      outTree[ifile]->Write("TB",TObject::kOverwrite);
      outTree_history[ifile]->Write("Run_history",TObject::kOverwrite);
    }
    outNtuple[ifile]->Close();
  }
  

}

void TBReader::TProfile_Maker(setup_config *SC,MakePlots *M){
  
  Init();
  //Run selection
  if( PID == 2 ){ cout << "Not e-/pi runs, skip Run "
		       << RunN << " for now." << endl; return; }
  

  bool double_fill = M->Check_Run(RunN);
  if(double_fill) {
    cout << "Run " << RunN << " already filled in "
	 << M-> TPro_output << " skip it!" << endl;
    return;  }  

  cout << "Looping evts" << endl;
  for(int ev = 0 ; ev < nevents ; ++ev){
    T_Rechit->GetEntry(ev);
    //T_DWC->GetEntry(ev);
    
    //Event selection
    //if(dwcReferenceType != 15) continue;    

    for(int ihit = 0 ; ihit < (int)rechit_amplitudeHigh->size() ; ++ihit){
      
      double HG,LG,TOT;
      int chip,ch,BD_order,moduleID;
      HG   = rechit_amplitudeHigh->at(ihit);
      LG   = rechit_amplitudeLow->at(ihit);
      TOT  = rechit_Tot->at(ihit);
      chip = (int)rechit_chip->at(ihit);
      ch   = (int)rechit_channel->at(ihit);
      ch   /= 2;
      moduleID = rechit_module->at(ihit);
      if ( SC->moduleID2BDorder.find(moduleID) == SC->moduleID2BDorder.end() ){
	continue; } // Not finding any moduleID in the mapping -> continue
      else{
	BD_order = SC->moduleID2BDorder.find(moduleID)->second;}
      if( LG < 5 ) continue;
      //      M->HG_LG[BD_order][chip][ch]->Fill(LG,HG,1);
      M->h_HG_LG[BD_order][chip][ch]->Fill(LG,HG,1);
      if( TOT < 5 ) continue;
      //      M->LG_TOT[BD_order][chip][ch]->Fill(TOT,LG,1);      
      //M->h_LG_TOT[BD_order][chip][ch]->Fill(TOT,LG,1);
    }    
  }
  cout << "End of looping evts" << endl;
}

/*
double TBReader::Calculate_medium(double medium){
}
*/

bool TBReader::Check_run_filled(TTree* tree){
  //Check if this run have already been filled
  int Filled_Runs = tree-> GetEntries();
  int Run_hist;
  tree->SetBranchAddress("RunNumber",&Run_hist);
  bool already_filled = false;
  for(int irun = 0 ; irun < Filled_Runs ; ++ irun){
    tree->GetEntry(irun);
    if(Run_hist == RunN){
      already_filled = true;
    }
  }
  if(!already_filled){
    Run_hist = RunN;
    tree->Fill();
  }
  return already_filled;
}


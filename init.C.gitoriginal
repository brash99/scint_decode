#include "includes.h"

int run; // run number

Int_t connected[NBARS];

TTree *t;
Int_t raw[NUMADCSLOTS][NUMCHANA], ped[NUMADCSLOTS][NUMCHANA], trail[NUMTDCSLOTS][NUMCHANT], lead[NUMTDCSLOTS][NUMCHANT];




void init(){
  cout << "Initializations done. To read a run:" << endl <<"read(RUN NUMBER)"<<endl;
}

void read(int runno){
	TString filename;
	filename.Form("sbs_%d_14.root",runno);
	TFile *_file0 = TFile::Open(filename);

	run=runno;

	// is bar connected? 1 if yes, 0 if no
	for (int i=0; i<NBARS; i++) connected[i]=1; // default = connected
	// any exception here, for example
	// connected[0] = 0; 
	// means first bar is not connected to any pixel.
	// See below for mappings


	TStyle *MyStyle = new TStyle("MyStyle","MyStyle");
	MyStyle->SetTitleFontSize(0.08);
	MyStyle->SetTitleX(0.15);
	MyStyle->SetTitleY(0.99);
	MyStyle->SetStatW(0.2);
	MyStyle->SetMarkerStyle(21);
	MyStyle->SetOptFit(0111);
	gROOT->SetStyle("MyStyle");

	// read tree
	t = (TTree *)_file0->Get("t"); 
	t->SetBranchAddress("adcraw", raw);
	t->SetBranchAddress("ped", ped);
	t->SetBranchAddress("tdct", trail);
	t->SetBranchAddress("tdcl", lead);

	cout << "Found " << t->GetEntries() <<" events"<<endl;
}

// scripts (yes, this goes after the previous declarations)
#include "mapping.C"
#include "Avert.C"
#include "checkpixels.C"
#include "verticalfit.C"
#include "NINOgains.C"


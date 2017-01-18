static const int NUMADCSLOTS = 4;
static const int NUMTDCSLOTS = 3;
static const int NUMCHANA = 64;
static const int NUMCHANT = 96;
static const float adc_charge = 50*1e-15; // 50 fC, corresponding to an ADC chan
static const float e = 1.6e-19; // C, electron charge
static const int NPMT = 14;
static const int NPIXELTOT = 16;
static const int NPIXEL = NPIXELTOT; // not connected pixels not known at date
static const int NBARS = NPMT*NPIXEL;

Double_t avg[NBARS];
Double_t bar[NBARS];

// root stuff
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "Riostream.h"
#include "TTree.h"
#include "TEfficiency.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"

// for fit on eta:
float Athresholds[NUMADCSLOTS][NUMCHANA];
float Awidths[NUMADCSLOTS][NUMCHANA];

// for fit on the vertical tracks:
float peak_coeff_range = 2;
float peak_coeff_mean_min = -2;
float peak_coeff_mean_max = +2;
float peak_coeff_rms_min = 0;
float peak_coeff_rms_max = 2;
float Peak_mean[NUMADCSLOTS][NUMCHANA]; // init&filling in plot_peak,plot_peak2
float Peak_rms[NUMADCSLOTS][NUMCHANA];  // init&filling in plot_peak,plot_peak2

// for noise frequency:
float last_tdc_cut;
TH1D *Noise_tdc[NUMTDCSLOTS][NUMCHANT]; // init&filling in plot_tdc_vs_adc

int run; // run number. Used in the titles

TTree *t;
Int_t raw[NUMADCSLOTS][NUMCHANA], ped[NUMADCSLOTS][NUMCHANA], trail[NUMTDCSLOTS][NUMCHANT], lead[NUMTDCSLOTS][NUMCHANT];


// mapping pmt,pixel -> adc/tdc slot,chan
Int_t handmapping_adc_slot(Int_t group);
Int_t handmapping_adc_chan(Int_t group, Int_t pixel);
Int_t handmapping_tdc_slot(Int_t group);
Int_t handmapping_tdc_chan(Int_t group, Int_t pixel);

#include "mapping.C"
//#include "Acuts.C"


//#include <fstream.h>

// For interactive use:
// reads root file corresponding to run runno, extracts the tree with
// raw data, define shortcuts for slot 0

// TODO? Should probably map the aliases in a way that corresponding
// TDC and ADC channels have the same logical slot and channel
// (hardware and software match naturally for slot 0 chans 0-47..not
// so for everything else!)

// IMPORTANT: AT THE MOMENT, ONLY SLOT 0 CHANS 0-47 HAVE SHORTCUTS!!!
// Note: in the general function section (way below) I added functions
// to convert slots and chans ADC from/to TDC. They are used in the
// plot_* routines, but I am not sure how to meaningfully use them for
// the following aliases
void init(int runno){
	TString filename;
	filename.Form("sbs_%d_14.root",runno);
	TFile *_file0 = TFile::Open(filename);

	run=runno;

	TStyle *MyStyle = new TStyle("MyStyle","MyStyle");
	//MyStyle->SetStatFont(42);
	MyStyle->SetTitleFontSize(0.08);
	MyStyle->SetTitleX(0.15);
	MyStyle->SetTitleY(0.99);
	MyStyle->SetStatW(0.4);
	MyStyle->SetMarkerStyle(6);
	gROOT->SetStyle("MyStyle");

	// read tree
	t = (TTree *)_file0->Get("t"); 
	t->SetBranchAddress("adcraw", raw);
	t->SetBranchAddress("ped", ped);
	t->SetBranchAddress("tdct", trail);
	t->SetBranchAddress("tdcl", lead);

	cout << "Found " << t->GetEntries() <<" events"<<endl;

	// aliases for slot0
	//==========
	// ADC slot 0 channels 0-47
	// t->SetAlias("a0", "adc[0][0]");
	// t->SetAlias("a1", "adc[0][1]");
	// t->SetAlias("a2", "adc[0][2]");
	// t->SetAlias("a3", "adc[0][3]");
	// t->SetAlias("a4", "adc[0][4]");
	// t->SetAlias("a5", "adc[0][5]");
	// t->SetAlias("a6", "adc[0][6]");
	// t->SetAlias("a7", "adc[0][7]");
	// t->SetAlias("a8", "adc[0][8]");
	// t->SetAlias("a9", "adc[0][9]");
	// t->SetAlias("a10","adc[0][10]");
	// t->SetAlias("a11","adc[0][11]");
	// t->SetAlias("a12","adc[0][12]");
	// t->SetAlias("a13","adc[0][13]");
	// t->SetAlias("a14","adc[0][14]");
	// t->SetAlias("a15","adc[0][15]");
	// t->SetAlias("a16","adc[0][16]");
	// t->SetAlias("a17","adc[0][17]");
	// t->SetAlias("a18","adc[0][18]");
	// t->SetAlias("a19","adc[0][19]");
	// t->SetAlias("a20","adc[0][20]");
	// t->SetAlias("a21","adc[0][21]");
	// t->SetAlias("a22","adc[0][22]");
	// t->SetAlias("a23","adc[0][23]");
	// t->SetAlias("a24","adc[0][24]");
	// t->SetAlias("a25","adc[0][25]");
	// t->SetAlias("a26","adc[0][26]");
	// t->SetAlias("a27","adc[0][27]");
	// t->SetAlias("a28","adc[0][28]");
	// t->SetAlias("a29","adc[0][29]");
	// t->SetAlias("a30","adc[0][30]");
	// t->SetAlias("a31","adc[0][31]");
	// t->SetAlias("a32","adc[0][32]");
	// t->SetAlias("a33","adc[0][33]");
	// t->SetAlias("a34","adc[0][34]");
	// t->SetAlias("a35","adc[0][35]");
	// t->SetAlias("a36","adc[0][36]");
	// t->SetAlias("a37","adc[0][37]");
	// t->SetAlias("a38","adc[0][38]");
	// t->SetAlias("a39","adc[0][39]");
	// t->SetAlias("a40","adc[0][40]");
	// t->SetAlias("a41","adc[0][41]");
	// t->SetAlias("a42","adc[0][42]");
	// t->SetAlias("a43","adc[0][43]");
	// t->SetAlias("a44","adc[0][44]");
	// t->SetAlias("a45","adc[0][45]");
	// t->SetAlias("a46","adc[0][46]");
	// t->SetAlias("a47","adc[0][47]");

	// // TDC slot 0 channels 0-47
	// // number of hits
	// t->SetAlias("n0", "tdcn[0][0]");
	// t->SetAlias("n1", "tdcn[0][1]");
	// t->SetAlias("n2", "tdcn[0][2]");
	// t->SetAlias("n3", "tdcn[0][3]");
	// t->SetAlias("n4", "tdcn[0][4]");
	// t->SetAlias("n5", "tdcn[0][5]");
	// t->SetAlias("n6", "tdcn[0][6]");
	// t->SetAlias("n7", "tdcn[0][7]");
	// t->SetAlias("n8", "tdcn[0][8]");
	// t->SetAlias("n9", "tdcn[0][9]");
	// t->SetAlias("n10","tdcn[0][10]");
	// t->SetAlias("n11","tdcn[0][11]");
	// t->SetAlias("n12","tdcn[0][12]");
	// t->SetAlias("n13","tdcn[0][13]");
	// t->SetAlias("n14","tdcn[0][14]");
	// t->SetAlias("n15","tdcn[0][15]");
	// t->SetAlias("n16","tdcn[0][16]");
	// t->SetAlias("n17","tdcn[0][17]");
	// t->SetAlias("n18","tdcn[0][18]");
	// t->SetAlias("n19","tdcn[0][19]");
	// t->SetAlias("n20","tdcn[0][20]");
	// t->SetAlias("n21","tdcn[0][21]");
	// t->SetAlias("n22","tdcn[0][22]");
	// t->SetAlias("n23","tdcn[0][23]");
	// t->SetAlias("n24","tdcn[0][24]");
	// t->SetAlias("n25","tdcn[0][25]");
	// t->SetAlias("n26","tdcn[0][26]");
	// t->SetAlias("n27","tdcn[0][27]");
	// t->SetAlias("n28","tdcn[0][28]");
	// t->SetAlias("n29","tdcn[0][29]");
	// t->SetAlias("n30","tdcn[0][30]");
	// t->SetAlias("n31","tdcn[0][31]");
	// t->SetAlias("n32","tdcn[0][32]");
	// t->SetAlias("n33","tdcn[0][33]");
	// t->SetAlias("n34","tdcn[0][34]");
	// t->SetAlias("n35","tdcn[0][35]");
	// t->SetAlias("n36","tdcn[0][36]");
	// t->SetAlias("n37","tdcn[0][37]");
	// t->SetAlias("n38","tdcn[0][38]");
	// t->SetAlias("n39","tdcn[0][39]");
	// t->SetAlias("n40","tdcn[0][40]");
	// t->SetAlias("n41","tdcn[0][41]");
	// t->SetAlias("n42","tdcn[0][42]");
	// t->SetAlias("n43","tdcn[0][43]");
	// t->SetAlias("n44","tdcn[0][44]");
	// t->SetAlias("n45","tdcn[0][45]");
	// t->SetAlias("n46","tdcn[0][46]");
	// t->SetAlias("n47","tdcn[0][47]");

	// // trailing edges
	// t->SetAlias("t0", "tdct[0][0]");
	// t->SetAlias("t1", "tdct[0][1]");
	// t->SetAlias("t2", "tdct[0][2]");
	// t->SetAlias("t3", "tdct[0][3]");
	// t->SetAlias("t4", "tdct[0][4]");
	// t->SetAlias("t5", "tdct[0][5]");
	// t->SetAlias("t6", "tdct[0][6]");
	// t->SetAlias("t7", "tdct[0][7]");
	// t->SetAlias("t8", "tdct[0][8]");
	// t->SetAlias("t9", "tdct[0][9]");
	// t->SetAlias("t10","tdct[0][10]");
	// t->SetAlias("t11","tdct[0][11]");
	// t->SetAlias("t12","tdct[0][12]");
	// t->SetAlias("t13","tdct[0][13]");
	// t->SetAlias("t14","tdct[0][14]");
	// t->SetAlias("t15","tdct[0][15]");
	// t->SetAlias("t16","tdct[0][16]");
	// t->SetAlias("t17","tdct[0][17]");
	// t->SetAlias("t18","tdct[0][18]");
	// t->SetAlias("t19","tdct[0][19]");
	// t->SetAlias("t20","tdct[0][20]");
	// t->SetAlias("t21","tdct[0][21]");
	// t->SetAlias("t22","tdct[0][22]");
	// t->SetAlias("t23","tdct[0][23]");
	// t->SetAlias("t24","tdct[0][24]");
	// t->SetAlias("t25","tdct[0][25]");
	// t->SetAlias("t26","tdct[0][26]");
	// t->SetAlias("t27","tdct[0][27]");
	// t->SetAlias("t28","tdct[0][28]");
	// t->SetAlias("t29","tdct[0][29]");
	// t->SetAlias("t30","tdct[0][30]");
	// t->SetAlias("t31","tdct[0][31]");
	// t->SetAlias("t32","tdct[0][32]");
	// t->SetAlias("t33","tdct[0][33]");
	// t->SetAlias("t34","tdct[0][34]");
	// t->SetAlias("t35","tdct[0][35]");
	// t->SetAlias("t36","tdct[0][36]");
	// t->SetAlias("t37","tdct[0][37]");
	// t->SetAlias("t38","tdct[0][38]");
	// t->SetAlias("t39","tdct[0][39]");
	// t->SetAlias("t40","tdct[0][40]");
	// t->SetAlias("t41","tdct[0][41]");
	// t->SetAlias("t42","tdct[0][42]");
	// t->SetAlias("t43","tdct[0][43]");
	// t->SetAlias("t44","tdct[0][44]");
	// t->SetAlias("t45","tdct[0][45]");
	// t->SetAlias("t46","tdct[0][46]");
	// t->SetAlias("t47","tdct[0][47]");

	// // leading edges
	// t->SetAlias("l0", "tdcl[0][0]");
	// t->SetAlias("l1", "tdcl[0][1]");
	// t->SetAlias("l2", "tdcl[0][2]");
	// t->SetAlias("l3", "tdcl[0][3]");
	// t->SetAlias("l4", "tdcl[0][4]");
	// t->SetAlias("l5", "tdcl[0][5]");
	// t->SetAlias("l6", "tdcl[0][6]");
	// t->SetAlias("l7", "tdcl[0][7]");
	// t->SetAlias("l8", "tdcl[0][8]");
	// t->SetAlias("l9", "tdcl[0][9]");
	// t->SetAlias("l10","tdcl[0][10]");
	// t->SetAlias("l11","tdcl[0][11]");
	// t->SetAlias("l12","tdcl[0][12]");
	// t->SetAlias("l13","tdcl[0][13]");
	// t->SetAlias("l14","tdcl[0][14]");
	// t->SetAlias("l15","tdcl[0][15]");
	// t->SetAlias("l16","tdcl[0][16]");
	// t->SetAlias("l17","tdcl[0][17]");
	// t->SetAlias("l18","tdcl[0][18]");
	// t->SetAlias("l19","tdcl[0][19]");
	// t->SetAlias("l20","tdcl[0][20]");
	// t->SetAlias("l21","tdcl[0][21]");
	// t->SetAlias("l22","tdcl[0][22]");
	// t->SetAlias("l23","tdcl[0][23]");
	// t->SetAlias("l24","tdcl[0][24]");
	// t->SetAlias("l25","tdcl[0][25]");
	// t->SetAlias("l26","tdcl[0][26]");
	// t->SetAlias("l27","tdcl[0][27]");
	// t->SetAlias("l28","tdcl[0][28]");
	// t->SetAlias("l29","tdcl[0][29]");
	// t->SetAlias("l30","tdcl[0][30]");
	// t->SetAlias("l31","tdcl[0][31]");
	// t->SetAlias("l32","tdcl[0][32]");
	// t->SetAlias("l33","tdcl[0][33]");
	// t->SetAlias("l34","tdcl[0][34]");
	// t->SetAlias("l35","tdcl[0][35]");
	// t->SetAlias("l36","tdcl[0][36]");
	// t->SetAlias("l37","tdcl[0][37]");
	// t->SetAlias("l38","tdcl[0][38]");
	// t->SetAlias("l39","tdcl[0][39]");
	// t->SetAlias("l40","tdcl[0][40]");
	// t->SetAlias("l41","tdcl[0][41]");
	// t->SetAlias("l42","tdcl[0][42]");
	// t->SetAlias("l43","tdcl[0][43]");
	// t->SetAlias("l44","tdcl[0][44]");
	// t->SetAlias("l45","tdcl[0][45]");
	// t->SetAlias("l46","tdcl[0][46]");
	// t->SetAlias("l47","tdcl[0][47]");
}

//====================================================================
//                     Multiplots section
//====================================================================

//--------------------------------------------------------------------
//                        ADC only
//--------------------------------------------------------------------
// For interactive use:
// plot canvas divided in 4x4 ADC plots of 16 channels starting from
// adc_chan_start, of given slot. ADC is plotted with and without the
// following cut: tdc value > tdc_cut
// By default, tdc starting chan and slot are calculated
// automatically; if a different values than -1 to
// tdc_chan_start/slot, those values will be used instead. May be
// useful during testing of ADC/TDC channels
TCanvas *plot_adc(int adc_slot, int adc_chan_start, int tdc_cut=1350, int tdc_slot = -1, int tdc_chan_start = -1){
	TString cut, draw, draw1, title;
	title.Form("run_%d_ADC",run);
	TCanvas *cADC= new TCanvas("cADC",title,850,850);
	int nbin=100;
	int min=-100, max=500;
	TH1D *htmp=new TH1D("htmp","htmp",nbin,min,max);
	TH1D *htmp1=new TH1D("htmp1","htmp1",nbin,min,max);
	TPad *current=0;
	TPaveStats *ps = 0;
	TList *list = 0;
	TText *tconst = 0;
	TLatex *myt = 0;
	TString tmpentry;
	MyStyle->SetStatX(0.9);
	MyStyle->SetStatY(0.6);
	last_tdc_cut = tdc_cut;

	cADC->Divide(4,4) ;
	// log scale in number of counts, so that the pedestal events
	// do not push the others to the lowest part of the color
	// spectrum
	cADC_1->SetLogy();
	cADC_2->SetLogy();
	cADC_3->SetLogy();
	cADC_4->SetLogy();
	cADC_5->SetLogy();
	cADC_6->SetLogy();
	cADC_7->SetLogy();
	cADC_8->SetLogy();
	cADC_9->SetLogy();
	cADC_10->SetLogy();
	cADC_11->SetLogy();
	cADC_12->SetLogy();
	cADC_13->SetLogy();
	cADC_14->SetLogy();
	cADC_15->SetLogy();
	cADC_16->SetLogy();

	// 16 channels to plot ( = 1 PMT)
	for (int i=adc_chan_start; i<adc_chan_start+16; i++){ 

	  if (tdc_chan_start == -1) tdc_chan_start=calc_tdc_chan(adc_slot,i);
	  if (tdc_slot == -1) tdc_slot = calc_tdc_slot(adc_slot,i);

	  int itdc = tdc_chan_start + i-adc_chan_start;
	  cADC->cd(i-adc_chan_start+1);
	  cout << "Doing plot " << i-adc_chan_start+1 << " / 16"<<endl;
	  htmp=new TH1D("htmp","htmp",nbin,min,max);
	  htmp1=new TH1D("htmp1","htmp1",nbin,min,max);
	  draw.Form("adc[%d][%d]>>htmp",adc_slot,i);
	  draw1.Form("adc[%d][%d]>>htmp1",adc_slot,i);
	  cut.Form("tdct[%d][%d]>%d",tdc_slot,itdc,tdc_cut);
	  htmp->SetLineColor(kBlue);
	  htmp1->SetLineColor(kRed);
	  title.Form("Run %d ADC slot %d chan %d tdc > %d",run,adc_slot, i, tdc_cut);
	  htmp->SetTitle(title);
	  htmp1->SetTitle(title);
	  
	  t->Draw(draw);
	  t->Draw(draw1,cut,"same");

	  cADC->Update();
	  // Retrieve the stat box
	  switch(i-adc_chan_start+1){
	  case 1:
	    current=cADC_1; break;
	  case 2:
	    current=cADC_2; break;
	  case 3:
	    current=cADC_3; break;
	  case 4:
	    current=cADC_4; break;
	  case 5:
	    current=cADC_5; break;
	  case 6:
	    current=cADC_6; break;
	  case 7:
	    current=cADC_7; break;
	  case 8:
	    current=cADC_8; break;
	  case 9:
	    current=cADC_9; break;
	  case 10:
	    current=cADC_10; break;
	  case 11:
	    current=cADC_11; break;
	  case 12:
	    current=cADC_12; break;
	  case 13:
	    current=cADC_13; break;
	  case 14:
	    current=cADC_14; break;
	  case 15:
	    current=cADC_15; break;
	  case 16:
	    current=cADC_16; break;
	  default:
	    cout << "Canvas number error" << endl;
	  }
	  int entries = htmp->GetEntries();
	  float mean = htmp->GetMean(1);
	  float RMS = htmp->GetRMS(1);

	  ps =(TPaveStats*)current->GetPrimitive("stats");
	  ps->SetName("mystats");
	  list = ps->GetListOfLines();
	  tconst = ps->GetLineWith("Entries"); list->Remove(tconst);
	  tconst = ps->GetLineWith("Mean"); list->Remove(tconst);
	  tconst = ps->GetLineWith("RMS"); list->Remove(tconst);

	  tmpentry.Form("Entries = %d",entries); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  tmpentry.Form("Mean = %.3f",mean); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  tmpentry.Form("RMS = %.3f",RMS); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);

	  htmp->SetStats(0);
	  current->Modified();
	};
	
	title.Form("run_%d_ADC_slot_%d_chan_%d_%d_tdc_cut_%d.png",
		   run,adc_slot,adc_chan_start,adc_chan_start+15,tdc_cut);
	cADC->Print(title);
	cADC->cd(0);
	return cADC;
}

//--------------------------------------------------------------------
//                        TDC only
//--------------------------------------------------------------------
// For interactive use:
// plotcanvas divided in 4x4 TDC plots of 16 channels starting from
// adc_chan_start, of given slot. TDC is plotted with and without the
// following cut: adc value > adc_cut
// By default, tdc starting chan = adc starting chan; if passed a
// different value than -1 to adc_chan_start, that value will be used
// instead. May be useful during testing of ADC/TDC channels
// Analogously for adc_slot
TCanvas *plot_tdc(int tdc_slot, int tdc_chan_start, int adc_cut, int adc_slot = -1, int adc_chan_start = -1){
	TString cut, draw, draw1, title;
	title.Form("run_%d_TDC",run);
	TCanvas *cTDC= new TCanvas("cTDC",title,850,850);
	int nbin=100;
	int min=-100, max=1500;
	TH1D *htmp=new TH1D("htmp","htmp",nbin,min,max);
	TH1D *htmp1=new TH1D("htmp1","htmp1",nbin,min,max);
	TString tmpentry;
	TPad *current=0;
	TPaveStats *ps = 0;
	TList *list = 0;
	TText *tconst = 0;
	TLatex *myt = 0;
	MyStyle->SetStatX(0.6);
	MyStyle->SetStatY(0.6);

	cTDC->Divide(4,4) ;
	// log scale in number of counts, so that the pedestal events
	// do not push the others to the lowest part of the color
	// spectrum
	cTDC_1->SetLogy();
	cTDC_2->SetLogy();
	cTDC_3->SetLogy();
	cTDC_4->SetLogy();
	cTDC_5->SetLogy();
	cTDC_6->SetLogy();
	cTDC_7->SetLogy();
	cTDC_8->SetLogy();
	cTDC_9->SetLogy();
	cTDC_10->SetLogy();
	cTDC_11->SetLogy();
	cTDC_12->SetLogy();
	cTDC_13->SetLogy();
	cTDC_14->SetLogy();
	cTDC_15->SetLogy();
	cTDC_16->SetLogy();

	// 16 channels to plot
	for (int i=tdc_chan_start; i<tdc_chan_start+16; i++){

	  if (adc_chan_start == -1) adc_chan_start = calc_adc_chan(tdc_slot,tdc_chan_start);
	  if (adc_slot == -1) adc_slot = adc_adc_slot(tdc_slot, tdc_chan_start);

	  cTDC->cd(i-tdc_chan_start+1);
	  int iadc = adc_chan_start + i-tdc_chan_start;
	  cout << "Doing plot " << i-tdc_chan_start+1 << " / 16"<<endl;
	  htmp=new TH1D("htmp","htmp",nbin,min,max);
	  htmp1=new TH1D("htmp1","htmp1",nbin,min,max);
	  draw.Form("tdct[%d][%d]>>htmp",tdc_slot,i);
	  draw1.Form("tdct[%d][%d]>>htmp1",tdc_slot,i);
	  cut.Form("adc[%d][%d]>%d",adc_slot,iadc,adc_cut);
	  htmp->SetLineColor(kBlue);
	  htmp1->SetLineColor(kRed);
	  title.Form("Run %d TDC slot %d chan %d adc > %d",run,tdc_slot, i, adc_cut);
	  htmp->SetTitle(title);
	  htmp1->SetTitle(title);
	  
	  t->Draw(draw);
	  t->Draw(draw1,cut,"same");

	  cTDC->Update();
	  // Retrieve the stat box
	  switch(i-tdc_chan_start+1){
	  case 1:
	    current=cTDC_1; break;
	  case 2:
	    current=cTDC_2; break;
	  case 3:
	    current=cTDC_3; break;
	  case 4:
	    current=cTDC_4; break;
	  case 5:
	    current=cTDC_5; break;
	  case 6:
	    current=cTDC_6; break;
	  case 7:
	    current=cTDC_7; break;
	  case 8:
	    current=cTDC_8; break;
	  case 9:
	    current=cTDC_9; break;
	  case 10:
	    current=cTDC_10; break;
	  case 11:
	    current=cTDC_11; break;
	  case 12:
	    current=cTDC_12; break;
	  case 13:
	    current=cTDC_13; break;
	  case 14:
	    current=cTDC_14; break;
	  case 15:
	    current=cTDC_15; break;
	  case 16:
	    current=cTDC_16; break;
	  default:
	    cout << "Canvas number error" << endl;
	  }
	  int entries = htmp->GetEntries();
	  float mean = htmp->GetMean(1);
	  float RMS = htmp->GetRMS(1);

	  ps = new TPaveText(0.5,0.8,0.9,0.9,"ndc");
	  ps->SetTextFont(42);
	  ps->SetTextSize(0.05);
	  ps->SetTextColor(kRed);
	  sprintf(tmpentry,"Total entries = %d",totentries); ps->AddText(tmpentry);



	  ps =(TPaveStats*)current->GetPrimitive("stats");
	  ps->SetName("mystats");
	  list = ps->GetListOfLines();
	  tconst = ps->GetLineWith("Entries"); list->Remove(tconst);
	  tconst = ps->GetLineWith("Mean"); list->Remove(tconst);
	  tconst = ps->GetLineWith("RMS"); list->Remove(tconst);

	  tmpentry.Form("Entries = %d",entries); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  tmpentry.Form("Mean = %.3f",mean); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  tmpentry.Form("RMS = %.3f",RMS); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);

	  htmp->SetStats(0);
	  current->Modified();
	};
	
	title.Form("run_%d_TDC_slot_%d_chan_%d_%d_adc_cut_%d.png",
		   run,tdc_slot,tdc_chan_start,tdc_chan_start+15,adc_cut);
	cTDC->Print(title);
	cTDC->cd(0);
	return cTDC;
}

//--------------------------------------------------------------------
//                        TDC vs ADC
//--------------------------------------------------------------------
// For interactive use: plot canvas divided in 4x4 TDC vs ADC plots of
// 16 channels starting from adc_chan_start, of given slot. No cuts.
// By default, tdc starting chan is calculated from adc starting chan;
// if passed a different value than -1 to tdc_chan_start, that value
// will be used instead. May be useful during testing of ADC/TDC
// channels. Fills global vector for noise freq. Optionally, a threshold on tdc for noise calculation purposes may be supplied
TCanvas *plot_tdc_vs_adc(int adc_slot, int adc_chan_start, int tdc_slot = -1, int tdc_chan_start = -1, int tdc_cut = 1350){
  TString cut, draw, draw1, title;
	title.Form("run_%d_TDC_vs_ADC",run);
	TCanvas *cTDCADC= new TCanvas("cTDCADC",title,850,850);
	int nbin_adc=50, nbin_tdc=50;
	int min_adc=-100, max_adc=1500;
	int min_tdc=-200, max_tdc=1500;
	TH2D *htmp=new TH2D("htmp","htmp",nbin_adc,min_adc,max_adc, nbin_tdc,min_tdc,max_tdc); // ok for TDC vs ADC

	//TString tmpentry;
	char tmpentry[25];
	TPad *current=0;
	//	TPaveStats *ps = 0;
	TList *list = 0;
	TText *tconst = 0;
	TLatex *myt = 0;
	MyStyle->SetStatX(0.9);
	MyStyle->SetStatY(0.6);

	cTDCADC->Divide(4,4); // 16 channels starting from adc_chan_start
	// log scale in number of counts, so that the pedestal events
	// do not push the others to the lowest part of the color
	// spectrum
	cTDCADC_1->SetLogz();
	cTDCADC_2->SetLogz();
	cTDCADC_3->SetLogz();
	cTDCADC_4->SetLogz();
	cTDCADC_5->SetLogz();
	cTDCADC_6->SetLogz();
	cTDCADC_7->SetLogz();
	cTDCADC_8->SetLogz();
	cTDCADC_9->SetLogz();
	cTDCADC_10->SetLogz();
	cTDCADC_11->SetLogz();
	cTDCADC_12->SetLogz();
	cTDCADC_13->SetLogz();
	cTDCADC_14->SetLogz();
	cTDCADC_15->SetLogz();
	cTDCADC_16->SetLogz();
	// 16 channels to plot
	for (int i=adc_chan_start; i<adc_chan_start+16; i++){
	  
	  if (tdc_chan_start == -1) tdc_chan_start=calc_tdc_chan(adc_slot,i);
	  if (tdc_slot == -1) tdc_slot = calc_tdc_slot(adc_slot,i);

	  cTDCADC->cd(i-adc_chan_start+1);
	  cout << "Doing plot " << i-adc_chan_start+1 << " / 16"<<endl;
	  
	  // TDC vs ADC:
	  htmp=new TH2D("htmp","htmp",nbin_adc,min_adc,max_adc, 
			nbin_tdc,min_tdc,max_tdc);
	  int itdc = tdc_chan_start + i-adc_chan_start;
	  draw.Form("tdct[%d][%d]:adc[%d][%d]>>htmp",tdc_slot,itdc,adc_slot,i);
	  htmp->GetXaxis()->SetTitle("ADC"); 
	  htmp->GetYaxis()->SetTitle("TDC");
	  
	  title.Form("Run %d TDC vs ADC slot %d chan %d",run,adc_slot, i);
	  htmp->SetTitle(title);

	  t->Draw(draw,"","colz");
	  cTDCADC->Update();

	  float Ncut = 0; float Nnocut = 0; 
	  int bincut = htmp->GetXaxis()->FindBin(tdc_cut);
	  float time = bincut *0.5e-9; // s 
	  // find the info for noise freq calculation
	  int binx1 = htmp->GetXaxis()->FindBin(min_adc);
	  int binx2 = htmp->GetXaxis()->FindBin(max_adc);
	  
	  //int biny1 = htmp->GetYaxis()->FindBin(min_tdc);
	  int biny1 = htmp->GetYaxis()->FindBin(0);
	  int biny2 = htmp->GetYaxis()->FindBin(max_tdc);
	  title.Form("Noise_tdc%d",itdc);
	  Noise_tdc[tdc_slot][itdc] = new TH1D(title,title, nbin_tdc,min_tdc,max_tdc);
	  for (int jj = 1; jj < nbin_tdc; jj++){
	    if (jj <= biny1) { // bin is below, or at,threshold set by biny1
	      int counts = 0;
	    } else {
	      int counts = htmp->Integral(binx1, binx2, biny1, jj);
	    }
	    Noise_tdc[tdc_slot][itdc]->SetBinContent(jj,counts);
	    Nnocut += counts;
	    if (jj < bincut) Ncut += counts;
	  }
	  if (Nnocut * time >0){
	    float noise = Ncut / Nnocut / time;
	  } else {
	    float noise = 0;
	  }

	  // Retrieve the stat box
	  switch(i-adc_chan_start+1){
	  case 1:
	    current=cTDCADC_1; break;
	  case 2:
	    current=cTDCADC_2; break;
	  case 3:
	    current=cTDCADC_3; break;
	  case 4:
	    current=cTDCADC_4; break;
	  case 5:
	    current=cTDCADC_5; break;
	  case 6:
	    current=cTDCADC_6; break;
	  case 7:
	    current=cTDCADC_7; break;
	  case 8:
	    current=cTDCADC_8; break;
	  case 9:
	    current=cTDCADC_9; break;
	  case 10:
	    current=cTDCADC_10; break;
	  case 11:
	    current=cTDCADC_11; break;
	  case 12:
	    current=cTDCADC_12; break;
	  case 13:
	    current=cTDCADC_13; break;
	  case 14:
	    current=cTDCADC_14; break;
	  case 15:
	    current=cTDCADC_15; break;
	  case 16:
	    current=cTDCADC_16; break;
	  default:
	    cout << "Canvas number error" << endl;
	  }
	  int entries = htmp->GetEntries();
	  float meanx = htmp->GetMean(1);
	  float meany = htmp->GetMean(2);
	  float RMSx = htmp->GetRMS(1);
	  float RMSy = htmp->GetRMS(2);

	  // ps =(TPaveStats*)current->GetPrimitive("stats");
	  // ps->SetName("mystats");
	  // list = ps->GetListOfLines();
	  // tconst = ps->GetLineWith("Entries"); list->Remove(tconst);
	  // tconst = ps->GetLineWith("Mean x"); list->Remove(tconst);
	  // tconst = ps->GetLineWith("Mean y"); list->Remove(tconst);
	  // tconst = ps->GetLineWith("RMS x"); list->Remove(tconst);
	  // tconst = ps->GetLineWith("RMS y"); list->Remove(tconst);

	  // tmpentry.Form("NoiseTDC<= %d",tdc_cut); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  // tmpentry.Form("Entries = %d",entries); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  // tmpentry.Form("mean x = %.3f",meanx); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  // tmpentry.Form("mean y = %.3f",meany); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  // tmpentry.Form("RMS x = %.3f",RMSx); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  // tmpentry.Form("RMS y = %.3f",RMSy); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	  // tmpentry.Form("Noise = %.3g kHz",noise*1e-3); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);

	  ps = new TPaveText(0.3,0.2,0.85,0.6,"ndc");
	  ps->SetTextFont(42);
	  ps->SetTextSize(0.05);
	  //ps->SetTextColor(kRed);
	  sprintf(tmpentry,"Entries = %d",entries); ps->AddText(tmpentry);
	  sprintf(tmpentry,"Mean x = %.3f",meanx); ps->AddText(tmpentry);
	  sprintf(tmpentry,"Mean y = %.3f",meany); ps->AddText(tmpentry);
	  sprintf(tmpentry,"RMS x = %.3f",RMSx); ps->AddText(tmpentry);
	  sprintf(tmpentry,"RMS y = %.3f",RMSy); ps->AddText(tmpentry);
	  sprintf(tmpentry,"Noise cut 0<TDC< %d",tdc_cut); ps->AddText(tmpentry);
	  sprintf(tmpentry,"      over: 0<TDC< %d",max_tdc); ps->AddText(tmpentry);
	  sprintf(tmpentry,"Noise = %.3f kHz",noise/1000); ps->AddText(tmpentry);
	  ps->Draw();

	  htmp->SetStats(0);
	  current->Modified();
	};
	
	title.Form("run_%d_TDC_vs_ADC_slot_%d_chan_%d_%d.png",run,adc_slot,
		   adc_chan_start,adc_chan_start+15);
	cTDCADC->Print(title);
	cTDCADC->cd(0);
	
	return cTDCADC;
}

//--------------------------------------------------------------------
//                        peak (vertical tracks)
//--------------------------------------------------------------------
// Interactive use:
// plots the peak for vertical tracks for 16 chans starting from given
// slot/channel imposing the adjacent channels have counts < value_lt,
// and the chosen channel at least value_gt counts. Fill the global vectors
TCanvas *plot_peak(int slot, int chan_start, int value_lt=5, int value_gt=150){
  TCanvas *cpeak = new TCanvas("cpeak","cpeak",850,850); 
  cpeak->Divide(4,4);
  int chan;

  TPad *current=0;
  TPaveText *ps = 0;
  TText *tconst = 0;
  TLatex *myt = 0;

  cpeak_1->SetLogy();cpeak_2->SetLogy();cpeak_3->SetLogy();cpeak_4->SetLogy();
  cpeak_5->SetLogy();cpeak_6->SetLogy();cpeak_7->SetLogy();cpeak_8->SetLogy();
  cpeak_9->SetLogy();cpeak_10->SetLogy();cpeak_11->SetLogy();cpeak_12->SetLogy();
  cpeak_13->SetLogy();cpeak_14->SetLogy();cpeak_15->SetLogy();cpeak_16->SetLogy();

  TString title, draw, draw1, cut;
  for (int ichan = 0; ichan < 16; ichan++){
    cpeak->Update();
    cpeak->cd(ichan+1);
    chan = chan_start+ichan;
    title.Form("Run %d slot %d channel %d",run,slot,chan);
    draw.Form("adc[%d][%d]>>htemp",slot,chan);
    draw1.Form("adc[%d][%d]>>htemp1",slot,chan);
    TH1D *htemp = new TH1D("htemp",title,150,-100,500);
    TH1D *htemp1 = new TH1D("htemp1",title,150,-100,500);
    
    // check for boundary conditions
    int tmpchan = chan % 16; // chan reduced in [0,15]
    switch(tmpchan){
    case 0: // first chan, so no previous one
      cut.Form("adc[%d][%d] < %d && adc[%d][%d]>%d", 
	       slot,chan+1,value_lt, slot,chan,value_gt);
      break;
    case 15: // last chan, so no next one
    cut.Form("adc[%d][%d] < %d && adc[%d][%d]>%d", 
	     slot,chan-1,value_lt, slot,chan,value_gt);
    break;
    default: // chan not on edge
      cut.Form("adc[%d][%d] < %d && adc[%d][%d] < %d && adc[%d][%d]>%d",
	       slot,chan-1,value_lt, slot,chan+1,value_lt, slot,chan,value_gt);
    }
    t->Draw(draw,"");
    t->Draw(draw1,cut);
    
    float counts, mean, rms;
    mean = htemp1->GetMean();
    rms = htemp1->GetRMS();
    counts = htemp1->GetMaximum();
    cout << "before: "<< mean << " " << rms << endl;
    TF1 *fpeak = new TF1("fpeak","gaus",mean-peak_coeff_range*rms,mean+peak_coeff_range*rms);
    fpeak->SetParameter(0,counts);
    fpeak->SetParameter(1,mean);
    fpeak->SetParameter(2,rms);
    fpeak->SetParLimits(0,0,5*counts);
    fpeak->SetParLimits(1,mean-peak_coeff_mean_min*rms,mean+peak_coeff_mean_max*rms);
    fpeak->SetParLimits(2,peak_coeff_rms_min*rms,peak_coeff_rms_max*rms);
    htemp1->Fit("fpeak","B");
    counts = fpeak->GetParameter(0);
    mean = fpeak->GetParameter(1);
    rms = fpeak->GetParameter(2);
    float chi2 = fpeak->GetChisquare();
    int ndf = fpeak->GetNDF();

    cout << "after: " << mean << " " << rms << endl;

    Peak_mean[slot][chan] = mean;
    Peak_rms[slot][chan] = rms;

    htemp->SetLineColor(kRed);
    htemp1->SetLineColor(kBlue);
    fpeak->SetLineColor(kGreen);

    htemp->Draw();
    htemp1->Draw("same");
    fpeak->Draw("same");
 
    // Retrieve the stat box
    switch(ichan+1){
    case 1:
      current=cpeak_1; break;
    case 2:
      current=cpeak_2; break;
    case 3:
      current=cpeak_3; break;
    case 4:
      current=cpeak_4; break;
    case 5:
      current=cpeak_5; break;
    case 6:
      current=cpeak_6; break;
    case 7:
      current=cpeak_7; break;
    case 8:
      current=cpeak_8; break;
    case 9:
      current=cpeak_9; break;
    case 10:
      current=cpeak_10; break;
    case 11:
      current=cpeak_11; break;
    case 12:
      current=cpeak_12; break;
    case 13:
      current=cpeak_13; break;
    case 14:
      current=cpeak_14; break;
    case 15:
      current=cpeak_15; break;
    case 16:
      current=cpeak_16; break;
    default:
      cout << "Canvas number error" << endl;
    }
    int totentries = htemp->GetEntries();
    int cutentries = htemp1->GetEntries();
    char tmpentry[25];
  
    ps = new TPaveText(0.5,0.8,0.9,0.9,"ndc");
    ps->SetTextFont(42);
    ps->SetTextSize(0.05);
    ps->SetTextColor(kRed);
    sprintf(tmpentry,"Total entries = %d",totentries); ps->AddText(tmpentry);
    
    ps1 = new TPaveText(0.5,0.6,0.9,0.8,"ndc");
    ps1->SetTextFont(42);
    ps1->SetTextSize(0.05);
    ps1->SetTextColor(kBlue);
    sprintf(tmpentry,"Cut: neighbours < %d",value_lt); ps1->AddText(tmpentry);
    sprintf(tmpentry,"Cut: this chan > %d",value_gt); ps1->AddText(tmpentry);
    sprintf(tmpentry,"Cut: entries %d",cutentries); ps1->AddText(tmpentry);
    sprintf(tmpentry,"%s"," "); ps1->AddText(tmpentry);
    
    ps2 = new TPaveText(0.5,0.38,0.9,0.63,"ndc");
    ps2->SetTextFont(42);
    ps2->SetTextSize(0.05);
    ps2->SetTextColor(kGreen);
    sprintf(tmpentry,"Peak mean = %d",mean); ps2->AddText(tmpentry);
    sprintf(tmpentry,"Peak rms = %d",rms); ps2->AddText(tmpentry);
    if (rms>0) {
      float phe = (mean/rms)**2;
      float gain = mean/phe*adc_charge/e;
    } else {
      float phe = 0;
      float gain = 0;
    }
    sprintf(tmpentry,"chi2/NDF = %.2f / %d",chi2,ndf); ps2->AddText(tmpentry);
    sprintf(tmpentry,"Nphe = %.2f",phe); ps2->AddText(tmpentry);
    sprintf(tmpentry,"gain = %.2g",gain); ps2->AddText(tmpentry);
    
    ps->Draw("same");
    ps1->Draw("same");
    ps2->Draw("same");

    htemp->SetStats(0);
    htemp1->SetStats(0);
    current->Modified();
  }
  cpeak->cd(0);

  title.Form("run_%d_peaks_slot_%d_chans_%d_%d_lt_%d_gt_%d.png",
	     run,slot,chan_start,chan_start+15,value_lt,value_gt);
  cpeak->Print(title);
  return cpeak;
}

//--------------------------------------------------------------------
//                        eta (efficiency)
//--------------------------------------------------------------------
// the fit function
Double_t _eta_func(Double_t *x, Double_t *par){
  Double_t A = x[0];
  Double_t Ath = par[0];
  Double_t Aw = par[1];
  static const Double_t pi = 4*atan(1);
  float diff;
  if (Aw > 0.) {
    diff = ((A - Ath) / Aw);

    // version atan
    //return atan( diff+diff**3 /5  ) / pi + 0.5; // between 0 and 1
    
    // version tanh
    // normalizing a*tanh(diff)+b in [0,1], we find a=b=1/2. The
    // resulting expression can be simplified in 1/(1+exp(-2diff))
    return 1./(1.+exp(-2*diff));
  } else {
    return -1;
  }
};

// For interactive use:
// plot 16 eta (calculated from histos) and its fit for a given slot
// and from a channel. The fit results are stored in global arrays
// Athresholds and Awidths
TCanvas *plot_eta(int adc_slot, int adc_chan_start, int tdc_cut=1350, int tdc_slot = -1, int tdc_chan_start = -1){
  const int n=200;
  const int bin_amplitude = 5;
  last_tdc_cut = tdc_cut;
  TString draw, cut, title;
  float ratio;
  TCanvas *ceta_multi = new TCanvas("ceta_multi","ceta_multi",850,850);
  ceta_multi->Divide(4,4);
  const int nbin = n/bin_amplitude;
  int Ncut, Nnocut;
  TH1D eta_raw_adctot = TH1D("eta_raw_adctot","eta_raw_adctot",nbin, 0, n);
  TH1D eta_raw_adccut = TH1D("eta_raw_adccut","eta_raw_adccut",nbin, 0, n);
  eta_raw_adctot.Sumw2();
  eta_raw_adccut.Sumw2();

  MyStyle->SetStatX(0.98);
  MyStyle->SetStatY(0.3);

  for (int chan = adc_chan_start; chan < adc_chan_start+16; chan++){
    ceta_multi->Update();

    if (tdc_chan_start == -1) tdc_chan_start=calc_tdc_chan(adc_slot,chan);
    if (tdc_slot == -1) tdc_slot = calc_tdc_slot(adc_slot,chan);

    ceta_multi->cd(chan-adc_chan_start+1);
    int itdc = tdc_chan_start + chan - adc_chan_start;

    draw.Form("adc[%d][%d]>>eta_raw_adctot",adc_slot,chan);
    t->Draw(draw,"");
    draw.Form("adc[%d][%d]>>eta_raw_adccut", adc_slot,chan);
    cut.Form("tdct[%d][%d]>%d",tdc_slot, itdc, tdc_cut);
    t->Draw(draw,cut);
    TEfficiency *pEff = 0;
    if (TEfficiency::CheckConsistency(eta_raw_adccut, eta_raw_adctot,"w") ){
      pEff = new TEfficiency(eta_raw_adccut, eta_raw_adctot);
    }
    float Ath_estimate = 70.;
    float Aw_estimate = 20.;
    //cout<<"Estimates: Ath "<< Ath_estimate <<" Aw "<< Aw_estimate <<endl;
    TF1 *eta = new TF1("eta",_eta_func,40.,140, 2);
    eta->SetParameter(0, Ath_estimate);
    eta->SetParameter(1, Aw_estimate);
    //    eta->SetParLimits(0, 20, Ath_estimate+50);
    //    eta->SetParLimits(1, 20, 1000);
    pEff->Draw("AP");
    pEff->Fit(eta,"R"); // preliminary fit
    float Ath = eta->GetParameter(0);
    float Aw = eta->GetParameter(1);
    float Ath_err = eta->GetParError(0);
    float Aw_err = eta->GetParError(1);
    float chi2 = eta->GetChisquare();
    // titles
    TH1D *htemp=(TH1D*)gPad->GetPrimitive("eta_raw_adctot_clone");
    title.Form("Run %d slot %d channel %d TDC > %d", run,adc_slot,chan,tdc_cut);
    htemp->SetTitle(title);
    delete eta;
    float Athmin = Ath-50<0?0:Ath-50;
    TF1 *eta = new TF1("eta",_eta_func,Athmin, Athmin + 100, 2);
    eta->SetLineColor(kRed);
    eta->SetParNames("Ath","Aw");
    eta->SetParameter(0, Ath);
    eta->SetParameter(1, Aw);
    eta->SetParLimits(0, Athmin, Athmin+100);
    //    eta->SetParLimits(1, 15, 50);
    pEff->Fit(eta,"R"); // better fit
    // update values
    Ath = eta->GetParameter(0);
    Aw = eta->GetParameter(1);
    Ath_err = eta->GetParError(0);
    Aw_err = eta->GetParError(1);
    chi2 = eta->GetChisquare();
    // fill global arrays
    Athresholds[adc_slot][chan] = Ath;
    Awidths[adc_slot][chan] = Aw;
    cout << "slot "<< adc_slot <<" chan "<< chan <<" Saved Ath = " << Ath << " Aw = " << Aw <<" in arrays Athresholds and Awidths"<<endl;
    //   ceta_multi->Modified();
  }

  title.Form("run_%d_eta_slot_%d_chans_%d_%d_cut_%d.png",
	     run,adc_slot,adc_chan_start,adc_chan_start+15,tdc_cut);
  ceta_multi->cd(0);
  ceta_multi->Print(title);
  return ceta_multi;
}

//====================================================================
//                        Ancillary routines
//     (normally there is no need to call them interactively)
//====================================================================

// Calculate all eta fits for slot, between chan_min and chan_max
// (inclusive), in order to populate the arrays. Used in save_eta()
void fit_eta_in_slot(int adc_slot=0, int adc_chan_min=0, int adc_chan_max=63,int tdc_cut=1350, int tdc_slot = -1, int tdc_chan = -1){
  last_tdc_cut = tdc_cut;
  for (int chan = adc_chan_min; chan < adc_chan_max+1; chan++){
    plot_eta(adc_slot, chan, tdc_cut, tdc_slot, tdc_chan);
  }
}

// find the mean of the peak for a given slot/channel from the global vectors
float peak_mean(int slot, int chan){
  return Peak_mean[slot][chan];
}
// Interactive use:
// find the RMS of the peak for a given slot/channel from the global vectors
float peak_rms(int slot, int chan){
  return Peak_rms[slot][chan];
}

// Calculate the number of phe for a given adc slot, chain and cut
float phe(int adc_slot, int adc_chan, int value_lt=5, int value_gt=100){
  float rms = peak_rms(adc_slot, adc_chan, value_lt, value_gt);
  if (rms > 0) {
    float mean = peak_mean(adc_slot,adc_chan,value_lt, value_gt);
    return (mean/rms)**2;
  } else {
    return 0;
  }
}

// calculate gain for a given adc slot, chan, cut
float gain(int adc_slot, int adc_chan, int value_lt=5, int value_gt=100){
  const float adc_charge = 50*1e-15; // 50 fC, corresponding to an ADC chan
  const float e = 1.6e-19; // C, electron charge

  float A = peak_mean(adc_slot,adc_chan,value_lt, value_gt);
  float phe = phe(adc_slot, adc_chan, value_lt, value_gt);

  float Qphe = A/phe*adc_charge; // total charge in phe

  return Qphe/e; // gain
}

// find the tdc slot corresponding to a given adc slot/channel
int calc_tdc_slot(int adc_slot, int adc_chan){
  int diff = NUMCHANT - NUMCHANA; // = 32 = NUMCHANA/2, so 3 ADC fill 2 TDC
  if (adc_slot == 0) return 0; // start filling tdc slot 0
  if (adc_slot == 1){
    if (adc_chan < diff) {
      return 0; // end filling tdc slot 0
    } else {
      return 1; // start filling tdc slot 1
    }
  }
  if (adc_slot == 2) return 1; // end filling tdc slot 1
  if (adc_slot == 3) return 2; // start filling tdc slot 2 
}
// find the tdc chan corresponding to a given adc slot/channel
int calc_tdc_chan(int adc_slot, int adc_chan){
  int diff = NUMCHANT - NUMCHANA; // = 32 = NUMCHANA/2, so 3 ADC fill 2 TDC
  if (adc_slot == 0) return adc_chan; // start filling tdc slot 0
  if (adc_slot == 1){
    if (adc_chan < diff) {
      return NUMCHANA + adc_chan; // end filling tdc slot 0
    } else {
      return adc_chan - diff; // start filling tdc slot 1
    }
  }
  if (adc_slot == 2) return adc_chan + diff; // end filling tdc slot 1
  if (adc_slot == 3) return adc_chan; // start filling tdc slot 2
}
// find the adc slot corresponding to a given tdc slot/chan
int calc_adc_slot(int tdc_slot, int tdc_chan){
  int diff = NUMCHANT - NUMCHANA; // = 32 = half ADC module
  if (tdc_slot == 0){
    if (tdc_chan < NUMCHANA) {
      return 0; // adc 0 is fully contained in tdc 0
    } else {
      return 1; // adc 1 is half contained in tdc 0
    }
  }
  if (tdc_slot == 1){
    if (tdc_chan < diff) {
      return 1; // this ends adc 1
    } else {
      return 2; // adc 2 is fully in tdc 1)
    }
  }
  if (tdc_slot == 2) return 3; // tdc 2 start === adc 3 start
}
int calc_adc_chan(int tdc_slot, int tdc_chan){
  int diff = NUMCHANT - NUMCHANA;
  if (tdc_slot == 0){
    if (tdc_chan < NUMCHANA) {
      return tdc_chan; // adc 0 is fully contained in tdc 0
    } else {
      return tdc_chan - NUMCHANA; // adc 1 is half contained in tdc 0
    }
  }
  if (tdc_slot == 1){
    if (tdc_chan < diff) {
      return tdc_chan + diff; // this ends adc 1
    } else {
      return tdc_chan - diff; // adc 2 is fully in tdc 1)
    }
  }
  if (tdc_slot == 2) return tdc_chan; // tdc 2 start === adc 3 start
}

calculate the noise frequency (in Hz) for given tdc slot, chan, cut
float noise_freq(int tdc_slot, int tdc_chan, int tdc_cut=last_tdc_cut){
  TString draw,zerocut, noisecut;
  float ntot, ncut;
  float time, freq;

  // TCanvas cnoise = TCanvas("cnoise","cnoise");
  // draw.Form("tdct[%d][%d]",tdc_slot,tdc_chan);
  // zerocut.Form("tdct[%d][%d]>0",tdc_slot,tdc_chan);
  // noisecut.Form("tdct[%d][%d]>0 && tdct[%d][%d]<%d",
  // 		tdc_slot,tdc_chan,tdc_slot,tdc_chan,tdc_cut);

  // ntot = t->Draw(draw,zerocut);
  // ncut = t->Draw(draw,noisecut);

  int nbins = Noise_tdc[tdc_slot][tdc_chan]->GetNbinsX();
  int bincut = Noise_tdc[tdc_slot][tdc_chan]->FindBin(tdc_cut);
  ntot = Noise_tdc[tdc_slot][tdc_chan]->Integral(0,nbins);
  ncut = Noise_tdc[tdc_slot][tdc_chan]->Integral(0,bincut);

  time = bincut * 0.5*1e-9; // in seconds, because each tdc chan = 0.5 ns

  cout<< ncut <<" events / "<< ntot <<" events / "<< time <<" s"<<endl;

  if (ntot * time > 0) {
    freq = ncut / ntot / time; // Hertz
  } else {
    freq = -1;
  }
  return freq;
}
  
//====================================================================
//                          Misc section
//====================================================================



//====================================================================
//                     Single plot versions section
//====================================================================

TCanvas *plot_adc1(int adc_slot, int adc_chan, int tdc_cut=1350, int tdc_slot = -1, int tdc_chan = -1){
	TString cut, draw, draw1, title;
	title.Form("run_%d_ADC",run);
	TCanvas *cADC1= new TCanvas("cADC1",title);
	int nbin=100;
	int min=-100, max=500;
	TH1D *htmp=new TH1D("htmp","htmp",nbin,min,max);
	TH1D *htmp1=new TH1D("htmp1","htmp1",nbin,min,max);
	TPad *current=0;
	TPaveStats *ps = 0;
	TList *list = 0;
	TText *tconst = 0;
	TLatex *myt = 0;
	TString tmpentry;
	MyStyle->SetStatX(0.9);
	MyStyle->SetStatY(0.6);
	last_tdc_cut = tdc_cut;

	// log scale in number of counts, so that the pedestal events
	// do not push the others to the lowest part of the color
	// spectrum
	cADC1->SetLogy();

	if (tdc_chan == -1) tdc_chan_start=calc_tdc_chan(adc_slot,adc_chan);
	if (tdc_slot == -1) tdc_slot = calc_tdc_slot(adc_slot,chan);

	draw.Form("adc[%d][%d]>>htmp",adc_slot,adc_chan);
	draw1.Form("adc[%d][%d]>>htmp1",adc_slot,adc_chan);
	cut.Form("tdct[%d][%d]>%d",tdc_slot,tdc_chan,tdc_cut);
	htmp->SetLineColor(kBlue);
	htmp1->SetLineColor(kRed);
	title.Form("Run %d ADC slot %d chan %d tdc > %d",run,adc_slot, adc_chan, tdc_cut);
	htmp->SetTitle(title);
	htmp1->SetTitle(title);
	  
	t->Draw(draw,"");
	t->Draw(draw1,cut,"same");

	//cADC1->Update();
	int entries = htmp->GetEntries();
	float mean = htmp->GetMean(1);
	float RMS = htmp->GetRMS(1);

	ps =(TPaveStats*)cADC1->GetPrimitive("stats");
	ps->SetName("mystats");
	list = ps->GetListOfLines();
	tconst = ps->GetLineWith("Entries"); list->Remove(tconst);
	tconst = ps->GetLineWith("Mean"); list->Remove(tconst);
	tconst = ps->GetLineWith("RMS"); list->Remove(tconst);

	tmpentry.Form("Entries = %d",entries); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	tmpentry.Form("Mean = %.3f",mean); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	tmpentry.Form("RMS = %.3f",RMS); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);

	htmp->SetStats(0);
	cADC1->Modified();
	
	title.Form("run_%d_ADC_slot_%d_chan_%d_tdc_cut_%d.png",
		   run,adc_slot,adc_chan,tdc_cut);
	// cADC1->Print(title);
	return cADC1;
}

TCanvas *plot_tdc1(int tdc_slot, int tdc_chan, int adc_cut, int adc_slot = -1, int adc_chan = -1){
	TString cut, draw, draw1, title;
	title.Form("run_%d_TDC",run);
	TCanvas *cTDC1= new TCanvas("cTDC1",title);
	int nbin=100;
	int min=-100, max=1500;
	TH1D *htmp=new TH1D("htmp","htmp",nbin,min,max);
	TH1D *htmp1=new TH1D("htmp1","htmp1",nbin,min,max);
	TString tmpentry;
	TPad *current=0;
	TPaveStats *ps = 0;
	TList *list = 0;
	TText *tconst = 0;
	TLatex *myt = 0;
	MyStyle->SetStatX(0.6);
	MyStyle->SetStatY(0.6);

	// log scale in number of counts, so that the pedestal events
	// do not push the others to the lowest part of the color
	// spectrum
	cTDC1->SetLogy();

	if (adc_chan == -1) adc_chan = calc_adc_chan(tdc_slot,tdc_chan);
	if (adc_slot == -1) adc_slot = adc_adc_slot(tdc_slot, tdc_chan);

	draw.Form("tdct[%d][%d]>>htmp",tdc_slot,tdc_chan);
	draw1.Form("tdct[%d][%d]>>htmp1",tdc_slot,tdc_chan);
	cut.Form("adc[%d][%d]>%d",adc_slot,adc_chan,adc_cut);
	htmp->SetLineColor(kBlue);
	htmp1->SetLineColor(kRed);
	title.Form("Run %d TDC slot %d chan %d adc > %d",run,tdc_slot, i, adc_cut);
	htmp->SetTitle(title);
	htmp1->SetTitle(title);
	
	t->Draw(draw);
	t->Draw(draw1,cut,"same");
	
	//cTDC1->Update();
	int entries = htmp->GetEntries();
	float mean = htmp->GetMean(1);
	float RMS = htmp->GetRMS(1);
	
	// ps =(TPaveStats*)cTDC1->GetPrimitive("stats");
	// ps->SetName("mystats");
	// list = ps->GetListOfLines();
	// tconst = ps->GetLineWith("Entries"); list->Remove(tconst);
	// tconst = ps->GetLineWith("Mean"); list->Remove(tconst);
	// tconst = ps->GetLineWith("RMS"); list->Remove(tconst);
	
	// tmpentry.Form("Entries = %d",entries); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	// tmpentry.Form("Mean = %.3f",mean); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
	// tmpentry.Form("RMS = %.3f",RMS); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);

	ps = new TPaveText(0.3,0.2,0.85,0.6,"ndc");
	ps->SetTextFont(42);
	ps->SetTextSize(0.05);
	//ps->SetTextColor(kRed);
	sprintf(tmpentry,"Entries = %d",entries); ps->AddText(tmpentry);
	sprintf(tmpentry,"Mean x = %.3f",meanx); ps->AddText(tmpentry);
	sprintf(tmpentry,"Mean y = %.3f",meany); ps->AddText(tmpentry);
	sprintf(tmpentry,"RMS x = %.3f",RMSx); ps->AddText(tmpentry);
	sprintf(tmpentry,"RMS y = %.3f",RMSy); ps->AddText(tmpentry);
	sprintf(tmpentry,"Noise cut 0<TDC<= %d",tdc_cut); ps->AddText(tmpentry);
	sprintf(tmpentry,"Noise = %.3f kHz",noise/1000); ps->AddText(tmpentry);
	ps->Draw();

	htmp->SetStats(0);
	cTDC1->Modified();
	
	title.Form("run_%d_TDC_slot_%d_chan_%d_adc_cut_%d.png",
		   run,tdc_slot,tdc_chan,adc_cut);
	// cTDC1->Print(title);
	return cTDC1;
}

//Optionally, a threshold on tdc for noise calculation purposes may be supplied
TCanvas *plot_tdc_vs_adc1(int adc_slot, int adc_chan, int tdc_slot = -1, int tdc_chan = -1, int tdc_cut){
  TString cut, draw, draw1, title;
  title.Form("run_%d_TDC_vs_ADC",run);
  TCanvas *cTDCADC1= new TCanvas("cTDCADC1",title);
  int nbin_adc=50, nbin_tdc=50;
  int min_adc=-100, max_adc=1500;
  int min_tdc=-200, max_tdc=1500;

  TH2D *htmp=new TH2D("htmp","htmp",nbin_adc,min_adc,max_adc, nbin_tdc,min_tdc,max_tdc); // TDC vs ADC

  Noise_tdc[tdc_slot][tdc_chan] = new TH1D(title,title, max_tdc,min_tdc,max_tdc);
  
  //TString tmpentry;
  char tmpentry[25];
  TPad *current=0;
  //TPaveStats *ps = 0;
  TList *list = 0;
  TText *tconst = 0;
  TLatex *myt = 0;
  MyStyle->SetStatX(0.9);
  MyStyle->SetStatY(0.6);
  
  // log scale in number of counts, so that the pedestal events
  // do not push the others to the lowest part of the color
  // spectrum
  cTDCADC1->SetLogz();
	  
  if (tdc_chan == -1) tdc_chan=calc_tdc_chan(adc_slot,adc_chan);
  if (tdc_slot == -1) tdc_slot = calc_tdc_slot(adc_slot,adc_chan);

  htmp=new TH2D("htmp","htmp",nbin_adc,min_adc,max_adc, 
		nbin_tdc,min_tdc,max_tdc);
  draw.Form("tdct[%d][%d]:adc[%d][%d]>>htmp",tdc_slot,tdc_chan,adc_slot,adc_chan);
  htmp->GetXaxis()->SetTitle("ADC"); 
  htmp->GetYaxis()->SetTitle("TDC");
	  
  title.Form("Run %d TDC vs ADC slot %d chan %d",run,adc_slot, adc_chan);
  htmp->SetTitle(title);
  
  t->Draw(draw,"","colz");
  cTDCADC1->Update();


  float Ncut = 0; float Nnocut = 0; 
  int bincut = htmp->GetXaxis()->FindBin(tdc_cut); // default tdc_cut
  float time = bincut *0.5e-9; // s 
  // find the info for noise freq calculation
  int binx1 = htmp->GetXaxis()->FindBin(min_adc);
  int binx2 = htmp->GetXaxis()->FindBin(max_adc);
  //int biny1 = htmp->GetYaxis()->FindBin(min_tdc);
  int biny1 = htmp->GetYaxis()->FindBin(0);
  int biny2 = htmp->GetYaxis()->FindBin(max_tdc);
  title.Form("Noise_tdc%d",tdc_chan);
  Noise_tdc[tdc_slot][tdc_chan] = new TH1D(title,title, nbin_tdc,min_tdc,max_tdc);
  for (int jj = 1; jj < nbin_tdc; jj++){
    if (jj <= biny1) { // bin is below, or at,threshold set by biny1
      int counts = 0;
    } else {
      int counts = htmp->Integral(binx1, binx2, biny1, jj);
    }
    Noise_tdc[tdc_slot][tdc_chan]->SetBinContent(jj,counts);
    Nnocut += counts;
    if (jj < bincut) Ncut += counts;
  }
  if (Nnocut * time >0){
    float noise = Ncut / Nnocut / time;
  } else {
    float noise = 0;
  }

  int entries = htmp->GetEntries();
  float meanx = htmp->GetMean(1);
  float meany = htmp->GetMean(2);
  float RMSx = htmp->GetRMS(1);
  float RMSy = htmp->GetRMS(2);

  // ps =(TPaveStats*)cTDCADC1->GetPrimitive("stats");
  // ps->SetName("mystats");
  // list = ps->GetListOfLines();
  // tconst = ps->GetLineWith("Entries"); list->Remove(tconst);
  // tconst = ps->GetLineWith("Mean x"); list->Remove(tconst);
  // tconst = ps->GetLineWith("Mean y"); list->Remove(tconst);
  // tconst = ps->GetLineWith("RMS x"); list->Remove(tconst);
  // tconst = ps->GetLineWith("RMS y"); list->Remove(tconst);
  
  // tmpentry.Form("Entries = %d",entries); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
  // tmpentry.Form("mean x = %.3f",meanx); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
  // tmpentry.Form("mean y = %.3f",meany); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
  // tmpentry.Form("RMS x = %.3f",RMSx); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
  // tmpentry.Form("RMS y = %.3f",RMSy); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
  // tmpentry.Form("Noise = %.3g MHz",noise*1e-6); myt = new TLatex(0,0,tmpentry); myt->SetTextFont(42); myt->SetTextSize(0.04); list->Add(myt);
  ps = new TPaveText(0.3,0.2,0.85,0.6,"ndc");
  ps->SetTextFont(42);
  ps->SetTextSize(0.05);
  //ps->SetTextColor(kRed);
  sprintf(tmpentry,"Entries = %d",entries); ps->AddText(tmpentry);
  sprintf(tmpentry,"Mean x = %.3f",meanx); ps->AddText(tmpentry);
  sprintf(tmpentry,"Mean y = %.3f",meany); ps->AddText(tmpentry);
  sprintf(tmpentry,"RMS x = %.3f",RMSx); ps->AddText(tmpentry);
  sprintf(tmpentry,"RMS y = %.3f",RMSy); ps->AddText(tmpentry);
  sprintf(tmpentry,"Noise cut: 0<TDC< %d",tdc_cut); ps->AddText(tmpentry);
  sprintf(tmpentry,"      over: 0<TDC< %d",max_tdc); ps->AddText(tmpentry);
  sprintf(tmpentry,"Noise = %.3f kHz",noise/1000); ps->AddText(tmpentry);
  ps->Draw();  
	  
  htmp->SetStats(0);
  cTDCADC1->Modified();

  title.Form("run_%d_TDC_vs_ADC_slot_%d_chan_%d.png",run,adc_slot, adc_chan);
  //	cTDCADC1->Print(title);
	
  return cTDCADC1;
}

// Interactive use:
// plots the peak for vertical tracks for given slot/channel imposing
// the adjacent channels have counts < value_lt, and the chosen
// channel at least value_gt counts. Fill the global vectors
TCanvas *plot_peak1(int slot, int chan, int value_lt=5, int value_gt=150){
  TCanvas *cpeak1 = new TCanvas("cpeak1","cpeak1"); 

  cpeak1->SetLogy();

  TString title, draw, draw1, cut;

  title.Form("Run %d slot %d channel %d",run,slot,chan);
  draw.Form("adc[%d][%d]>>htemp",slot,chan);
  draw1.Form("adc[%d][%d]>>htemp1",slot,chan);
  TH1D *htemp = new TH1D("htemp",title,150,-100,500);
  TH1D *htemp1 = new TH1D("htemp1",title,150,-100,500);

  // check for boundary conditions
  int tmpchan = chan % 16; // chan reduced in [0,15]
  switch(tmpchan){
  case 0: // first chan, so no previous one
    cut.Form("adc[%d][%d] < %d && adc[%d][%d]>%d", 
	     slot,chan+1,value_lt, slot,chan,value_gt);
    break;
  case 15: // last chan, so no next one
    cut.Form("adc[%d][%d] < %d && adc[%d][%d]>%d", 
	     slot,chan-1,value_lt, slot,chan,value_gt);
    break;
  default: // chan not on edge
    cut.Form("adc[%d][%d] < %d && adc[%d][%d] < %d && adc[%d][%d]>%d",
	     slot,chan-1,value_lt, slot,chan+1,value_lt, slot,chan,value_gt);
  }
  t->Draw(draw,"");
  t->Draw(draw1,cut);
    
  float counts, mean, rms;
  mean = htemp1->GetMean();
  rms = htemp1->GetRMS();
  counts = htemp1->GetMaximum();
  cout << "before: "<< mean << " " << rms << endl;
  TF1 *fpeak = new TF1("fpeak","gaus",mean-peak_coeff_range*rms,mean+peak_coeff_range*rms);
  fpeak->SetParameter(0,counts);
  fpeak->SetParameter(1,mean);
  fpeak->SetParameter(2,rms);
  fpeak->SetParLimits(0,0,5*counts);
  fpeak->SetParLimits(1,mean-peak_coeff_mean_min*rms,mean+peak_coeff_mean_max*rms);
  fpeak->SetParLimits(2,peak_coeff_rms_min*rms,peak_coeff_rms_max*rms);
  htemp1->Fit("fpeak","B");
  counts = fpeak->GetParameter(0);
  mean = fpeak->GetParameter(1);
  rms = fpeak->GetParameter(2);
  float chi2 = fpeak->GetChisquare();
  int ndf = fpeak->GetNDF();

  cout << "after: " << mean << " " << rms << endl;

  Peak_mean[slot][chan] = mean;
  Peak_rms[slot][chan] = rms;
  
  htemp->SetLineColor(kRed);
  htemp1->SetLineColor(kBlue);
  fpeak->SetLineColor(kGreen);

  htemp->Draw();
  htemp1->Draw("same");
  fpeak->Draw("same");
 
  int totentries = htemp->GetEntries();
  int cutentries = htemp1->GetEntries();
  char tmpentry[25];
 
  ps = new TPaveText(0.5,0.8,0.9,0.9,"ndc");
  ps->SetTextFont(42);
  ps->SetTextSize(0.05);
  ps->SetTextColor(kRed);
  sprintf(tmpentry,"Total entries = %d",totentries); ps->AddText(tmpentry);

  ps1 = new TPaveText(0.5,0.6,0.9,0.8,"ndc");
  ps1->SetTextFont(42);
  ps1->SetTextSize(0.05);
  ps1->SetTextColor(kBlue);
  sprintf(tmpentry,"Cut: neighbours < %d",value_lt); ps1->AddText(tmpentry);
  sprintf(tmpentry,"Cut: this chan > %d",value_gt); ps1->AddText(tmpentry);
  sprintf(tmpentry,"Cut: entries %d",cutentries); ps1->AddText(tmpentry);
  sprintf(tmpentry,"%s"," "); ps1->AddText(tmpentry);

  ps2 = new TPaveText(0.5,0.38,0.9,0.63,"ndc");
  ps2->SetTextFont(42);
  ps2->SetTextSize(0.05);
  ps2->SetTextColor(kGreen);
  sprintf(tmpentry,"Peak mean = %d",mean); ps2->AddText(tmpentry);
  sprintf(tmpentry,"Peak rms = %d",rms); ps2->AddText(tmpentry);
  if (rms>0) {
    float phe = (mean/rms)**2;
    float gain = mean/phe*adc_charge/e;
  } else {
    float phe = 0;
  }
  sprintf(tmpentry,"chi2 / ndf = %.2f / %d",chi2,ndf); ps2->AddText(tmpentry);
  sprintf(tmpentry,"Nphe = %.2f",phe); ps2->AddText(tmpentry);
  sprintf(tmpentry,"gain = %.2g",gain); ps2->AddText(tmpentry);
  
  ps->Draw("same");
  ps1->Draw("same");
  ps2->Draw("same");
  htemp->SetStats(0);
  htemp1->SetStats(0);
  
  return cpeak1;
}

For interactive use:
plot eta (calculated from histos) and its fit for a given slot and
channel.
Starting values for threshold and width may be suggested.
The fit results are stored in global arrays Athresholds and Awidths
TH1D *plot_eta1(int adc_slot, int adc_chan, int tdc_cut=1350, int tdc_slot = -1, int tdc_chan = -1){
  const int n=200;
  const int bin_amplitude = 5;
  last_tdc_cut = tdc_cut;
  TString draw, cut, title;
  float ratio;

  TCanvas *ceta_raw = new TCanvas("ceta_raw","ceta_raw");

  const int nbin = n/bin_amplitude;

  int Ncut, Nnocut;

  if (tdc_chan == -1) tdc_chan = calc_tdc_chan(adc_slot,adc_chan);
  if (tdc_slot == -1) tdc_slot = calc_tdc_slot(adc_slot,adc_chan);

  TH1D eta_raw_adctot = TH1D("eta_raw_adctot","eta_raw_adctot",nbin, 0, n);
  TH1D eta_raw_adccut = TH1D("eta_raw_adccut","eta_raw_adccut",nbin, 0, n);
  eta_raw_adctot.Sumw2();
  eta_raw_adccut.Sumw2();

  MyStyle->SetStatX(0.98);
  MyStyle->SetStatY(0.3);

  draw.Form("adc[%d][%d]>>eta_raw_adctot",adc_slot,adc_chan);
  t->Draw(draw,"");
  draw.Form("adc[%d][%d]>>eta_raw_adccut", adc_slot,adc_chan);
  cut.Form("tdct[%d][%d]>%d",tdc_slot, tdc_chan, tdc_cut);
  t->Draw(draw,cut);

  TEfficiency *pEff = 0;
  if (TEfficiency::CheckConsistency(eta_raw_adccut, eta_raw_adctot,"w") ){
    pEff = new TEfficiency(eta_raw_adccut, eta_raw_adctot);
  }
  float Ath_estimate = 70.;
  float Aw_estimate = 20.;
  //cout<<"Estimates: Ath "<< Ath_estimate <<" Aw "<< Aw_estimate <<endl;
  TF1 *eta = new TF1("eta",_eta_func,40.,140, 2);
  eta->SetParameter(0, Ath_estimate);
  eta->SetParameter(1, Aw_estimate);
  eta->SetParLimits(0, 20, Ath_estimate+50);
  //    eta->SetParLimits(1, 20, 1000);
  pEff->Draw("AP");
  pEff->Fit(eta,"R"); // preliminary fit
  float Ath = eta->GetParameter(0);
  float Aw = eta->GetParameter(1);
  float Ath_err = eta->GetParError(0);
  float Aw_err = eta->GetParError(1);
  float chi2 = eta->GetChisquare();
  delete eta;
  float Athmin = Ath-50<0?0:Ath-50;
  TF1 *eta = new TF1("eta",_eta_func,Athmin, Athmin + 100, 2);
  eta->SetLineColor(kRed);
  eta->SetParNames("Ath","Aw");
  eta->SetParameter(0, Ath);
  eta->SetParameter(1, Aw);
  eta->SetParLimits(0, Athmin, Athmin+100);
  //    eta->SetParLimits(1, 15, 50);
  pEff->Draw("AP");
  pEff->Fit(eta,"R"); // better fit
  ceta_raw->Modified();
  // titles
  TH1D *htemp=(TH1D*)gPad->GetPrimitive("eta_raw_adctot_clone");
  title.Form("Run %d slot %d channel %d TDC > %d", run,adc_slot,adc_chan,tdc_cut);
  htemp->SetTitle(title);
  // update values
  Ath = eta->GetParameter(0);
  Aw = eta->GetParameter(1);
  Ath_err = eta->GetParError(0);
  Aw_err = eta->GetParError(1);
  chi2 = eta->GetChisquare();
  // fill global arrays
  Athresholds[adc_slot][adc_chan] = Ath;
  Awidths[adc_slot][adc_chan] = Aw;
  cout << "slot "<< adc_slot <<" chan "<< adc_chan <<" Saved Ath = " << Ath << " Aw = " << Aw <<" in arrays Athresholds and Awidths"<<endl;

  return htemp;
}

//====================================================================
//                    Work in progress section
//                                or,
//               the crashy, the buggy and the ugly (c)
//====================================================================


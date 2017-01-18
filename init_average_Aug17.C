#include <stdio.h>

static const int NUMADCSLOTS = 4;
static const int NUMTDCSLOTS = 3;
static const int NUMCHANA = 64;
static const int NUMCHANT = 96;
static const float adc_charge = 50*1e-15; // 50 fC, corresponding to an ADC chan
static const float e = 1.6e-19; // C, electron charge
static const int NPMT = 14;
static const int NPIXEL = 16;
static const int NBARS = NPMT*NPIXEL;
static const int NNEIGH = 4; // number of neighbour bars to consider in cuts.
// to set unconnected bar, see init()

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

Int_t connected[NBARS];

int run; // run number

TTree *t;
Int_t raw[NUMADCSLOTS][NUMCHANA], ped[NUMADCSLOTS][NUMCHANA], trail[NUMTDCSLOTS][NUMCHANT], lead[NUMTDCSLOTS][NUMCHANT];

void init(int runno){
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
}

// MAPPING
// PMT POS BARS    ADC    TDC   NINO
//              SLOT CH SLOT CH  ID
//==================================
//  1  1    0..15  0  0    0  0  18
//  2  8   16..31  1 48    1 16   3
//  3  2   32..47  0 16    0 16   2
//  4  9   48..63  2  0    1 32  17
//  5  3   64..79  0 32    0 32   4
//  6 10   80..95  2 16    1 48  19
//  7  4  96..111  0 48    0 48   5
//  8 11 112..127  2 32    1 64   8
//  9  5 128..143  1  0    0 64   7
// 10 12 144..159  2 48    1 80  11
// 11  6 160..175  3 32    0 80  13
// 12 13 176..191  3  0    2  0  12
// 13  7 192..207  3 48    1  0   6
// 14 14 208..223  3 16    2 16  14



Int_t aslot(Int_t bar){
  Int_t slot;
  if ( (bar >= 0 && bar <= 15) || 
       (bar >= 32 && bar <= 47) || 
       (bar >= 64 && bar <= 79) ||
       (bar >= 96 && bar <= 111) ) slot = 0;
  if ( (bar >= 16 && bar <= 31) ||
       (bar >= 128 && bar <= 143) ) slot = 1;
  if ( (bar >= 48 && bar <= 63) ||
       (bar >= 80 && bar <= 95) ||
       (bar >= 112 && bar <= 127) ||
       (bar >= 144 && bar <= 159) ) slot = 2;
  if ( (bar >= 160 && bar <= 223) ) slot = 3;
  return slot;
}

Int_t achan(Int_t bar){
  Int_t chan = -1;
  if (bar >= 0 && bar <= 15) chan = bar;
  if (bar >= 16 && bar <= 31) chan = bar + 32;
  if (bar >= 32 && bar <= 47) chan = bar - 16;
  if (bar >= 48 && bar <= 63) chan = bar - 48;
  if (bar >= 64 && bar <= 79) chan = bar - 32;
  if (bar >= 80 && bar <= 95) chan = bar - 64;
  if (bar >= 96 && bar <= 111) chan = bar - 48;
  if (bar >= 112 && bar <= 127) chan = bar - 80;
  if (bar >= 128 && bar <= 143) chan = bar - 128;
  if (bar >= 144 && bar <= 159) chan = bar - 96;
  if (bar >= 160 && bar <= 175) chan = bar - 128;
  if (bar >= 176 && bar <= 191) chan = bar - 176;
  if (bar >= 192 && bar <= 207) chan = bar - 144;
  if (bar >= 208 && bar <= 223) chan = bar - 192;
  return chan;
}

TCanvas * Avert(Int_t thresh_current_bar, Int_t thres_neigh_bars){
  Int_t slot, chan, bar, obar;
  TH1D htmp ("htmp","htmp",2200,-200,2000);
  TString buffer, cut, drawme;

  Double_t average, x[NBARS], y[NBARS];

  TCanvas * cAvert = new TCanvas("cAvert", "cAvert");

  memset(x, 0, NBARS * sizeof(Int_t));
  memset(y, 0, NBARS * sizeof(Int_t));

  cout << "# of neighbours considered: "<<  NNEIGH << endl;

  for (bar = 0; bar < NBARS; bar++){
    if ( connected[bar] ){
      slot = aslot(bar);
      chan = achan(bar);
      drawme.Form("adc[%d][%d]>>htmp", slot, chan);
      cut.Form("adc[%d][%d] > %d", slot, chan, thresh_current_bar); // cut on current

      // let's build the cut on neighbours
      for( obar = bar - NNEIGH; obar < bar + NNEIGH + 1; obar++){
	buffer = "";
	if (obar >= 0 && obar < NBARS && obar != bar){ // ok, valid bar
	  if ( connected[obar] ){
	    slot = aslot(obar);
	    chan = achan(obar);
	    // for vertical tracks:
	    buffer.Form(" && adc[%d][%d] < %d", slot, chan, thres_neigh_bars);
	    // for horizontal tracks:
	    //buffer.Form("adc[%d][%d] > %d", slot, chan, thres_neigh_bars);
	  }
	}
	cut = cut + buffer; // append cut of current neighbour
      }

      // apply cuts, store resulting histogram in htmp
      t->Draw(drawme, cut);

      // here we can just extract mean, or fit htmp, or whatever
      average = htmp.GetMean();
    } else { // current bar not connected
      average = 0;
    }

    cout << "bar "<< bar;
    // cout <<" is supposed to be";
    // if ( connected[bar] == 0) {
    //   cout<<" NOT connected" <<endl;
    // } else {
    //   cout << " connected to slot "<< aslot(bar) <<" chan "<< achan(bar) << endl;
    //   cout <<" cut applied: "<< cut << endl;
    // }

    // t->Draw(drawme);
    // cout <<" average before = "<< htmp.GetMean();

    cout <<" average = "<< average << endl;

    // save stuff in arrays
    x[bar] = bar;
    y[bar] = average;
  }

  // time to plot
  TGraph * gAvert = new TGraph( NBARS, x, y);
  gAvert->SetTitle("<A-vert> vs pixel number; pixel; <ADC>");
  gAvert->Draw("AP");
  return cAvert;
}

//====================================================================
//                    Work in progress section
//                                or,
//               the crashy, the buggy and the ugly (c)
//====================================================================

// TCanvas * display_single(Int_t event, Int_t threshold = 100){
//   TCanvas *cdisplays = new TCanvas("cdisplays","cdisplays");

//   int debug = 1;
//   if (debug) cout << "Entering display_single" <<endl;

//   Double_t bar[NBARS], yped[NBARS], yraw[NBARS], delta[NBARS];

//   memset( bar, 0, NBARS );
//   memset( yped, 0, NBARS);
//   memset( yraw, 0, NBARS);
//   memset( delta, 0, NBARS);

//   if (debug) cout << "memset ok" <<endl;

//   t->GetEntry(event);

//   if (debug) cout << "getentry ok" <<endl;

//   Int_t max = 0;
  
//   for(int pmt=1; pmt<NPMT+1; pmt++){
//     for(int pixel=1; pixel<NPIXEL+1; pixel++){
 
//       int global = (pmt-1)*NPIXEL + (pixel-1); // dummy mapping

//       int slota = handmapping_adc_slot( pmt );
//       int chana = handmapping_adc_chan( pmt, pixel);
//       int slott = handmapping_tdc_slot( pmt );
//       int chant = handmapping_tdc_chan( pmt, pixel);

//       if (global < NBARS){
// 	bar[global] = global +1;   // dummy mapping
// 	yped[global] = ped[slota][chana];
// 	yraw[global] = raw[slota][chana];
// 	delta[global] = trail[slott][chant]-lead[slott][chant];

// 	if (yped[global] > max) max = yped[global];
// 	if (yraw[global] > max) max = yraw[global];

// 	if (yped[global] - yped[global] > threshold)
// 	  cout << "PMT "<< pmt << " pixel "<< pixel << " bar "<< global <<" ADC "<< yped[global] << " ped "<< yped[global] << " ADC-ped "<< yped[global] - yped[global] << endl;
//       }
//     }
//   }

//   if (debug) cout << "loop ok" <<endl;

//   TGraph *gped = new TGraph(NBARS, bar, yped);
//   TGraph *graw = new TGraph(NBARS, bar, yraw);
//   TGraph *gdelta = new TGraph(NBARS, bar, delta);
//   Float_t size = 0.3;
//   Int_t style = 21;

//   TString buffer;
//   buffer.Form("ADC Run %d, event %d; Bar", run, event);

//   if (debug) cout << "tgraph created ok" <<endl;

//   gped->SetTitle(buffer);
//   gped->SetMarkerColor(kRed);
//   gped->SetMarkerStyle(style);
//   gped->SetMarkerSize(size);

//   graw->SetTitle(buffer);
//   graw->SetMarkerColor(kGreen);
//   graw->SetMarkerStyle(style);
//   graw->SetMarkerSize(size);

//   if (debug) cout << "tgraph props ok" <<endl;

//   TMultiGraph *mg1 = new TMultiGraph();
//   mg1->Add(gped,"p");
//   mg1->Add(graw,"p");
//   mg1->SetTitle(buffer);

//   if (debug) cout << "mg1 ok" <<endl;

//   buffer.Form("TDC Run %d, event %d; Bar", run, event);

//   gdelta->SetTitle(buffer);
//   gdelta->SetMarkerColor(kBlack);
//   gdelta->SetMarkerStyle(style);
//   gdelta->SetMarkerSize(size);

//   if (debug) cout << "delta ok" <<endl;

//   TMultiGraph *mg2 = new TMultiGraph();
//   mg2->Add(gdelta);
//   mg2->SetTitle(buffer);

//   if (debug) cout << "mg2 ok" <<endl;

//   TLegend* leg1 = new TLegend(.9, .45, .99, .65);
//   leg1->AddEntry(gped,"Pedestal","p");
//   leg1->AddEntry(graw,"ADC","p");
//   TLegend* leg2 = new TLegend(.9, .45, .99, .65);
//   leg2->AddEntry(gdelta,"TDC length","p");

//   if (debug) cout << "legend ok" <<endl;

//   //  TPad * pad = new TPad();
//   cdisplays->Divide(1,2);
//   //pad->Divide(1,2);

//   if (debug) cout << "pad divide ok" <<endl;

//   cdisplays->cd(1);
//   //pad->cd(1);

//   mg1->Draw("ap");
//   cdisplays->Modified(); cdisplays->Update();
//   //pad->Modified(); pad->Update();
//   mg1->GetXaxis()->SetLimits(0,NPMT*NPIXEL);
//   mg1->SetMaximum(1400);
//   //gPad->Modified();

//   //gped->Draw("AP");
//   //graw->Draw("AP same");

//   leg1->Draw();

//   if (debug) cout << "sub pad 1 ok" <<endl;

//   cdisplays->cd(2); 
//   //pad->cd(2);

//   mg2->Draw("ap");
//   cdisplays->Modified(); cdisplays->Update();
//   //pad->Modified(); pad->Update();
//   mg2->GetXaxis()->SetLimits(0,NPMT*NPIXEL);
//   mg2->SetMaximum(1400);
//   //gPad->Modified();
//   // gdelta->Draw("ap");

//   leg2->Draw();

//   if (debug) cout << "sub pad 2 ok" <<endl;

//   //  return pad;
//   return cdisplays;
// }


// void display(Int_t event, Int_t thres = 0){
//   int endme=0;
//   int nentries = t->GetEntries();
//   char choice;
//   //  TPad * pad;
//   TCanvas *pad;

//   pad = display_single(event, thres);

//   while (endme == 0){
//     pad->Draw();
//     cout << "(n)ext event, (p)revious event, (g)oto event, (q)uit"<<endl;

//     //choice = getchar();
//     system("stty raw");//seting the terminal in raw mode
//     while(1) {
//       choice=getchar();
//       if(choice=='q'){ 
// 	system("stty cooked");
	
//       }
//       //  printf("you pressed %c\n ",choice);
//     }


//     switch(choice){
//     case 'p': 
//       if (event > 1) {
// 	pad = display_single(event -1, thres);
//       } else {
// 	cout << "Already at first event" << endl;
//       };
//       break;
//     case 'n':
//       if (event < nentries) {
// 	pad = display_single(event +1, thres);
//       } else {
// 	cout << "Already at last event"<<endl;
//       };
//       break;
//     case 'q':
//       endme = 1;
//       break;
//     case 'g':
//       cin >> event;
//       if ( event >= 1 && event <= nentries ) {
// 	pad = display_single(event, thres);
//       } else {
// 	cout << "Event number not valid" << endl;
//       };
//       break;
//     default:
//       cout << "Invalid choice"<< endl;
//   }
// }
// }


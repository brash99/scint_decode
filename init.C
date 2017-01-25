#include "TROOT.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TList.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include "mapping.C"

using namespace std;

static const Int_t NUMADCSLOTS = 4;
static const Int_t NUMTDCSLOTS = 3;
static const Int_t NUMCHANA = 64;
static const Int_t NUMCHANT = 96; // in the root file, we have 2*3 times this. read Fastbus_main1.C
static const float adc_charge = 50*1e-15; // 50 fC, corresponding to an ADC chan
static const float e = 1.6e-19; // C, electron charge
static const Int_t xcanvas = 800; // width of canvases
static const Int_t ycanvas = 800; // height of canvases

Int_t run; // run number. Used in the titles

// pedestal
Int_t ped[NUMADCSLOTS][NUMCHANA];

// generic use
float buffer[10000];
TStyle *MyStyle = new TStyle("MyStyle","MyStyle");
Int_t currentpad;

// trees
TTree *t; // main one, read from file, raw adc & tdc data
TTree *tfriend; // friend of t, mapped, gain-corrected data...well, not yet

Int_t adc[NUMADCSLOTS][NUMCHANA];
Int_t adcraw[NUMADCSLOTS][NUMCHANA];
Int_t tdct[NUMTDCSLOTS][NUMCHANT];
Int_t tdcl[NUMTDCSLOTS][NUMCHANT];
Int_t tdcn[NUMTDCSLOTS][NUMCHANT];

// NOTE: mapping is done in Fastbus_main1.C
// The routine is copied here for reference and checks
// MAPPING
// We count pmt and pixel from 1, as opposed to 0, I prefer to just
// add +1 to the constants here as opposed to every loop, and in the
// definition of the branchs
static const Int_t MAXPMTIDS = 1 + 400; // number of user ID's
static const Int_t MAXPMTS = 1 + NUMADCSLOTS * NUMCHANA / 16; // max number of PMT's which can be read by DAQ
static const Int_t NPIXELS = 1 + NPIXEL;

Int_t map_nino_id[MAXPMTS];
Int_t map_pmt_id[MAXPMTIDS]; // index = PMT "user ID"
// pmt, pixel -> slot, chan
Int_t map_adc_slot[MAXPMTS][NPIXELS];
Int_t map_adc_chan[MAXPMTS][NPIXELS];
Int_t map_tdc_slot[MAXPMTS][NPIXELS];
Int_t map_tdc_chan[MAXPMTS][NPIXELS];
Int_t map_adc_pmt[NUMADCSLOTS][NUMCHANA];
Int_t map_adc_pixel[NUMADCSLOTS][NUMCHANA];
Int_t map_tdc_pmt[NUMTDCSLOTS][NUMCHANT];
Int_t map_tdc_pixel[NUMTDCSLOTS][NUMCHANT];

const char *mapfile = "sbs_mapping.txt"; // config file
Int_t map_len; // number of read pmt config lines

// just to make cint compiler happy
Int_t mapping_read(const char*);

// reads root file corresponding to run runno, extracts the tree
void init(Int_t runno){
  TString filename;
  filename.Form("sbs_%d_14.root",runno);
  TFile *_file0 = TFile::Open(filename);

  run=runno;
	
  //MyStyle->SetStatFont(42);
  MyStyle->SetTitleFontSize(0.08);
  MyStyle->SetTitleX(0.15);
  MyStyle->SetTitleY(0.99);
  MyStyle->SetStatW(0.9);
  MyStyle->SetMarkerStyle(6);
  gStyle->SetCanvasDefH(xcanvas);
  gStyle->SetCanvasDefW(ycanvas);
  gROOT->SetStyle("MyStyle");

  // read tree
  // TTree *t = (TTree *)_file0->Get("t"); 
  t = (TTree *)_file0->Get("t"); 
  t->SetBranchAddress("adc",&adc);
  t->SetBranchAddress("adcraw",&adcraw);
  t->SetBranchAddress("tdct",&tdct);
  t->SetBranchAddress("tdcl",&tdcl);
  t->SetBranchAddress("tdcn",&tdcn);

  Int_t nentries = t->GetEntries();
  cout << "Found " << nentries << " events"<<endl;

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
// automatically; if values different than -1 is passed to
// tdc_chan_start/slot, those values will be used instead. May be
// useful during testing of ADC/TDC channels
TCanvas *plot_adc(Int_t pmt=1, Int_t tdc_min=1300, Int_t tdc_width=200){
	
	Int_t adc_slot = handmapping_adc_slot(pmt);
	Int_t adc_chan_start = handmapping_adc_chan(pmt,1);
	Int_t tdc_slot = handmapping_tdc_slot(pmt);
	Int_t tdc_chan_start = handmapping_tdc_chan(pmt,1);
	Int_t pixel1 = handmapping_pmt_pixel1(pmt);
	Int_t pixel2 = handmapping_pmt_pixel2(pmt);
	
	TString cut, draw, draw1, title;
	title.Form("run_%d_ADC",run);
	TCanvas *cADC= new TCanvas("cADC",title,xcanvas,ycanvas);
	
	TH1D *htmpa[16];//=new TH1D("htmp","htmp",nbin,min,max);
	TH1D *htmpb[16];//=new TH1D("htmp1","htmp1",nbin,min,max);
	
	TString tmpentry;
	MyStyle->SetStatX(0.9);
	MyStyle->SetStatY(0.6);
	MyStyle->SetStatW(0.4);

        Int_t nbin=600;
        Int_t min=-100, max=1000;
        for(Int_t icounter=1;icounter<17;icounter++){
                tmpentry.Form("htmpa%d", icounter);
                htmpa[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpb%d", icounter);
                htmpb[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                htmpa[icounter - 1]->SetLineColor(kBlue);
                htmpb[icounter - 1]->SetLineColor(kRed);
	  	title.Form("Run %d ADC slot %d chan %d: %d < tdc < %d",run,adc_slot, icounter, tdc_min,tdc_min+tdc_width);
                htmpa[icounter - 1]->SetTitle(title);
                htmpb[icounter - 1]->SetTitle(title);
        }

        Int_t nentries=t->GetEntries();

        for (Int_t id=10;id<nentries;id++){
        t->GetEntry(id);
	
	 // fill histos
	 for (Int_t i=adc_chan_start; i<adc_chan_start+16; i++){ 

	  currentpad = i - adc_chan_start + 1;

	  Int_t itdc = tdc_chan_start + currentpad-1;
	  Int_t iadc = i;

	  htmpa[currentpad-1]->Fill(adc[adc_slot][iadc]);
          if(tdcl[tdc_slot][itdc]>tdc_min&&tdcl[tdc_slot][itdc]<tdc_min+tdc_width){
                htmpb[currentpad-1]->Fill(adc[adc_slot][iadc]);
          }

	 }
	}

	cADC->Clear();
	cADC->Divide(4,4) ;

	//plot histos
	Int_t icount=0;
	for (Int_t i=0; i<16; i++){

	  if(i != pixel1-1 && i != pixel2-1) {
	  	cout<<"into loop 2, i = " << i << endl;

	  	cADC->cd( icount + 1 );
	  	gPad->SetLogy();

	  	//cADC->Update();

	  	Int_t entries = htmpa[i]->GetEntries();
	  	float mean = htmpa[i]->GetMean(1);
	  	float RMS = htmpa[i]->GetRMS(1);

	  	cout << entries <<" "<< mean <<" "<< RMS <<endl;

	  	htmpa[i]->SetStats(0);
	  	// current->Modified();

	  	htmpa[i]->Draw();
	 	htmpb[i]->Draw("same");
	  
		icount++;
	  }


	};
	
	title.Form("run_%d_ADC_slot_%d_chan_%d_%d_tdc_min_%d_max_%d.png",
		   run,adc_slot,adc_chan_start,adc_chan_start+15,tdc_min,tdc_min+tdc_width);
	cADC->Print(title);
	cADC->cd(0);
	return cADC;
}

//--------------------------------------------------------------------
//                        ADC Ratio Plots
//--------------------------------------------------------------------
// useful during testing of ADC/TDC channels
TCanvas *plot_ratio(Int_t pmt=1, Int_t tdc_min=1300, Int_t tdc_width=200){
	
	Int_t adc_slot = handmapping_adc_slot(pmt);
	Int_t adc_chan_start = handmapping_adc_chan(pmt,1);
	Int_t tdc_slot = handmapping_tdc_slot(pmt);
	Int_t tdc_chan_start = handmapping_tdc_chan(pmt,1);
	Int_t pixel1 = handmapping_pmt_pixel1(pmt);
	Int_t pixel2 = handmapping_pmt_pixel2(pmt);
	
	TString cut, draw, draw1, title;
	title.Form("run_%d_ADC",run);
	TCanvas *cRATIO= new TCanvas("cRATIO",title,xcanvas,ycanvas);
	
	TH1D *htmpa[16];//=new TH1D("htmp","htmp",nbin,min,max);
	TH1D *htmpb[16];//=new TH1D("htmp1","htmp1",nbin,min,max);
	TH1D *htmpc[16];//=new TH1D("htmp1","htmp1",nbin,min,max);
	
	TString tmpentry;
	MyStyle->SetStatX(0.9);
	MyStyle->SetStatY(0.6);
	MyStyle->SetStatW(0.4);

        Int_t nbin=25;
        Int_t min=0, max=100;
        for(Int_t icounter=1;icounter<17;icounter++){
                tmpentry.Form("htmpa%d", icounter);
                htmpa[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpb%d", icounter);
                htmpb[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpc%d", icounter);
                htmpc[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                htmpa[icounter - 1]->SetLineColor(kBlue);
                htmpb[icounter - 1]->SetLineColor(kRed);
                htmpc[icounter - 1]->SetLineColor(kRed);
                title.Form("Run %d ADCRATIO slot %d chan %d: %d < tdc < %d",run,adc_slot, icounter, tdc_min,tdc_min+tdc_width);
                htmpa[icounter - 1]->SetTitle(title);
                htmpb[icounter - 1]->SetTitle(title); 
                htmpc[icounter - 1]->SetTitle(title); 
        }

        Int_t nentries=t->GetEntries();

        for (Int_t id=10;id<nentries;id++){
        t->GetEntry(id);
	
	 // fill histo 
	 for (Int_t i=adc_chan_start; i<adc_chan_start+16; i++){ 

	  currentpad = i - adc_chan_start + 1;

	  Int_t itdc = tdc_chan_start + currentpad-1;
	  Int_t iadc = i;

          htmpa[currentpad-1]->Fill(adc[adc_slot][itdc]);
          if(tdcl[tdc_slot][itdc]>tdc_min&&tdcl[tdc_slot][itdc]<tdc_min+tdc_width){
                htmpb[currentpad-1]->Fill(adc[adc_slot][iadc]);
          }

	 }
	}

	cRATIO->Clear();
	cRATIO->Divide(4,4) ;

	//plot histos
	Int_t icount=0;
	TF1 *myfit = new TF1("myfit","1.0-0.5*ROOT::Math::erfc((x-[0])/[1])",0,1);
	myfit->SetParName(0,"Mean");
	myfit->SetParName(1,"Width");

	for (Int_t i=0; i<16; i++){

	  if(i != pixel1-1 && i != pixel2-1) {
	  	cout<<"into loop 2, i = " << i << endl;

	  	cRATIO->cd( icount + 1 );
	  	//gPad->SetLogy();

	  	//cRATIO->Update();

	  	Int_t entries = htmpa[i]->GetEntries();
	  	float mean = htmpa[i]->GetMean(1);
	  	float RMS = htmpa[i]->GetRMS(1);

	  	cout << entries <<" "<< mean <<" "<< RMS <<endl;

	  	htmpb[i]->SetStats(0);
	  	// current->Modified();

	 	//htmpb[i]->Draw();
		myfit->SetParameter(0,40.0);
		myfit->SetParameter(1,10.0);
	  	htmpc[i] = (TH1D*)htmpb[i]->Clone();
	  	htmpc[i]->Divide(htmpa[i]);
		htmpc[i]->Fit("myfit");
	  
		icount++;
	  }


	};
	
	title.Form("run_%d_RATIO_slot_%d_chan_%d_%d_tdc_min_%d_max_%d.png",
		   run,adc_slot,adc_chan_start,adc_chan_start+15,tdc_min,tdc_min+tdc_width);
	cRATIO->Print(title);
	cRATIO->cd(0);
	return cRATIO;
}

//--------------------------------------------------------------------
//                        ADC Occupancy Plot
//--------------------------------------------------------------------
// useful during testing of ADC/TDC channels
TCanvas *plot_occupancy(Int_t adc_cut=40, Int_t multiplicity_cut = 12, Int_t tdc_min=1300, Int_t tdc_width=200){

  	Int_t nentries = t->GetEntries();
	TString cut, title;
	title.Form("run_%d_Occupancy",run);
	Int_t nbin=196;
	Int_t min=1, max=197;
	Int_t nbinm=11;
	Int_t minm=-1, maxm=10;
	TH1D *hoccupancy = new TH1D("hoccupancy","hoccupancy",nbin,min,max);
	TH1D *hmultiplicity = new TH1D("hmultiplicity","hmultiplicity",nbinm,minm,maxm);
	TCanvas *cOCCUPANCY= new TCanvas("cOCCUPANCY",title,xcanvas,ycanvas);
        
	for (Int_t id=1;id<nentries;id++){
	//for (Int_t id=10;id<51;id++){
		t->GetEntry(id);

		Int_t nmultiplicity=0;
		Int_t good_paddle[100];
		for(Int_t icount=0;icount<100;icount++){good_paddle[icount]=-1;}

		for (Int_t pmt=1; pmt<15; pmt++){
	
			Int_t adc_slot = handmapping_adc_slot(pmt);
			Int_t adc_chan_start = handmapping_adc_chan(pmt,1);
			Int_t tdc_slot = handmapping_tdc_slot(pmt);
			Int_t tdc_chan_start = handmapping_tdc_chan(pmt,1);
			Int_t pixel1 = handmapping_pmt_pixel1(pmt);
			Int_t pixel2 = handmapping_pmt_pixel2(pmt);
	
			TString tmpentry;
			MyStyle->SetStatX(0.9);
			MyStyle->SetStatY(0.9);
			MyStyle->SetStatW(0.3);
			// fill histos
	
        		Int_t currentpixel=0;
			for (Int_t i=adc_chan_start; i<adc_chan_start+16; i++){ 

	  			currentpad = i - adc_chan_start + 1;
	  			if(currentpad != pixel1 && currentpad != pixel2) {						
					currentpixel++;
	  
	  				Int_t paddle = (pmt-1)*14+(15-currentpixel);

	  				Int_t itdc = tdc_chan_start + currentpad-1;
	  				Int_t iadc = i;

	  				//if((adc[adc_slot][iadc]>50||(tdcl[tdc_slot][itdc]>1300&&tdcl[tdc_slot][itdc]<1500))) cout << " Event = " << id << "PMT = " << pmt << " Pixel = " << currentpad << " ADC = " << adc[adc_slot][iadc] << " TDC = " << tdcl[tdc_slot][itdc] <<endl;
	  				if (tdct[tdc_slot][itdc] > tdc_min && tdcl[tdc_slot][itdc] < tdc_min+tdc_width) {
						nmultiplicity++;
						if (adc[adc_slot][iadc] > adc_cut) {
							good_paddle[nmultiplicity-1]=paddle;
						}
	  				}
	  			}		
			}
		}
		if(nmultiplicity>0&&nmultiplicity<=multiplicity_cut){
			for(Int_t icount=0;icount<nmultiplicity;icount++){hoccupancy->Fill(good_paddle[icount]);}
		}
		hmultiplicity->Fill(nmultiplicity);
	}

	cOCCUPANCY->Clear();
	cOCCUPANCY->Divide(1,2) ;

	//plot histos
	
	title.Form("run_%d_OCCUPANCY_tdc_min_%d_max_%d.png",
		   run,tdc_min,tdc_min+tdc_width);
	cOCCUPANCY->Print(title);
	cOCCUPANCY->cd(1);
	hoccupancy->Draw();
	hoccupancy->GetXaxis()->SetNdivisions(14,14,0,0);
	hoccupancy->SetLineColor(kBlue);
	gPad->SetGridx();
	cOCCUPANCY->cd(2);
	hmultiplicity->Draw();
	hmultiplicity->SetLineColor(kBlue);
	return cOCCUPANCY;
}
//                        TDC hits
//--------------------------------------------------------------------
//
TCanvas *plot_tdc_hits(Int_t pmt=1, Int_t adc_cut=40){
	
        Int_t adc_slot = handmapping_adc_slot(pmt);
	Int_t adc_chan_start = handmapping_adc_chan(pmt,1);
	Int_t tdc_slot = handmapping_tdc_slot(pmt);
	Int_t tdc_chan_start = handmapping_tdc_chan(pmt,1);
	Int_t pixel1 = handmapping_pmt_pixel1(pmt);
	Int_t pixel2 = handmapping_pmt_pixel2(pmt);
	
	TString cut, draw, draw1, title;
	title.Form("run_%d_TDCHITS",run);
	TCanvas *cTDCHITS= new TCanvas("cTDCHITS",title,xcanvas,ycanvas);
	
	TH1D *htmpa[16];//=new TH1D("htmpa","htmpa",nbin,min,max);
	TH1D *htmpb[16];//=new TH1D("htmpb","htmpb",nbin,min,max);
	
	TString tmpentry;
	MyStyle->SetStatW(0.4);
	MyStyle->SetStatX(0.6);
	MyStyle->SetStatY(0.9);
	
        Int_t nbin=11;
	Int_t min=-1, max=10;
        for(Int_t icounter=1;icounter<17;icounter++){
		tmpentry.Form("htmpa%d", icounter);
		htmpa[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
		tmpentry.Form("htmpb%d", icounter);
		htmpb[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
	  	htmpa[icounter - 1]->SetLineColor(kBlue);
	  	htmpb[icounter - 1]->SetLineColor(kRed);
	  	title.Form("Run %d TDC HITS slot %d chan %d adc > %d",run,tdc_slot, icounter, adc_cut);
	  	htmpa[icounter - 1]->SetTitle(title);
	  	htmpb[icounter - 1]->SetTitle(title);
	}

	Int_t nentries=t->GetEntries();

	for (Int_t id=10;id<nentries;id++){
	t->GetEntry(id);

	 // fill histos
	 for (Int_t i=tdc_chan_start; i<tdc_chan_start+16; i++){

	  currentpad = i - tdc_chan_start + 1;

	  Int_t itdc = i;
	  Int_t iadc = adc_chan_start + currentpad-1;

	  htmpa[currentpad-1]->Fill(tdcn[tdc_slot][itdc]);
	  if(adc[adc_slot][iadc]>adc_cut){
		htmpb[currentpad-1]->Fill(tdcn[tdc_slot][itdc]);
	  }
	  
	 }
	}

	cTDCHITS->Clear();
	cTDCHITS->Divide(4,4) ;

	//plot histos
	Int_t icount=0;
	for (Int_t i=0; i<16; i++){

	  if(i != pixel1-1 && i != pixel2-1) {
	  cout<<"into loop 2, i = " << i << endl;

	  cTDCHITS->cd( icount + 1 );
	  gPad->SetLogy();

	  //cADC->Update();

	  Int_t entries = htmpa[i]->GetEntries();
	  float mean = htmpa[i]->GetMean(1);
	  float RMS = htmpa[i]->GetRMS(1);

	  cout << entries <<" "<< mean <<" "<< RMS <<endl;

	  htmpa[i]->SetStats(0);
	  // current->Modified();

	  htmpa[i]->Draw();
	  htmpb[i]->Draw("same");

	  icount++;
	  }
	};
	
	title.Form("run_%d_TDCHITS_slot_%d_chan_%d_%d_adc_cut_%d.png",
		   run,tdc_slot,tdc_chan_start,tdc_chan_start+15,adc_cut);
	cTDCHITS->Print(title);
	cTDCHITS->cd(0);
	return cTDCHITS;

}


//                        TDC only
//--------------------------------------------------------------------
// For interactive use:
// plotcanvas divided in 4x4 TDC plots of 16 channels starting from
// adc_chan_start, of given slot. TDC is plotted with and without the
// following cut: adc value > adc_cut
// By default, tdc starting chan = adc starting chan; if a value
// different than -1 is passed to adc_chan_start, that value will be used
// instead. May be useful during testing of ADC/TDC channels
// Analogously for adc_slot
TCanvas *plot_tdc(Int_t pmt=1, Int_t adc_cut=40, Int_t tdc_min=1300, Int_t tdc_width=200){
	
        Int_t adc_slot = handmapping_adc_slot(pmt);
	Int_t adc_chan_start = handmapping_adc_chan(pmt,1);
	Int_t tdc_slot = handmapping_tdc_slot(pmt);
	Int_t tdc_chan_start = handmapping_tdc_chan(pmt,1);
	Int_t pixel1 = handmapping_pmt_pixel1(pmt);
	Int_t pixel2 = handmapping_pmt_pixel2(pmt);
	
	TString cut, draw, draw1, title;
	title.Form("run_%d_TDC",run);
	TCanvas *cTDC= new TCanvas("cTDC",title,xcanvas,ycanvas);
	
	TH1D *htmpa[16];//=new TH1D("htmpa","htmpa",nbin,min,max);
	TH1D *htmpb[16];//=new TH1D("htmpb","htmpb",nbin,min,max);
	//TH1D *htmp = new TH1D("htmp","htmp",nbin,min,max);
	
	TString tmpentry;
	MyStyle->SetStatW(0.4);
	MyStyle->SetStatX(0.6);
	MyStyle->SetStatY(0.6);

        Int_t nbin=tdc_width;
        Int_t min=tdc_min, max=tdc_min+tdc_width;
        for(Int_t icounter=1;icounter<17;icounter++){
                tmpentry.Form("htmpa%d", icounter);
                htmpa[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpb%d", icounter);
                htmpb[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                htmpa[icounter - 1]->SetLineColor(kBlue);
                htmpb[icounter - 1]->SetLineColor(kRed);
                title.Form("Run %d TDC slot %d chan %d adc > %d",run,tdc_slot, icounter, adc_cut);
                htmpa[icounter - 1]->SetTitle(title);
                htmpb[icounter - 1]->SetTitle(title);
        }

        Int_t nentries=t->GetEntries();

        for (Int_t id=10;id<nentries;id++){
        t->GetEntry(id);

	 // fill histos
	 for (Int_t i=tdc_chan_start; i<tdc_chan_start+16; i++){

	  currentpad = i - tdc_chan_start + 1;

	  Int_t itdc = i;
	  Int_t iadc = adc_chan_start + currentpad-1;

	  htmpa[currentpad-1]->Fill(tdcl[tdc_slot][itdc]);
	  if(adc[adc_slot][iadc]>adc_cut){
		htmpb[currentpad-1]->Fill(tdcl[tdc_slot][itdc]);
	  }
	  
	 }
	}

	cTDC->Clear();
	cTDC->Divide(4,4) ;

	//plot histos
	Int_t icount=0;
	for (Int_t i=0; i<16; i++){

	  if(i != pixel1-1 && i != pixel2-1) {
	  cout<<"into loop 2, i = " << i << endl;

	  cTDC->cd( icount + 1 );
	  //gPad->SetLogy();

	  //cADC->Update();

	  Int_t entries = htmpa[i]->GetEntries();
	  float mean = htmpa[i]->GetMean(1);
	  float RMS = htmpa[i]->GetRMS(1);

	  cout << entries <<" "<< mean <<" "<< RMS <<endl;

	  htmpa[i]->SetStats(0);
	  // current->Modified();

	  htmpa[i]->Draw();
	  htmpb[i]->Draw("same");

	  icount++;
	  }
	};
	
	title.Form("run_%d_TDC_slot_%d_chan_%d_%d_adc_cut_%d.png",
		   run,tdc_slot,tdc_chan_start,tdc_chan_start+15,adc_cut);
	cTDC->Print(title);
	cTDC->cd(0);
	return cTDC;

}


//====================================================================
//                        Ancillary routines
//     (normally there is no need to call them interactively)
//====================================================================

Int_t main(Int_t argc, char* argv[]){
  if (argc != 1) return 0;
  Int_t run = atoi(argv[1]);
  init(run);
  return 0;
}


//====================================================================
//                    Work in progress section
//                                or,
//               the crashy, the buggy and the ugly (c)
//====================================================================

Int_t mapping_read(const char* mapfilein){

  // pre-processing, let's remove all non-data lines; the result is
  // saved to a temp file
  char file[100];
  sprintf(file,"%s.tmp", mapfilein);

  gROOT->ProcessLine(Form(".! awk '/^[ \t]*[0123456789]/' %s > %s",mapfilein,file));
  
  // now we read the temp file
  Int_t pmtID, ADCslot, ADCchan, TDCslot, TDCchan, NINOID;

  ifstream in(file);

  // resetting maps
  for (Int_t ii = 0; ii < MAXPMTS; ii++){
      map_pmt_id[ ii ] = -1;
      map_nino_id[ ii ] = -1;
      for (Int_t jj = 0; jj < NPIXELS; jj++){
	map_adc_slot[ ii ][ jj ] = -1;
	map_tdc_slot[ ii ][ jj ] = -1;
	map_adc_chan[ ii ][ jj ] = -1;
	map_tdc_chan[ ii ][ jj ] = -1;
      }
  }
  for (Int_t ii = 0; ii < NUMADCSLOTS; ii++){
    for (Int_t jj = 0; jj < NUMCHANA; jj++){
      map_adc_pmt[ ii ][ jj ] = -1;
      map_adc_pixel[ ii ][ jj ] = -1;
    }
  }
  for (Int_t ii = 0; ii < NUMTDCSLOTS; ii++){
    for (Int_t jj = 0; jj < NUMCHANT; jj++){
      map_tdc_pmt[ ii ][ jj ] = -1;
      map_tdc_pixel[ ii ][ jj ] = -1;
    }
  }

  Int_t nlines = 0;
  while ( !in.eof() ){
    nlines++; // we start at 1

    in >> pmtID;
    in >> ADCslot;
    in >> ADCchan;
    in >> TDCslot;
    in >> TDCchan;
    in >> NINOID;

    //cout<< "read: " << pmtID <<" "<< ADCslot <<" "<<  ADCchan <<" "<<  TDCslot <<" "<<  TDCchan <<" "<<  NINOID <<endl;

    map_nino_id[ nlines ] = NINOID;
    map_pmt_id[ pmtID ] = nlines;
 
    // pmt, pixel -> slot, chan
    //-------------------------
    // unroll channel mappings
    for (Int_t i = 1; i < NPIXELS; i++){
      // slots do not change with pixel
      map_adc_slot[ nlines ][ i ] = ADCslot;
      map_tdc_slot[ nlines ][ i ] = TDCslot;
      // channels do
      map_adc_chan[ nlines ][ i ] = ADCchan + i;
      map_tdc_chan[ nlines ][ i ] = TDCchan + i;
    }

    // slot, chan -> pmt, pixel
    //-------------------------
    for (Int_t pixel = 1; pixel < NPIXELS; pixel++){
      map_adc_pmt  [ ADCslot ][ ADCchan + pixel - 1] = pmtID;
      map_adc_pixel[ ADCslot ][ ADCchan + pixel - 1] = pixel;  

      map_tdc_pmt  [ TDCslot ][ TDCchan + pixel - 1] = pmtID;
      //for (Int_t j=0; j<6; j++) map_tdc_pixel[ TDCslot ][ TDCchan + pixel + j ]= pixel;
      map_tdc_pixel[ TDCslot ][ TDCchan + pixel - 1]= pixel;
    }
  }
  in.close();
  if (nlines > MAXPMTS) cout << "WARNING: number of pmt config lines are > MAXPMTS = " << MAXPMTS << endl;

  if (0){
    cout << "Test pmt, pixel -> slot, chan" <<endl;
    for (Int_t i=0; i<nlines+1; i++){
      for (Int_t j=0; j<16; j++){
	cout << "PMT "<< i <<" pixel "<< j 
	     <<" adc slot = " << map_adc_slot[ i ][ j ]
	     <<" tdc slot = " << map_tdc_slot[ i ][ j ]
	     <<" adc chan = " << map_adc_chan[ i ][ j ]
	     <<" tdc chan = " << map_tdc_chan[ i ][ j ]
	     << endl;
      }
    }
    cout << "Test slot, chan -> pmt, pixel" <<endl;
    cout << "Now ADC:"<<endl;
    for (Int_t i=0; i<NUMADCSLOTS; i++){
      for (Int_t j=0; j<NUMCHANA; j++){
	cout << "adc slot "<< i <<" chan "<< j 
	     <<" pmt = " << map_adc_pmt[ i ][ j ]
	     <<" pixel = " << map_adc_pixel[ i ][ j ]
	     << endl;
      }
    }
    cout << "Now TDC:"<<endl;
    for (Int_t i=0; i<NUMTDCSLOTS; i++){
      for (Int_t j=0; j<NUMCHANT; j++){
	cout << "tdc slot "<< i <<" chan "<< j 
	     <<" pmt = " << map_tdc_pmt[ i ][ j ]
	     <<" pixel = " << map_tdc_pixel[ i ][ j ]
	     << endl;
      }
    }
  }
  
  map_len = nlines;
  return nlines;
}

TCanvas *check_mapping(Int_t adc_thres){
  Int_t global_chan = 0;
  Int_t counter_sc = 0; // sc = slot & chan
  Int_t counter_pp = 0; // pp = pmt & pixel
  Int_t counter_bm = 0; // bm = before mapping
  Int_t counter_NC_adc = 0; // NC = not connected adc channel
  Int_t counter_NC_tdc = 0; // NC = not connected ydc channel
  float mean;
  // float rms;
  TString plotme, cut;
  Int_t ndata;
  Int_t tdcslot, tdcchan, pmt, pixel;

  Int_t badvalue = -10;
  Int_t nentries = t->GetEntries();

  static const Int_t numchans = NUMADCSLOTS * NUMCHANA;
  Int_t adc_chans_sc[ numchans ];
  Int_t tdc_chans_sc[ numchans ];
  Int_t adc_chans_pp[ numchans ];
  Int_t tdc_chans_pp[ numchans ];
  Int_t adc_chans_bm[ numchans ];
  Int_t tdc_chans_bm[ numchans ];

  memset( adc_chans_sc, badvalue, numchans * sizeof(Int_t) );
  memset( tdc_chans_sc, badvalue, numchans * sizeof(Int_t) );
  memset( adc_chans_pp, badvalue, numchans * sizeof(Int_t) );
  memset( tdc_chans_pp, badvalue, numchans * sizeof(Int_t) );
  memset( adc_chans_bm, badvalue, numchans * sizeof(Int_t) );
  memset( tdc_chans_bm, badvalue, numchans * sizeof(Int_t) );

  TH1D htmp("htmp","htmp",2200, -200, 2000);

  mapping_read(mapfile);

  cout << "Building slot & channel -> pmt & pixel, ALSO Before Mapping versions" << endl;
  // SLOT & CHAN version
  // cout << "...ADC"<<endl;
  // // ADC
  // for (Int_t slot = 0; slot < NUMADCSLOTS; slot++){
  //   for (Int_t chan = 0; chan < NUMCHANA; chan++){
  //     global_chan = slot * NUMCHANA + chan; 
  //     pmt = map_adc_pmt[ slot ][ chan ];
  //     pixel = map_adc_pixel[ slot ][ chan ];

  //     // check whether the current channel is plugged to a pixel
  //     if (pmt < 0 || pixel < 0){
  // 	cout << "WARNING: ADC slot "<< slot <<" channel "<< chan <<" would map to pmt "<< pmt <<" pixel "<< pixel <<", which means it is not connected" << endl;
  // 	cout << "...we don't mind here, because we are building the Before Mapping plot" << endl;
  // 	counter_NC_adc++;
  // 	//	continue;
  //     }

  //     // Before mapping version
  //     plotme.Form("adc[%d][%d]>>htmp", slot, chan);
  //     cut.Form("adc[%d][%d] > %d", slot, chan, adc_thres);

  //     ndata = t->Draw(plotme, cut, "goff");
  //     mean = htmp.GetMean();
  //     if (ndata > nentries) cout << "NM slot = "<< slot <<" chan = "<< chan << " mean adc = " << mean << " ndata = " << ndata << endl;
  //     //      rms = htmp->GetRMS();
      
  //     if (ndata > 0){
  // 	counter_nm++;
  // 	adc_chans_nm[ global_chan ] = mean;
  //     }

  //     // slot & chan -> pmt & pixel
  //     plotme.Form("a[%d][%d]>>htmp", pmt, pixel);
  //     cut.Form("a[%d][%d] > %d", pmt, pixel, adc_thres);

  //     ndata = t->Draw(plotme, cut, "goff");
  //     mean = htmp.GetMean();

  //     // check whether the current channel is plugged to a pixel
  //     if (pmt < 0 || pixel < 0){
  // 	// A warning has been already issued during the Before Mapping
  // 	// version, so we can just skip the channel at this point
  // 	continue;
  //     }

  //     if (ndata > nentries) cout << "SC slot = "<< slot <<" chan = "<< chan << " mean adc = " << mean << " ndata = " << ndata << endl;
  //     //      rms = htmp->GetRMS();
      
  //     if (ndata > 0 && mean > adc_thres){
  // 	counter_sc++;
  // 	adc_chans_sc[ global_chan ] = mean;

  // 	// TDC only done if ADC is good
  // 	pmt = map_tdc_pmt[ slot ][ chan ];
  // 	pixel = map_tdc_pixel[ slot ][ chan ];

  // 	// check whether the current channel is plugged to a pixel
  // 	if (pmt < 0 || pixel < 0){
  // 	  cout << "ERROR: TDC slot "<< slot <<" channel "<< chan <<" maps to pmt "<< pmt <<" pixel "<< pixel <<", which means it is not connected" << endl;
  // 	  continue;
  // 	}

  // 	plotme.Form("t[%d][%d]>>htmp", pmt, pixel);
  // 	cut.Form("a[%d][%d] > %d", pmt, pixel, adc_thres);

  // 	ndata = t->Draw(plotme, cut, "goff");
  // 	mean = htmp.GetMean();
  // 	if (ndata > nentries) cout << "Warning: SC slot = "<< slot <<" chan = "<< chan << " pmt " << pmt << " pixel "<< pixel << " mean tdc = " << mean << " ndata = " << ndata << endl;
  // 	//      rms = htmp->GetRMS();

  // 	if (ndata > 0){
  // 	  tdc_chans_sc[ global_chan ] = mean;
  // 	}
  //     }
  //   }
  // }

  // TDC
  // cout << "...TDC"<<endl;
  
  // for (Int_t slot = 0; slot < NUMTDCSLOTS; slot++){
  //   for (Int_t chan = 0; chan < NUMCHANT; chan++){
  //     global_chan = slot * NUMCHANT + chan; 
  //     pmt = map_tdc_pmt[ slot ][ chan ];
  //     pixel = map_tdc_pixel[ slot ][ chan ];

  //     // check whether the current channel is plugged to a pixel
  //     if (pmt < 0 || pixel < 0){
  // 	cout << "WARNING: TDC slot "<< slot <<" channel "<< chan <<" would map to pmt "<< pmt <<" pixel "<< pixel <<", which means it is not connected" << endl;
  // 	cout << "...we don't mind here, because we are building the Before Mapping plot" << endl;
  // 	counter_NC_tdc++;
  // 	//	continue;
  //     }

  //     // Before mapping version
  //     plotme.Form("tdc[%d][%d]>>htmp", slot, chan);
  //     cut.Form("adc[%d][%d] > %d", slot, chan, adc_thres);

  //     ndata = t->Draw(plotme, cut, "goff");
  //     mean = htmp.GetMean();
  //     if (ndata > nentries) cout << "NM slot = "<< slot <<" chan = "<< chan << " mean adc = " << mean << " ndata = " << ndata << endl;
  //     //      rms = htmp->GetRMS();
      
  //     if (ndata > 0){
  // 	counter_nm++;
  // 	adc_chans_nm[ global_chan ] = mean;
  //     }

  //     // slot & chan -> pmt & pixel
  //     plotme.Form("a[%d][%d]>>htmp", pmt, pixel);
  //     cut.Form("a[%d][%d] > %d", pmt, pixel, adc_thres);

  //     ndata = t->Draw(plotme, cut, "goff");
  //     mean = htmp.GetMean();

  //     // check whether the current channel is plugged to a pixel
  //     if (pmt < 0 || pixel < 0){
  // 	// A warning has been already issued during the Before Mapping
  // 	// version, so we can just skip the channel at this point
  // 	continue;
  //     }

  //     if (ndata > nentries) cout << "SC slot = "<< slot <<" chan = "<< chan << " mean adc = " << mean << " ndata = " << ndata << endl;
  //     //      rms = htmp->GetRMS();
      
  //     if (ndata > 0 && mean > adc_thres){
  // 	counter_sc++;
  // 	adc_chans_sc[ global_chan ] = mean;

  // 	// TDC only done if ADC is good
  // 	pmt = map_tdc_pmt[ slot ][ chan ];
  // 	pixel = map_tdc_pixel[ slot ][ chan ];

  // 	// check whether the current channel is plugged to a pixel
  // 	if (pmt < 0 || pixel < 0){
  // 	  cout << "ERROR: TDC slot "<< slot <<" channel "<< chan <<" maps to pmt "<< pmt <<" pixel "<< pixel <<", which means it is not connected" << endl;
  // 	  continue;
  // 	}

  // 	plotme.Form("t[%d][%d]>>htmp", pmt, pixel);
  // 	cut.Form("a[%d][%d] > %d", pmt, pixel, adc_thres);

  // 	ndata = t->Draw(plotme, cut, "goff");
  // 	mean = htmp.GetMean();
  // 	if (ndata > nentries) cout << "Warning: SC slot = "<< slot <<" chan = "<< chan << " pmt " << pmt << " pixel "<< pixel << " mean tdc = " << mean << " ndata = " << ndata << endl;
  // 	//      rms = htmp->GetRMS();

  // 	if (ndata > 0){
  // 	  tdc_chans_sc[ global_chan ] = mean;
  // 	}
  //     }
  //   }
  // }


  //     // slot & chan -> pmt & pixel
  //     pmt = map_adc_pmt[ slot ][ chan ];
  //     pixel = map_adc_pixel[ slot ][ chan ];
  //     plotme.Form("a[%d][%d]>>htmp", pmt, pixel);
  //     cut.Form("a[%d][%d] > %d", pmt, pixel, adc_thres);

  //     ndata = t->Draw(plotme, cut, "goff");
  //     mean = htmp.GetMean();

  //     if (ndata > nentries) cout << "SC slot = "<< slot <<" chan = "<< chan << " mean adc = " << mean << " ndata = " << ndata << endl;
  //     //      rms = htmp->GetRMS();
      
  //     if (ndata > 0 && mean > adc_thres){
  // 	counter_sc++;
  // 	adc_chans_sc[ global_chan ] = mean;

  // 	// TDC only done if ADC is good
  // 	pmt = map_tdc_pmt[ slot ][ chan ];
  // 	pixel = map_tdc_pixel[ slot ][ chan ];

  // 	// check whether the current channel is plugged to a pixel
  // 	if (pmt < 0 || pixel < 0){
  // 	  cout << "ERROR: TDC slot "<< slot <<" channel "<< chan <<" maps to pmt "<< pmt <<" pixel "<< pixel <<", which means it is not connected" << endl;
  // 	  continue;
  // 	}

  // 	plotme.Form("t[%d][%d]>>htmp", pmt, pixel);
  // 	cut.Form("a[%d][%d] > %d", pmt, pixel, adc_thres);

  // 	ndata = t->Draw(plotme, cut, "goff");
  // 	mean = htmp.GetMean();
  // 	if (ndata > nentries) cout << "SC slot = "<< slot <<" chan = "<< chan << " mean tdc = " << mean << " ndata = " << ndata << endl;
  // 	//      rms = htmp->GetRMS();

  // 	if (ndata > 0){
  // 	  tdc_chans_sc[ global_chan ] = mean;
  // 	}
  //     }
  //   }
  // }

  cout << "Building pmt & pixel -> slot & chan version" << endl;
  // PMT & PIXEL version
  for (Int_t pmt = 1; pmt < map_len; pmt++){
    for (Int_t pixel = 1; pixel < NPIXELS; pixel++){
      global_chan = pmt * (NPIXELS-1) + (pixel-1);  // -1 because...read its def

      // ADC
      plotme.Form("a[%d][%d]>>htmp", pmt, pixel);
      cut.Form("a[%d][%d] > %d", pmt, pixel, adc_thres);

      ndata = t->Draw(plotme, cut, "goff");
      if (ndata > 0) counter_pp++;
      mean = htmp.GetMean();
      //      rms = htmp->GetRMS();

      if (ndata > nentries) cout << "PP pmt = "<< pmt <<" pixel = "<< pixel << " mean adc = " << mean << " ndata = " << ndata << endl;

      if (ndata > 0 && mean > adc_thres){
	counter_pp++;
	adc_chans_pp[ global_chan ] = mean;
      }

      // TDC
      plotme.Form("t[%d][%d]>>htmp", pmt, pixel);
      cut.Form("a[%d][%d] > %d", pmt, pixel, adc_thres);

      ndata = t->Draw(plotme, cut, "goff");
      mean = htmp.GetMean();
      //      rms = htmp->GetRMS();

      if (ndata > nentries) cout << "PP pmt = "<< pmt <<" pixel = "<< pixel << " mean tdc = " << mean << " ndata = " << ndata << endl;

      if (ndata > 0){
	tdc_chans_pp[ global_chan ] = mean;
      }

    }
  }

  // TGraph *g_nm = new TGraph( numchans, tdc_chans_nm, adc_chans_nm);
  // g_nm->SetTitle("ADC vs TDC, NOT MAPPED slot&chan version");
 
  // TGraph *g_sc = new TGraph( numchans, tdc_chans_sc, adc_chans_sc);
  // g_sc->SetTitle("ADC vs TDC mapping check, slot&chan -> pmt&pixel");

  TGraph *g_pp = new TGraph( numchans, tdc_chans_pp, adc_chans_pp);
  g_pp->SetTitle("ADC vs TDC mapping check, pmt&pixel -> slot&chan");


  TCanvas *ccheck_map = new TCanvas("ccheck_map", "Mapping checks", xcanvas, ycanvas);
  ccheck_map->Divide(2,2);
  // ccheck_map->cd(1); g_nm->Draw("ap");
  // ccheck_map->cd(3); g_sc->Draw("ap");
  ccheck_map->cd(4); g_pp->Draw("ap");

  cout << "Counters:" << endl;
  // cout << " Not mapped = " << counter_nm << endl;
  // cout << " slot&chan -> pmt&pixel = " << counter_sc << endl;
  cout << " pmt&pixel -> slot&chan = " << counter_pp << endl;

  return ccheck_map;
}

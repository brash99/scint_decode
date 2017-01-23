// A simple decoder and histogrammer for fastbus data.
// R. Michaels, Oct 24, 2014.
// developed later by Dasuni, et al.

#define MAXROC     32
#define MAXHITS    100
#define MAXEVENTS  1000000
// constants used to distinguish which case the bit 16 of TDC data word is referring to:
#define TRAILING_EDGE 0
#define LEADING_EDGE  1

#define NUMCHANA   64 // ADC chans per slot
#define NUMCHANT   96 //*2*3 // (TDC chans per slot) * (leading+trailing edge) * (number of pulses recorded)

// If the following values change, (un)comment as well the relevant initial values of
//     nbad_adc[NUMADCSLOTS], nbad_tdc[NUMTDCSLOTS], slotindTDC[NUMTDCSLOTS]
// which are below
#define NUMADCSLOTS 4   // ADC # of slots. =4 for test crate, =4 ? for the other ones
#define NUMTDCSLOTS 3   // TDC # of slots. =3 for test crate, =9 for the other ones

int quick;

#include <iostream>
#include <string>
#include <vector>
#include "THaCodaFile.h"
#include "THaEtClient.h"
#include "TString.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include <vector>
#include "Riostream.h"
#include "TTree.h"

using namespace std;

void usage();
//void decode(int* data, TTree *tree_adc, TTree *tree_tdc);
void decode(int* data, TTree *tree);
void analysis();
void clear();
void pedsup();
void readpedsup();

// MAPPING
// Int_t handmapping_adc_slot(Int_t);
// Int_t handmapping_adc_chan(Int_t, Int_t);
// Int_t handmapping_tdc_slot(Int_t);
// Int_t handmapping_tdc_chan(Int_t, Int_t);
#include "mapping.C"

// for tests
Int_t match_tdc_slot_for_adc_slot_chan(Int_t, Int_t);
Int_t match_tdc_chan_for_adc_slot_chan(Int_t, Int_t);
//
// Global data 
// -----------
Int_t evlen, evtype, evnum;
Int_t *irn, *rocpos, *roclen;
Int_t myroc, myslot[4], mychan[4];
Int_t ievent = 0;
Int_t fHasHeader=1;   // =1 for 1877 and 1881, but =0 for 1875
Int_t fWordSeen;
Int_t debug=0;
Int_t numLocal = 0;
Int_t branchnum = 999;

//ADC
Int_t numslothitsADC[NUMADCSLOTS];
Int_t adcdat[NUMADCSLOTS][NUMCHANA];
Int_t adcchan[NUMADCSLOTS][NUMCHANA];
Int_t nbad_adc[NUMADCSLOTS]={0,0,0,0}; // for test crate. Should be fine for the other ones also
Int_t slotindADC[NUMADCSLOTS];
Int_t nbadevents_adc = 0;

// mapping
Int_t a[NPMT+1][NPIXEL+1];
Int_t t[NPMT+1][NPIXEL+1];
Int_t l[NPMT+1][NPIXEL+1];


#define MAXPED 2000 // pedestal is searched up to this position
Int_t PEDSUP; // 0 = read pedestal file, 1 = calculate and write file. Read from argv
Int_t pedestal_threshold[NUMADCSLOTS][NUMCHANA];

TH1F  *hadc1;
TH1F  *hadc2[4];
TH2F  *hadc3[NUMADCSLOTS];
TH2F  *h2adc;
TH1F  *hadcped[NUMADCSLOTS][NUMCHANA];
TH1F  *htdc[NUMTDCSLOTS][NUMCHANT];
TH1F  *hbranch;

//TDC
Int_t numslothitsTDC[NUMTDCSLOTS];
Int_t tdctdat[NUMTDCSLOTS][NUMCHANT];
Int_t tdcldat[NUMTDCSLOTS][NUMCHANT];
Int_t tdcchan[NUMTDCSLOTS][NUMCHANT];
// for test crate:
Int_t nbad_tdc[NUMTDCSLOTS]={0,0,0};
Int_t slotindTDC[NUMTDCSLOTS]={14,15,16};
// for the other ones:
//Int_t nbad_tdc[NUMTDCSLOTS]={0,0,0,0,0,0,0,0,0};
//Int_t slotindTDC[NUMTDCSLOTS]={11,12,13,14,15,16,17,18,19} // for the other ones
Int_t nbadevents_tdc = 0;

TH1F  *htdc1;
TH1F  *hleadingtime;
TH1F  *htdcpulsewidth;

// Trees variables
Int_t tree_adc_slot, tree_adc_chan;
Int_t tree_adc_hits[NUMADCSLOTS][NUMCHANA];
Int_t tree_adc_hits_minus_pedestal[NUMADCSLOTS][NUMCHANA];
Int_t tree_adc_pedestal_chan[NUMADCSLOTS][NUMCHANA];
Int_t tree_tdc_slot, tree_tdc_chan;
Int_t tree_tdc_hits[NUMTDCSLOTS][NUMCHANT];
Int_t tree_tdc_nhits[NUMTDCSLOTS][NUMCHANT];
Int_t tree_tdc_le_hits[NUMTDCSLOTS][NUMCHANT];
Int_t isleadingedge;
bool is_good; // true if event has both one leading and one trailing edge
// -----------

int main(int argc, char* argv[])
{
  //============ Initializations ============
  THaCodaData *coda;      
  char ctitle[100];
  char dtitle[100];
  char rtitle[100];
  char htitle[100];
  char cpedtitle[100];
  char hpedtitle[100];
  Int_t maxevent = 10000;
  Int_t runno = 1000;
  Int_t istatus;
  Int_t lslot;

  int choice1 = 1; // =1 CODA file, else ET connection
  irn = new Int_t[MAXROC];
  rocpos = new Int_t[MAXROC];
  roclen = new Int_t[MAXROC];

  //default args
  PEDSUP = 0;
  maxevent = MAXEVENTS;

  // ----- argc, argv -----
  if (argc >= 3) {
    runno = atoi(argv[1]); // required arg: run#
    myroc = atoi(argv[2]); // required arg: roc#
  }
  if (argc > 3) // optional arg: maxevent
    maxevent = atoi(argv[3]);
  if (argc > 4)  // optional arg: pedsup
    PEDSUP = atoi(argv[4]);
  if (argc > 5 || argc < 3) { // wrong # of args
    usage();
    return 1;
  }  

  cout << "args elaborated" << endl;

  if (myroc == 14){
    lslot = 20;
  }
  else if (myroc == 5){
    lslot = 17;
  }
  else if (myroc == 6){
    lslot = 11;
  }
  else if (myroc == 7){
    lslot = 17;
  }
  else if (myroc == 15){
    lslot = 5;
  }
  else if (myroc == 16){
    lslot = 14;
  }
  else if (myroc == 17){
    lslot = 8;
  }
  else{
    lslot = 11;
  }

  myslot[0]= lslot; myslot[1]= lslot+1; myslot[2]= lslot+2; myslot[3]= lslot+3;
  mychan[0] = 0;  mychan[1] = 0; mychan[2] = 0; mychan[3] = 0; // defaults

  cout << "my roc is " << myroc << endl;

  for (int i = 0; i < NUMADCSLOTS; i++){
    slotindADC[i] = lslot;
    lslot++;
  }

  cout << "Fastbus analysis "<<endl;
  cout << "Events to process "<<maxevent<<endl;

  // ----- init root -----
  sprintf(rtitle,"sbs_%d_%i.root",runno,myroc);
  sprintf(dtitle,"../data/scint_%d.dat",runno);
  TROOT fbana("fbroot","Hall A SCINT analysis for Fastbus");
  TFile hfile(rtitle,"RECREATE","SBS data");

  // Trees
  //TTree *tree_adc;
  //TTree *tree_tdc;

  char tempstr[16];
  // tree_adc = new TTree("adc","Raw ADC data tree");
  // tree_adc->Branch("event"	    ,&ievent,        "event/I");
  // tree_adc->Branch("slot"	    ,&tree_adc_slot, "slot/I");
  // tree_adc->Branch("chan" 	    ,&tree_adc_chan, "chan/I");
  // sprintf(tempstr,"value[%d][%d]/I" ,NUMADCSLOTS, NUMCHANA);
  // tree_adc->Branch("vr"	            ,tree_adc_hits, tempstr); // raw values
  // tree_adc->Branch("v"	            ,tree_adc_hits_minus_pedestal, tempstr);
  // sprintf(tempstr,"pedestal[%d][%d]/I", NUMADCSLOTS, NUMCHANA);
  // tree_adc->Branch("p"              ,tree_adc_pedestal_chan,tempstr);
 
  // tree_tdc = new TTree("tdc", "Raw TDC data tree");
  // tree_tdc->Branch("event"	   ,&ievent,        "event/I");
  // tree_tdc->Branch("slot"	   ,&tree_tdc_slot, "slot/I");
  // tree_tdc->Branch("chan" 	   ,&tree_tdc_chan, "chan/I");
  // sprintf(tempstr,"value[%d][%d]/I",NUMTDCSLOTS, NUMCHANT);
  // tree_tdc->Branch("v"		   ,tree_tdc_hits, tempstr);
  // tree_tdc->Branch("lev"	   ,tree_tdc_le_hits, tempstr);
  // tree_tdc->Branch("is_leadingedge",&isleadingedge, "is_leadingedge/I");

  TTree *tree = new TTree("t","ADC & TDC data tree");
  tree->Branch("event"	    ,&ievent,        "event/I");
  sprintf(tempstr,"value[%d][%d]/I" ,NUMADCSLOTS, NUMCHANA);
  tree->Branch("adcraw"	            ,tree_adc_hits, tempstr); // adc raw values
  tree->Branch("adc"	            ,tree_adc_hits_minus_pedestal, tempstr); // adc values - pedestal
  sprintf(tempstr,"pedestal[%d][%d]/I", NUMADCSLOTS, NUMCHANA);
  tree->Branch("ped"              ,tree_adc_pedestal_chan,tempstr); //pedestal
  sprintf(tempstr,"value[%d][%d]/I",NUMTDCSLOTS, NUMCHANT);
  tree->Branch("tdcn"		   ,tree_tdc_nhits, tempstr); // tdc number of leading and trailing edges
  tree->Branch("tdct"		   ,tdctdat, tempstr); // tdc trailing edge
  tree->Branch("tdcl"         	   ,tdcldat, tempstr); //tdc leading edge
  // tree->Branch("is_good"           ,&is_good,"is_good/B");

  // mapping
  sprintf(tempstr,"value[%d][%d]/I", NPMT+1, NPIXEL+1);
  tree->Branch("a", a, tempstr);
  tree->Branch("t", t, tempstr);
  tree->Branch("l", l, tempstr);

  // ----- init output -----
  // ADC
  sprintf(ctitle,"ADC hits per slot ");
  hadc1 = new TH1F("hadc1",ctitle,31,-1,30.);
  sprintf(ctitle,"ADC data on slot %d channel %d",myslot[0],mychan[0]);
  hadc2[0] = new TH1F("hadc2_0",ctitle,200,0,4000.);
  sprintf(ctitle,"ADC data on slot %d channel %d",myslot[1],mychan[1]);
  hadc2[1] = new TH1F("hadc2_1",ctitle,200,0,4000.);
  sprintf(ctitle,"ADC data on slot %d channel %d",myslot[2],mychan[2]);
  hadc2[2] = new TH1F("hadc2_2",ctitle,200,0,4000.);
  sprintf(ctitle,"ADC data on slot %d channel %d",myslot[3],mychan[3]);
  hadc2[3] = new TH1F("hadc2_3",ctitle,200,0,4000.);

  sprintf(ctitle,"ADC data on 2 channnels");
  h2adc = new TH2F("h2adc",ctitle,200,0,4000.,200,0,4000.);
  sprintf(ctitle,"Branch #");
  hbranch = new TH1F("hbranch",ctitle,10,0,10.);

  for (int i = 0; i < NUMADCSLOTS;i++) {
    sprintf(ctitle,"ADC data on all channels of slot %d",slotindADC[i]);
    sprintf(htitle,"ADC_chan_slot_%d",slotindADC[i]);
    hadc3[i] = new TH2F(htitle,ctitle,200,0,2000.,64,0,64.);
  }
  for (int i = 0; i < NUMADCSLOTS;i++) {
    for (int j = 0; j < NUMCHANA;j++) {
      sprintf(cpedtitle,"myADC data on chan %d of slot %d", j, slotindADC[i]);
      sprintf(hpedtitle,"myADC_chan_%d_slot_%d",j,slotindADC[i]);
      hadcped[i][j] = new TH1F(hpedtitle,cpedtitle,MAXPED,0,MAXPED);
    }
  }
  for (int i = 0; i < NUMTDCSLOTS;i++) {
    for (int j = 0; j < NUMCHANT;j++) {
      sprintf(cpedtitle,"TDC data on chan %d of slot %d", j, slotindTDC[i]);
      sprintf(hpedtitle,"TDC_chan_%d_slot_%d",j,slotindTDC[i]);
      htdc[i][j] = new TH1F(hpedtitle,cpedtitle,1000,0,1000.);
    }
  }

  // TDC
  sprintf(ctitle,"TDC hits per slot ");
  htdc1 = new TH1F("htdc1",ctitle,31,-1,30.);
  htdc1->GetXaxis()->SetTitle("slot");
  htdc1->GetYaxis()->SetTitle("counts");
  sprintf(ctitle,"TDC leading edge times ");
  hleadingtime = new TH1F("hleadingtime",ctitle,1000,0,1000.);
  hleadingtime->GetXaxis()->SetTitle("channel (0.5 ns each)");
  hleadingtime->GetYaxis()->SetTitle("counts");
  sprintf(ctitle,"TDC pulse width ");
  htdcpulsewidth = new TH1F("htdcpulsewidth",ctitle,1000,0,1000.);
  htdcpulsewidth->GetXaxis()->SetTitle("channel (0.5 ns each)");
  htdcpulsewidth->GetYaxis()->SetTitle("counts");

  // ----- init input -----
  if (choice1 == 1) {  // CODA File
    // CODA file "run.dat" may be a link to CODA datafile on disk
    TString filename(dtitle);

    coda = new THaCodaFile();
    if (coda->codaOpen(filename) != 0) {
      cout << "ERROR: Cannot open CODA data" << endl;
      goto end1;
    }
  } else {         // Online ET connection
    int mymode = 1;
    TString mycomputer("sbs1");
    TString mysession("sbsfb1");

    coda = new THaEtClient();
    if (coda->codaOpen(mycomputer, mysession, mymode) != 0) {
      cout << "ERROR:  Cannot open ET connection" << endl;
      goto end1;
    }
  }

  // ----- Pedestal file -----
  if (PEDSUP == 0) {
    cout << "Reading pedestal file because PEDSUP = 0 " <<endl;
    readpedsup();
  };

  // ============ End of initializations ============

  cout <<"Ready to go! starting event loop" << endl;
  // Loop over events

  for (int iev = 0; iev < maxevent; iev++)  {//the loop over the events
    if (iev > 0 && ((iev%1000)==0) ) printf("%d events\n",iev);

    clear();
    istatus = coda->codaRead();  

    if (istatus != 0) {  // EOF or no data
      if ( istatus == -1) {
        if (choice1 == 1) {
          cout << "End of CODA file. Bye bye." << endl;
        }
        if (choice1 == 2) cout << "CODA/ET not running. Bye bye." << endl;
      } else {
        cout << "ERROR: codaRread istatus = " << hex << istatus << endl;
      }
      goto end1;
    } else {   // we have data ...
      ievent++;
      //      decode( coda->getEvBuffer(), tree_adc, tree_tdc);
      decode( coda->getEvBuffer(), tree);
      if (evtype < 10) 
        analysis();
    }
  }

 end1:

  cout << endl;
  if (PEDSUP==1) {
    cout << "CALCULATING pedestal because PEDSUP = 1" << endl;
    pedsup();
  }

  printf("\n");
  printf("Number of events analyzed = %d \n",ievent);
  for (int i=0;i<NUMADCSLOTS;i++){
    printf("In Slot %d, Number of events with bad adc = %d \n",slotindADC[i],nbad_adc[i]);
  }

  // TDC
  for (int i=0;i<NUMTDCSLOTS;i++){
    printf("In Slot %d, Number of events with bad tdc = %d \n",slotindTDC[i],nbad_tdc[i]);
  }
  cout << endl;

  cout<<"Events with no adc = "<<  nbadevents_adc <<endl;
  cout<<"Events with no tdc = "<<  nbadevents_tdc <<endl;


  if (debug){
  for (int slot=0; slot<NUMADCSLOTS; slot++)
    for (int chan=0; chan<NUMCHANA; chan++)
      if (pedestal_threshold[slot][chan] < 1)
	cout << "WARNING: pedestal of ADC slot "<< slot <<" channel "<< chan <<" is "<< pedestal_threshold[slot][chan]<<endl;
  }

  coda->codaClose();

  hfile.Write();
  hfile.Close();


  // for (int pmt = 1; pmt < 15; pmt++)
  //   for (int pixel = 1; pixel<17; pixel++){
  //     int as = handmapping_adc_slot(pmt);
  //     int ac = handmapping_adc_chan(pmt, pixel);
  //     int ts = handmapping_tdc_slot(pmt);
  //     int tc = handmapping_tdc_chan(pmt, pixel);
  //     printf("pmt %2d pix %2d aslot %2d achan %2d tslot %2d tchan %2d\n", pmt, pixel, as, ac, ts, tc);
  //   }



  return 0;
}; //end of main function

void usage() {  
  cout << "Usage:  'fbana [runno] [roc#] [maxevents] [PEDSUP] ' " << endl;
  cout << "Need [runno] [roc#] arguments "<<endl;
  cout << "[maxev] optional argument, default "<< MAXEVENTS << endl;
  cout << "[PEDSUPl] optional argument. 0 means use pedestal saved in file ped_test.dat."<<endl;
  cout << "                             1 means calculate it and save it." <<endl;
  cout << "                             Default = 0" <<endl;
  cout << "Can have at most four arguments "<<endl;
}; 

void clear() {
  fWordSeen = 0;
  //ADC
  memset (adcdat, 0, NUMADCSLOTS*NUMCHANA*sizeof(Int_t));
  memset (adcchan, 0, NUMADCSLOTS*NUMCHANA*sizeof(Int_t));
  memset (numslothitsADC, 0, NUMADCSLOTS*sizeof(Int_t));

  //TDC
  memset (numslothitsTDC, 0, NUMTDCSLOTS*sizeof(Int_t));
  memset (tdctdat, 0, NUMTDCSLOTS*NUMCHANT*sizeof(Int_t));
  memset (tdcldat, 0, NUMTDCSLOTS*NUMCHANT*sizeof(Int_t));
  memset (tdcchan, 0, NUMTDCSLOTS*NUMCHANT*sizeof(Int_t));

  // mapping
  memset (a, 0, (NPMT+1)*(NPIXEL+1)*sizeof(Int_t));
  memset (t, 0, (NPMT+1)*(NPIXEL+1)*sizeof(Int_t));
  memset (l, 0, (NPMT+1)*(NPIXEL+1)*sizeof(Int_t));
}

//void decode (int* data, TTree *tree_adc, TTree *tree_tdc) {
void decode (int* data, TTree *tree) {
  // ----- init -----
  Int_t ichan = 0, rdata = 0;
  evlen = data[0] + 1;
  evtype = data[1]>>16;
  evnum = data[4];
  static int dodump = 0;  // dump the raw data

  Int_t edgetype = 0;
  int slotnew=0;

  if (evtype > 10) return;
  //if (evnum < 10 ) return; 
// Sometimes, the first few events are bad. No idea why. Let's skip them

  // init all chans at -100 (non physical value)
  for (int j=0;j<NUMADCSLOTS;j++){
    for (int ii=0;ii<NUMCHANA;ii++){
      tree_adc_hits[j][ii] = -100;
    }
  }
  for (int j=0;j<NUMTDCSLOTS;j++){
    for (int ii=0;ii<NUMCHANT;ii++){
      tree_tdc_nhits[j][ii] = 0;
      tree_tdc_hits[j][ii] = -100;
      tree_tdc_le_hits[j][ii] = -100;
    }
  }

  // dump event?
  if (dodump) {
    cout << "\n\n Event number " << dec << evnum;
    cout << " length " << evlen << " type " << evtype << endl;
    int ipt = 0;
    for (int j = 0; j < (evlen/5); j++) {
      cout << dec << "\n evbuffer[" << ipt << "] = ";
      for (int k=j; k<j+5; k++) {
        cout << hex << data[ipt++] << " ";
      }
      cout << endl;
    }
    if (ipt < evlen) {
      cout << dec << "\n evbuffer[" << ipt << "] = ";
      for (int k=ipt; k<evlen; k++) {
        cout << hex << data[ipt++] << " ";
      }
      cout << endl;
    }
  }

  //printf("#localTrig = %d \n",numLocal);

  // ----- loop to find header of ADC block read -----
  int index=0;
  int indexlast=0;
  if (debug) cout << " START LOOKING AT  at index = " << index << " last = "<< indexlast << " " << roclen[myroc]<< endl;
  while (((data[index]&0xffffffff)!=0x0da000011)) {
    index++;
  }
  while (((data[indexlast]&0xffffffff)!=0x0da000022)) {
    indexlast++;
  }
  if (debug) cout << " ADC header at index = " << index << " last = "<< indexlast << endl;

  // ----- loop to find header of TDC block read. Start at indexlast, end at indexend -----
  int indexend=indexlast;
  while (((data[indexend]&0xffffffff)!=0x0da000033)) {
    indexend++;
  }
  if (debug) cout << " TDC header at index = " << indexlast << " last = "<< indexend << endl;

  // ----- main ADC loop -----
  slotnew = 0;
  for (int j = index+1; j < indexlast; j++) {
    if (debug) printf("data[%d] = 0x%x = (dec)  \n",j,data[j]);

    int slot = (data[j]&0xf8000000)>>27; 
    int slot_ndata = 0;
    Int_t slotindex=0;

    if (slot!=slotnew) {
      slotindex=0;
      slotnew=slot;
      slot_ndata= (data[j]&0x7f) - 1;
      if (debug) cout << slot_ndata << " " << slot << endl;
      for (int jj=0;jj<NUMADCSLOTS;jj++) {
        if (slot==slotindADC[jj]) slotindex=jj;
      }
      if (slot_ndata > NUMCHANA) {
        printf("*** ADC bad evnum= %d,slot_ndata = %d, data[%d] = 0x%x = (dec)  \n",evnum,slot_ndata,j,data[j]);
        nbad_adc[slotindex]++;
        return;
      }
    }
    if (debug) cout << "ADC slot "<<slot<< " slot index = " << slotindex << " slot indata = " << slot_ndata << endl;

    tree_adc_slot = slotindex;
    //
    for (int ii=0;ii<slot_ndata;ii++){
      j++;
      if (debug==2) printf("data[%d] = 0x%x = (dec) %d %d  \n",j,data[j],(data[j]&0xfe0000)>>17,data[j]&0x3fff);
      ichan = (data[j]&0xfe0000)>>17; // 1881
      rdata = data[j]&0x3fff;  // 1881
      //	   ichan = (data[j]&0xfe0000)>>17; // 1877
      //rdata = data[j]&0xffff;  // 1877
      if (ichan >= 0 && ichan < NUMCHANA && rdata >= 0) {
        adcdat[slotindex][ii] = rdata; 
        adcchan[slotindex][ii] = ichan; 

        tree_adc_chan = ichan;
        tree_adc_hits[tree_adc_slot][tree_adc_chan] = rdata;
	tree_adc_pedestal_chan[tree_adc_slot][tree_adc_chan] = pedestal_threshold[tree_adc_slot][tree_adc_chan];

	tree_adc_hits_minus_pedestal[tree_adc_slot][tree_adc_chan] = rdata - pedestal_threshold[tree_adc_slot][tree_adc_chan];

        //	cout << "************* CHECK PEDESTAL slot "<< slot<<" slotindex " << slotindex << " chan "<< tree_adc_chan <<" tree_adc_hit " << tree_adc_hits <<" pedestal_chan " << 	tree_adc_pedestal_chan << endl;
      }
    }
    numslothitsADC[slotindex]=slot_ndata;
  } // end of ADC loop

  Int_t nbadevents_counter = 0;
  for (int i=0; i<NUMADCSLOTS; i++)  nbadevents_counter += numslothitsTDC[i];  
  if (nbadevents_counter > 0)
    nbadevents_adc++;

  // ----- main TDC loop -----

  //adapted from ADC loop
  slotnew=0;

  is_good = false; // overall flag for the event
  for (int j = indexlast+1; j < indexend; j++) {
    if (debug) printf("data[%d] = 0x%x = (dec)  \n",j,data[j]);

    int slot = (data[j]&0xf8000000)>>27; 
    int slot_ndata = 0;
    Int_t slotindex=0;

    if (slot!=slotnew) {
      slotindex=0;
      slotnew=slot;
      //      slot_ndata= (data[j]&0x7f) - 1;
      slot_ndata= (data[j]&0x7ff) - 1;
      if (debug) cout << slot_ndata << " " << slot << endl;
      for (int jj=0;jj<NUMTDCSLOTS;jj++) {
        if (slot==slotindTDC[jj]) slotindex=jj;
      }
      if (slot_ndata > NUMCHANT) {
        printf("*** TDC bad evnum= %d,slot_ndata = %d, data[%d] = 0x%x = (dec)  \n",evnum,slot_ndata,j,data[j]);
        nbad_tdc[slotindex]++;
        return;
      }
    }
    if (debug) cout << "TDC slot "<<slot<< " slot index = " << slotindex << " slot indata = " << slot_ndata << endl;

    tree_tdc_slot = slotindex;

    bool has_leading = false;
    bool has_trailing = false;
    bool more_than_1 = false;

    numslothitsTDC[slotindex]=slot_ndata;

    if (slot_ndata > 0){ // we have TDC data!
      is_good = true;
      for (int ii=0;ii<slot_ndata;ii++){
	j++;
	if (debug==2) printf("data[%d] = 0x%x = (dec) %d %d  \n",j,data[j],(data[j]&0xfe0000)>>17,data[j]&0xffff);
      ichan = (data[j]&0xfe0000)>>17; // 1877
      rdata = data[j]&0xffff;  // 1877
      if (ichan >= 0 && ichan < NUMCHANT && rdata >= 0) {
        tree_tdc_chan = ichan;
        //tree_tdc_hits[tree_tdc_slot][tree_tdc_chan] = rdata;
	edgetype = (data[j]&0x10000)>>16; // the type of the edge is in bit 16
	// if (rdata >= tree_tdc_hits[tree_tdc_slot][tree_tdc_chan] && edgetype==TRAILING_EDGE) tree_tdc_hits[tree_tdc_slot][tree_tdc_chan] = rdata;
	// if (rdata >= tree_tdc_le_hits[tree_tdc_slot][tree_tdc_chan] && edgetype==LEADING_EDGE) tree_tdc_le_hits[tree_tdc_slot][tree_tdc_chan] = rdata;

	if (edgetype == LEADING_EDGE)
	  tdcldat[slotindex][ichan] = rdata;

	if (edgetype == TRAILING_EDGE)
	  tdctdat[slotindex][ichan] = rdata;
	
   	// check if event is good
  	if (edgetype == LEADING_EDGE) {
  	  if (!has_leading) {
  	    // first leading edge encountered for this event
  	    has_leading = true;
  	  } else {
  	    // what?!? more than one leading edge?
  	    more_than_1 = true;
  	  }
  	}
  	if (edgetype == TRAILING_EDGE) {
  	  if (!has_trailing) {
  	    // first leading edge encountered for this event
  	    has_trailing = true;
  	  } else {
  	    // what?!? more than one trailing edge?
  	    more_than_1 = true;
  	  }
  	}

	//is_good = (is_good) && (has_leading) && (has_trailing);// && (!more_than_1);
	
        tree_tdc_nhits[tree_tdc_slot][tree_tdc_chan]++;
	// tdcdat[edgetype][slotindex][ii] = rdata; 
        tdcchan[slotindex][ii] = ichan; 
      } // matches if (ichan >= 0 && ichan < NUMCHANT && rdata >= 0) {


      } // matches for (int ii=0;ii<slot_ndata;ii++){
    }// matches if (slot_ndata > 0){

  } // end of TDC loop, matches for (int j = indexlast+1; j < indexend; j++) 

  
  
  // //old
  // is_good = true;
  // //
  // for (int j = indexlast+1; j < indexend; ) { // j will be manually updated
  //   if (debug) printf("TDC data[%d] = 0x%x = (dec)  \n",j,data[j]);

  //   int slot = (data[j]&0xf8000000)>>27; 
  //   int slot_ndata = 0;
  //   Int_t slotindex=0;
 
  //   if (slot!=slotnew) {
  //     slotindex=0;
  //     slotnew=slot;
  //     slot_ndata= (data[j]&0xff) - 1;

  //     if (debug) cout << slot_ndata << " " << slot << endl;
  //     for (int jj=0;jj<NUMTDCSLOTS;jj++) {
  //       if (slot==slotindTDC[jj]) slotindex=jj;
  //     }
  //     if (slot_ndata > NUMCHANT) {
  //       printf("***TDC bad evnum= %d,slot_ndata = %d, data[%d] = 0x%x = (dec)  \n",evnum,slot_ndata,j,data[j]);
  //       nbad_tdc[slotindex]++;
  //       return;
  //     }
  //   }
  //   if (debug) cout << "TDC slot "<<slot<< " slot index = " << slotindex << " slot indata = " << slot_ndata << endl;
    
  //   int how_many_datawords = slot_ndata;

  //   tree_tdc_slot = slotindex;

  //   bool has_leading = false;
  //   bool has_trailing = false;
  //   bool more_than_1 = false;
  //   int oldchan = (data[j+1]&0xfe0000)>>17; // init as the first chan
 
  //   for (int ii=0;ii<how_many_datawords;ii++){
  //     j++;
  //     ichan = (data[j]&0xfe0000)>>17; // 1881
  //     rdata = data[j]&0xffff;  // 1881
  //     //	   ichan = (data[j]&0xfe0000)>>17; // 1877
  //     //rdata = data[j]&0xffff;  // 1877
  //     edgetype = (data[j]&0x10000)>>16; // the type of the edge is in bit 16
  //     if (debug==2) printf("TDC data[%d] = 0x%x = (dec) %d %d %d \n",j,data[j],(data[j]&0xfe0000)>>17,data[j]&0x3fff,edgetype);

  //     if (ichan != oldchan) {
  // 	oldchan = ichan;
  // 	has_leading = false;
  // 	has_trailing = false;
  // 	more_than_1 = false;
  //     }
     
  //     if (ichan >= 0 && ichan < NUMCHANT && rdata >= 0) {
  //       tdcdat[edgetype][slotindex][ii] = rdata; 
  //       tdcchan[slotindex][ii] = ichan; 

  //       tree_tdc_chan = ichan;
  // 	// trailing edge is one furtherest from the common stop, reverse the normal idea of trailing and leading edge
  // 	if (edgetype==TRAILING_EDGE)tree_tdc_hits[tree_tdc_slot][tree_tdc_chan] = rdata;
  // 	if (edgetype==LEADING_EDGE)tree_tdc_le_hits[tree_tdc_slot][tree_tdc_chan] = rdata;
  // 	//tree_tdc_hits[tree_tdc_slot][tree_tdc_chan] = rdata;
  //       isleadingedge = edgetype;
  // 	// check if event is good
  // 	if (edgetype == LEADING_EDGE) {
  // 	  if (!has_leading) {
  // 	    // first leading edge encountered for this event
  // 	    has_leading = true;
  // 	  } else {
  // 	    // what?!? more than one leading edge?
  // 	    more_than_1 = true;
  // 	  }
  // 	}
  // 	if (edgetype == TRAILING_EDGE) {
  // 	  if (!has_trailing) {
  // 	    // first leading edge encountered for this event
  // 	    has_trailing = true;
  // 	  } else {
  // 	    // what?!? more than one trailing edge?
  // 	    more_than_1 = true;
  // 	  }
  // 	}
  //     }
  //     numslothitsTDC[slotindex]=slot_ndata;
  //     is_good = (is_good) && (has_leading) && (has_trailing) && (!more_than_1);
  //   }
  //   j++; // manual update after all datawords have been read
  // } // end of TDC loop
  
    // if (tree_tdc_hits[0][1] >= -100. ) cout << tree_tdc_hits[0][1] << " ADC data = " << tree_adc_hits[0][1] << " " << pedestal_threshold[0][1] << endl;
  
    // ----- fill -----
  // force event to be good?
  // is_good = true;

  // if (is_good) {
  //   tree->Fill();
  // } else {
  //    if (debug) cout << "*** NO associated TDC data, only slot headers. Skipping event "<< evnum <<endl;
  // }



  // check on first two groups
  // Int_t nbadevents_tdc_counter = 0;
  // for (int i=0; i<NUMTDCSLOTS; i++)  nbadevents_tdc_counter += numslothitsTDC[i];  
  // if (nbadevents_tdc_counter == 0)
  //   nbadevents_tdc++;

  // for (int pixel = 1; pixel<17; pixel++){
  //   a[1][pixel] = tree_adc_hits_minus_pedestal[handmapping_adc_slot(1)][handmapping_adc_chan(1, pixel)];
  //   a[2][pixel] = tree_adc_hits_minus_pedestal[handmapping_adc_slot(2)][handmapping_adc_chan(2, pixel)];

  //   t[1][pixel] = tree_tdc_hits[handmapping_tdc_slot(1)][handmapping_tdc_chan(1, pixel)];
  //   t[2][pixel] = tree_tdc_hits[handmapping_tdc_slot(2)][handmapping_tdc_chan(2, pixel)];
  // }

  // if (1==-1) {
  //   cout <<"Ev. "<< evnum << " " << endl;; 
  // if (adcdat[0][0 ] >  pedestal_threshold[0][0]+30) {
  //   printf(" top scin ADC   %4d %4d \n", pedestal_threshold[0][0],adcdat[0][0]);
  // }
  // if (adcdat[0][1 ] >  pedestal_threshold[0][1]+30) {
  //   printf("bottom scint ADC  %4d %4d \n", pedestal_threshold[0][1],adcdat[0][1]);
  // }
  // if (tdcldat[2][48 ] >  0) {
  //   printf("top scint TDC  %4d %4d\n",tdcldat[2][48],tdctdat[2][48]);
  // }
  // if (tdcldat[2][49 ] >  0) {
  //   printf("bottom scint TDC  %4d %4d \n",tdcldat[2][49],tdctdat[2][49]);
  // }
  // }
  // // quick display of ADC values for group 1
  // Double_t test;
  // test=0;
  // quick=1;

  // int slot_a = 1;
  // int chan_a = 2;
  // int slot_t = 0;
  // int chan_t = 66;
  // int adc_more_than = 30;

  // quick = 0;
  // if (quick == 1){
  //   for (int slot=0; slot<4; slot++) {
  //     for (int chan=0; chan<64; chan++) {       
  // 	if (adcdat[slot][ chan ] >  pedestal_threshold[slot][ chan ]+ adc_more_than && !(slot==1 && chan >15 && chan < 60)) {
  // 	  if ( //slot == slot_a && chan ==chan_a && 
  // 	      1) {
  // 	    cout <<"Ev. "<< evnum << " "; 
  // 	    //cout<<"ADC-a ";
  // 	    //cout<<endl;
  // 	    printf("ADC  %4d %4d %4d %4d match TDC slot %4d chan %4d", slot,chan, pedestal_threshold[slot][ chan ],adcdat[slot][ chan ],match_tdc_slot_for_adc_slot_chan(slot,chan),match_tdc_chan_for_adc_slot_chan(slot,chan));
  // 	    test=1;
  // 	    cout<<endl;
  // 	  }
  // 	} // matches "if (adcdat[slot][ chan ] >  pedestal_threshold[slot][ chan ]..."
  //     }
  //   }
  //   //
  // 	  if (tdcldat[2][48] == 0) {
  // 	       cout <<"Problem tdc Ev. "<< evnum << " ";
  // 	  }
  //    //
  //   if ( test ==1) {
  //     cout<<"TDC-a ";
  //     cout<<endl;
      
  //     for (int slott=0; slott<3; slott++) {
  // 	for (int chant=0; chant<96; chant++) {         
  // 	  if (tdcldat[slott][ chant ] > 0) {
  // 	    printf(" %4d %4d %4d %4d", slott,chant,tdcldat[slott][ chant ],tdctdat[slott][ chant ]);
  // 	    if (slott == slot_t && chant == chan_t) cout<<"                 <==";
  // 	    cout << endl;
  // 	  }
  // 	}
  //     }
  //   }
  // }


  /*
  cout<<"TDC-l ";
  for (int chan=65; chan<65+16; chan++)
    printf(" %4d", tdcldat[0][ chan ]);
  cout<<endl;
  cout<<"TDC-t ";
  for (int chan=65; chan<65+16; chan++)
    printf(" %4d", tdctdat[0][ chan ]);
  cout<<endl;
  }
  */

  // mapping
  Int_t slot, chan;
  for (int pmt = 1; pmt < NPMT+1; pmt++){
    for (int pixel = 1; pixel < NPIXEL+1; pixel++){
      slot = handmapping_adc_slot(pmt);
      chan = handmapping_adc_chan(pmt, pixel);
      a[pmt][pixel] = adcdat[slot][chan];

      slot = handmapping_tdc_slot(pmt);
      chan = handmapping_tdc_chan(pmt, pixel);
      t[pmt][pixel] = tdctdat[slot][chan];
      l[pmt][pixel] = tdcldat[slot][chan];
    }
  }


  tree->Fill();
} // end of decode()

void analysis() {
  Int_t islot, ichan,ihit, rdata;
  Int_t rawtimes[2];
  if (debug) cout << " analysis " << endl;

  // ADC
  //if (numLocal >= 1){
  for (islot = 0; islot < NUMADCSLOTS; islot++) {
    hadc1->Fill(slotindADC[islot] ,numslothitsADC[islot]);

    if (debug) cout << " slot " << slotindADC[islot] << " " << numslothitsADC[islot] << endl;
    for (ihit = 0; ihit < numslothitsADC[islot] ; ihit++) {
      rdata = adcdat[islot][ihit];
      ichan = adcchan[islot][ihit];
      hbranch->Fill(branchnum);

      if (ichan == mychan[0] &&  slotindADC[islot]== myslot[0]) hadc2[0]->Fill(rdata);
      if (ichan == mychan[1] &&  slotindADC[islot]== myslot[1]) hadc2[1]->Fill(rdata);
      if (ichan == mychan[2] &&  slotindADC[islot]== myslot[2]) hadc2[2]->Fill(rdata);
      if (ichan == mychan[3] &&  slotindADC[islot]== myslot[3]) hadc2[3]->Fill(rdata);

      hadc3[islot]->Fill(rdata,float(ichan)); // all channels of this slot
      hadcped[islot][ichan]->Fill(rdata); // for reading peadestals
    }
  } // end of ADC

    // TDC
    //   if (numLocal >= 1){
  Int_t ngoodpulses = 0; // stores the number of hit couples, that is: a leading edge and a trailing edge in successive hits
  for (islot = 0; islot < NUMTDCSLOTS; islot++) {
    if (debug) cout << "TDC  analysis: islot " << islot << " slotindTDC[islot] " << slotindTDC[islot] << " numslothitsTDC[islot] " << numslothitsTDC[islot] << endl;
    htdc1->Fill(slotindTDC[islot] ,numslothitsTDC[islot]);
    ihit=0;
    while (ihit < numslothitsTDC[islot] ) { 
      if (tdcldat[islot][ihit] > 0 && tdctdat[islot][ihit+1] > 0) {
	ngoodpulses++;
	rawtimes[TRAILING_EDGE] = tdctdat[islot][ihit+1];
	rawtimes[LEADING_EDGE]  = tdcldat[islot][ihit];
	ichan = tdcchan[islot][ihit];
	if (debug) cout << " TDC (before filling) slot " << slotindTDC[islot] << " " << numslothitsTDC[islot] << " " << ichan << " " << rawtimes[0]<< " " << rawtimes[1]<< " " << rawtimes[TRAILING_EDGE] - rawtimes[LEADING_EDGE] << endl;
	hleadingtime->Fill(float(rawtimes[LEADING_EDGE]));
	htdcpulsewidth->Fill(float(rawtimes[TRAILING_EDGE] - rawtimes[LEADING_EDGE]));
	htdc[islot][ichan]->Fill(rawtimes[LEADING_EDGE]); // for reading pedestals
      }
      ihit++;
    }

    if (debug) cout << "    TDC analysis: event = "<< ievent <<" slot = " << islot << " numslothitsTDC = " << numslothitsTDC[islot] <<", found " << ngoodpulses << " good pulses " << endl;
  } // end of TDC
    //	}
}

// calculate and save pedestal channels
void pedsup(){
  ofstream outfile("ped_test.dat");
  Double_t maxbin, pedx;
  Int_t pedend, content;
  Int_t thres = 0;
  for (Int_t i = 0; i < NUMADCSLOTS;i++) {
    outfile << "slot=" << slotindADC[i] << endl;
    for (Int_t j = 0; j < NUMCHANA;j++) {
      maxbin =  hadcped[i][j]->GetMaximumBin();
      pedx =  hadcped[i][j]->GetBinCenter(maxbin);

      //thres = pedx;   
      // OR
      int maxx = pedx + 50;
      int minx = pedx - 30; 
      TF1* gausfit=new TF1("gausfit","gaus",0,MAXPED);
      hadcped[i][j]->Fit("gausfit","","",minx,maxx);
      float mean  = gausfit->GetParameter(1);
      //float sigma = gausfit->GetParameter(2);
      thres = int(mean);
   
      pedend = 0;
      for (Int_t ichan=MAXPED; ichan>thres; ichan--){
	content = hadcped[i][j]->GetBinContent(ichan);
	if (content > 0 && pedend == 0) pedend = ichan;
      }
      outfile <<  thres;
      //    outfile << " " << pedend;
      outfile << endl;
    }
  }
  outfile.close();
}

void old_pedsup(){
  ofstream outfile("ped_test.dat");
  Double_t maxbin, pedx, maxx, minx, mean = 0, sigma=0; 
  Int_t thres = 0;
  TF1* gausfit=new TF1("gausfit","gaus",0,2000);
  for (int i = 0; i < NUMADCSLOTS;i++) {
    outfile << "slot=" << slotindADC[i] << endl;
    for (int j = 0; j < NUMCHANA;j++) {
      maxbin =  hadcped[i][j]->GetMaximumBin();
      pedx =  hadcped[i][j]->GetBinCenter(maxbin);
      maxx = pedx + 50;
      minx = pedx - 50; 
      hadcped[i][j]->Fit("gausfit","","",minx,maxx);
      mean  = gausfit->GetParameter(1);
      sigma = gausfit->GetParameter(2);
      thres = int(mean + 5*sigma);
      outfile <<  thres << endl;

      cout<<"slot "<<slotindADC[i]<<" index "<< j <<" mean,sigma,thres "<< mean <<" "<< sigma <<" "<< thres <<endl;

    }
  }
  outfile.close();
}

// read out pedestal channels from file
void readpedsup(){
  ifstream infile("ped_test.dat");
  //char thres[4], endp[4];
  char thres[4];
  Int_t slot = 0;
  char header[8]; // will contain "slot=.." where .. are the slot digits
  char slotstr[2]; // helps parsing header
  debug=0;
  for (int i = 0; i < NUMADCSLOTS;i++) {
    infile >> header; // slot=..
    //cout << "readpedsup(): read line: " << header << endl;
    slotstr[0] = header[5]; // first slot digit
    slotstr[1] = header[6]; // second slot digit
    slot = atoi(slotstr);   // convert header substring (saved in slotstr) to slot#
    cout << "readpedsup(): reading pedestals of ADC slot " <<slot << "..." << endl;
    for (int j = 0; j < NUMCHANA;j++) {
      infile >>  thres;// >> endp;
      pedestal_threshold[i][j] = atoi(thres);
      if (debug) cout << "readpedsup(): chan = "<< j<<" thres = '" << thres << "' value = " << pedestal_threshold[i][j] <<endl;

      if ( debug && (pedestal_threshold[i][j] < 1 || pedestal_threshold[i][j] > MAXPED) ) cout << "WARNING: pedestal of physical slot "<< slot <<" channel "<< j <<" is "<< pedestal_threshold[i][j] <<" which is out of boundary!"<<endl;

    }
  }
  debug=0;
  infile.close();
  cout << "readpedsup(): done! "<<endl;
}

Int_t match_tdc_slot_for_adc_slot_chan(Int_t aslot, Int_t achan){
  if ( aslot == 0 && achan <64) return 0;
  if ( aslot == 1 && achan >=0 && achan < 16) return 0;
  if ( aslot == 1 && achan >=48 && achan < 64) return 1;
  if ( aslot == 1 && achan >=16 && achan < 48) return -1;
  if ( aslot == 2 && achan >=0 && achan < 64) return 1;
  if ( aslot == 3 && achan >=16 && achan < 32) return 2;
  if ( aslot == 3 && achan >=32 && achan < 48) return 0;
  if ( aslot == 3 && achan >=48 && achan < 64) return 1;
  return -1;
}
Int_t match_tdc_chan_for_adc_slot_chan(Int_t aslot, Int_t achan){
  if ( aslot == 0 && achan <64) return achan;
  if ( aslot == 1 && achan >=0 && achan < 16) return 0;
  if ( aslot == 1 && achan >=48 && achan < 64) return achan+64;
  if ( aslot == 1 && achan >=16 && achan < 48) return achan-32;
  if ( aslot == 2 && achan >=0 && achan < 64) return achan+32;
  if ( aslot == 3 && achan >=16 && achan < 32) return achan;
  if ( aslot == 3 && achan >=32 && achan < 48) return achan+48;
  if ( aslot == 3 && achan >=48 && achan < 64) return achan-48;
  return -1;
}
// Int_t handmapping_adc_slot(Int_t group){
//   switch (group){
//   case 1:
//   case 2:
//   case 3:
//   case 4:
//     return 0;
//     break;
//   case 5:
//   case 8:
//     return 1;
//     break;
//   case 6:
//   case 7:
//   case 13:
//   case 14:
//     return 3;
//     break;
//   case 9:
//   case 10:
//   case 11:
//   case 12:
//     return 2;
//     break;
//   default:
//     return -1;
//   }
// }


// Int_t handmapping_adc_chan(Int_t group, Int_t pixel){
//   pixel--;
//   switch(group){
//   case 1:
//   case 5:
//   case 9:
//   case 13:
//     return pixel;
//     break;
//   case 2:
//   case 10:
//   case 14:
//     return 16+pixel;
//     break;
//   case 3:
//   case 6:
//   case 11:
//     return 32+pixel;
//     break;
//   case 4:
//   case 7:
//   case 8:
//   case 12:
//     return 48+pixel;
//     break;
//   default:
//     return -1;
//   }
// }



// Int_t handmapping_tdc_slot(Int_t group){
//   switch (group){
//   case 1:
//   case 2:
//   case 3:
//   case 4:
//   case 5:
//   case 6:
//     return 0;
//     break;
//   case 7:
//   case 8:
//   case 9:
//   case 10:
//   case 11:
//   case 12:
//     return 1;
//     break;
//   case 13:
//   case 14:
//     return 2;
//     break;
//   default:
//     return -1;
//   }
// }



// Int_t handmapping_tdc_chan(Int_t group, Int_t pixel){
//   pixel--;
//   switch(group){
//   case 1:
//   case 7:
//   case 13:
//     return pixel;
//     break;
//   case 2:
//   case 8:
//   case 14:
//     return 16+pixel;
//     break;
//   case 3:
//   case 9:
//     return 32+pixel;
//     break;
//   case 4:
//   case 10:
//     return 48+pixel;
//     break;
//   case 5:
//   case 11:
//     return 64+pixel;
//     break;
//   case 6:
//   case 12:
//     return 80+pixel;
//   default:
//     return -1;
//   }
// }

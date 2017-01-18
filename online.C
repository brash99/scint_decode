// A simple decoder and histogrammer for fastbus data.
// R. Michaels, Oct 24, 2014.
// developed later by Dasuni, et al.

#define MAXROC     32
#define MAXHITS    100
#define MAXEVENTS  1000000
// constants used to distinguish which case the bit 16 of TDC data
// word is referring to. Note that these values are SWAPPED with
// respect to the TDC specifications, because we are working in common
// stop mode
#define TRAILING_EDGE 1
#define LEADING_EDGE  0

#define NUMCHANA   64 // ADC chans per slot
#define NUMCHANT   96 //*2*3 // (TDC chans per slot) * (leading+trailing edge) * (number of pulses recorded)

#define NUMADCSLOTS 4   // # of ADC slots
#define NUMTDCSLOTS 3   // # of TDC slots

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
void decode(Int_t* data, TTree *tree);
void analysis();
void clear();
void pedsup();
void readpedsup();
Int_t mapping_read(const char*);

// Global data 
// -----------
Int_t evlen, evtype, evnum;
Int_t *irn, *rocpos, *roclen;
Int_t myroc, myslot[4], mychan[4];
Int_t ievent = 0;
Int_t fHasHeader=1;   // =1 for 1877 and 1881, but =0 for 1875
Int_t fWordSeen;
Int_t debug=0;                     // *** DEBUG is here ***
Int_t numLocal = 0;
Int_t branchnum = 999;

//ADC
Int_t numslothitsADC[NUMADCSLOTS];
Int_t adcrawdat[NUMADCSLOTS][NUMCHANA];
Int_t adcchan[NUMADCSLOTS][NUMCHANA];
Int_t adcmped[NUMADCSLOTS][NUMCHANA];
//Int_t nbad_adc[NUMADCSLOTS]={0,0,0,0}; // for test crate. Should be fine for the other ones also
Int_t slotindADC[NUMADCSLOTS];
// Pedestal
Int_t PEDSUP; // 0 = read pedestal file, 1 = calculate and write file. Read from argv
Int_t ped_threshold[NUMADCSLOTS][NUMCHANA];
Int_t ped_ends[NUMADCSLOTS][NUMCHANA];

// TH1F  *hadc1;
// TH1F  *hadc2[4];
// TH2F  *hadc3[NUMADCSLOTS];
// TH2F  *h2adc;
TH1F  *hadcped[NUMADCSLOTS][NUMCHANA];
// TH1F  *htdc[NUMTDCSLOTS][NUMCHANT];
TH1F  *hbranch;

//TDC
Int_t numslothitsTDC[NUMTDCSLOTS];
Int_t tdctdat[NUMTDCSLOTS][NUMCHANT]; // trailing edge
Int_t tdcldat[NUMTDCSLOTS][NUMCHANT]; // leading edge
Int_t tdcchan[NUMTDCSLOTS][NUMCHANT];
//Int_t nbad_tdc[NUMTDCSLOTS]={0,0,0};
Int_t slotindTDC[NUMTDCSLOTS]={14,15,16};


// TH1F  *htdc1;
// TH1F  *hleadingtime;
// TH1F  *htdcpulsewidth;

// // MAPPING
// // We count pmt and pixel from 1, as opposed to 0, I prefer to just
// // add +1 to the constants here as opposed to every loop, and in the
// // definition of the branchs
// static const Int_t MAXPMTIDS = 1 + 400; // number of user ID's
// static const Int_t MAXPMTS = 1 + NUMADCSLOTS * NUMCHANA / 16; // max number of PMT's which can be read by DAQ
// static const Int_t NPIXELS = 1 + 16;
// Int_t map_nino_id[MAXPMTS];
// Int_t map_pmt_id[MAXPMTIDS]; // index = PMT "user ID"
// // pmt, pixel -> slot, chan
// Int_t map_adc_slot[MAXPMTS][NPIXELS];
// Int_t map_adc_chan[MAXPMTS][NPIXELS];
// Int_t map_tdc_slot[MAXPMTS][NPIXELS];
// Int_t map_tdc_chan[MAXPMTS][NPIXELS];
// Int_t map_adc_pmt[NUMADCSLOTS][NUMCHANA];
// Int_t map_adc_pixel[NUMADCSLOTS][NUMCHANA];
// Int_t map_tdc_pmt[NUMTDCSLOTS][NUMCHANT];
// Int_t map_tdc_pixel[NUMTDCSLOTS][NUMCHANT];
// Int_t map_len; // number of read pmt config lines
// const char *mapfile = "sbs_mapping.txt"; // config file

// Trees variables
// Int_t tree_adc_slot, tree_adc_chan;
// Int_t tree_adc_hits[NUMADCSLOTS][NUMCHANA];
// Int_t tree_adc_hits_m_ped[NUMADCSLOTS][NUMCHANA];
// Int_t tree_adc_pedestal_chan[NUMADCSLOTS][NUMCHANA];
// Int_t tree_tdc_slot, tree_tdc_chan;
// Int_t tree_tdc_hits[NUMTDCSLOTS][NUMCHANT];
// Int_t tree_tdc_nhits[NUMTDCSLOTS][NUMCHANT];
// Int_t tree_tdc_le_hits[NUMTDCSLOTS][NUMCHANT];
// Int_t isleadingedge;
// bool is_good; // true if event has both one leading and one trailing edge

// Int_t adcmapped[MAXPMTIDS][NPIXELS];
// Int_t adcrawmapped[MAXPMTIDS][NPIXELS];
// Int_t tdctmapped[MAXPMTIDS][NPIXELS];
// Int_t tdclmapped[MAXPMTIDS][NPIXELS];
// -----------

Int_t main(Int_t argc, char* argv[])
{
  //============ Initializations ============
  THaCodaData *coda;      
  char ctitle[100];
  char dtitle[100];
  char rtitle[100];
  //char htitle[100];
  char cpedtitle[100];
  char hpedtitle[100];
  Int_t maxevent = 10000;
  Int_t runno = 1000;
  Int_t istatus;
  Int_t lslot;

  Int_t choice1 = 1; // =1 CODA file, else ET connection
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

  // myslot[0]= lslot; myslot[1]= lslot+1; myslot[2]= lslot+2; myslot[3]= lslot+3;
  // mychan[0] = 0;  mychan[1] = 0; mychan[2] = 0; mychan[3] = 0; // defaults

  cout << "my roc is " << myroc << endl;

  for (Int_t i = 0; i < NUMADCSLOTS; i++){
    slotindADC[i] = lslot;
    lslot++;
  }

  cout << "Fastbus analysis "<<endl;
  cout << "Events to process "<<maxevent<<endl;

  // Int_t nPMTs = mapping_read(mapfile);
  // cout << "Mapping read from "<< mapfile <<endl;
  // cout << " Found "<< nPMTs <<" PMT config lines"<< endl;

  // ----- init root -----
  sprintf(rtitle,"sbs_%d_%i.root",runno,myroc);
  sprintf(dtitle,"../data/scint_%d.dat",runno);
  TROOT fbana("fbroot","Hall A SCINT analysis for Fastbus");
  TFile hfile(rtitle,"RECREATE","SBS data");

  char tempstr[16];

  TTree *tree = new TTree("t","ADC & TDC data tree");
  //  tree->Branch("event"	    ,&ievent,        "event/I");
  
  sprintf(tempstr,"value[%d][%d]/I" ,NUMADCSLOTS, NUMCHANA);
  tree->Branch("adcraw"	            ,adcrawdat, tempstr); // adc raw values
  tree->Branch("adc"	            ,adcmped, tempstr); // adc values - pedestal

  // sprintf(tempstr,"pedestal[%d][%d]/I", NUMADCSLOTS, NUMCHANA);
  // tree->Branch("ped"              ,tree_adc_pedestal_chan,tempstr); //pedestal

  sprintf(tempstr,"value[%d][%d]/I",NUMTDCSLOTS, NUMCHANT);
  // tree->Branch("tdcn"		   ,tree_tdc_nhits, tempstr); // tdc number of leading and trailing edges
  tree->Branch("tdct"		   ,tdctdat, tempstr); // tdc trailing edge
  tree->Branch("tdcl"         	   ,tdcldat, tempstr); //tdc leading edge
  // tree->Branch("is_good"           ,&is_good,"is_good/B");

  // sprintf(tempstr,"value[%d][%d]/I" ,MAXPMTS, NPIXELS);
  // tree->Branch("a",adcmapped,tempstr);
  // tree->Branch("ar",adcrawmapped,tempstr);
  // tree->Branch("t",tdctmapped,tempstr);
  // tree->Branch("l",tdclmapped,tempstr);

  // ----- init output -----
  // ADC
  // sprintf(ctitle,"ADC hits per slot ");
  // hadc1 = new TH1F("hadc1",ctitle,31,-1,30.);
  // sprintf(ctitle,"ADC data on slot %d channel %d",myslot[0],mychan[0]);
  // hadc2[0] = new TH1F("hadc2_0",ctitle,200,0,4000.);
  // sprintf(ctitle,"ADC data on slot %d channel %d",myslot[1],mychan[1]);
  // hadc2[1] = new TH1F("hadc2_1",ctitle,200,0,4000.);
  // sprintf(ctitle,"ADC data on slot %d channel %d",myslot[2],mychan[2]);
  // hadc2[2] = new TH1F("hadc2_2",ctitle,200,0,4000.);
  // sprintf(ctitle,"ADC data on slot %d channel %d",myslot[3],mychan[3]);
  // hadc2[3] = new TH1F("hadc2_3",ctitle,200,0,4000.);

  // sprintf(ctitle,"ADC data on 2 channnels");
  // h2adc = new TH2F("h2adc",ctitle,200,0,4000.,200,0,4000.);
  // sprintf(ctitle,"Branch #");
  hbranch = new TH1F("hbranch",ctitle,10,0,10.);

  // for (Int_t i = 0; i < NUMADCSLOTS;i++) {
  //   sprintf(ctitle,"ADC data on all channels of slot %d",slotindADC[i]);
  //   sprintf(htitle,"ADC_chan_slot_%d",slotindADC[i]);
  //   hadc3[i] = new TH2F(htitle,ctitle,200,0,2000.,64,0,64.);
  // }
  for (Int_t i = 0; i < NUMADCSLOTS;i++) {
    for (Int_t j = 0; j < NUMCHANA;j++) {
      sprintf(cpedtitle,"myADC data on chan %d of slot %d", j, slotindADC[i]);
      sprintf(hpedtitle,"myADC_chan_%d_slot_%d",j,slotindADC[i]);
      hadcped[i][j] = new TH1F(hpedtitle,cpedtitle,1000,0,1000.);
    }
  }
  // for (Int_t i = 0; i < NUMTDCSLOTS;i++) {
  //   for (Int_t j = 0; j < NUMCHANT;j++) {
  //     sprintf(cpedtitle,"TDC data on chan %d of slot %d", j, slotindTDC[i]);
  //     sprintf(hpedtitle,"TDC_chan_%d_slot_%d",j,slotindTDC[i]);
  //     htdc[i][j] = new TH1F(hpedtitle,cpedtitle,1000,0,1000.);
  //   }
  // }

  // TDC
  // sprintf(ctitle,"TDC hits per slot ");
  // htdc1 = new TH1F("htdc1",ctitle,31,-1,30.);
  // htdc1->GetXaxis()->SetTitle("slot");
  // htdc1->GetYaxis()->SetTitle("counts");
  // sprintf(ctitle,"TDC leading edge times ");
  // hleadingtime = new TH1F("hleadingtime",ctitle,1000,0,1000.);
  // hleadingtime->GetXaxis()->SetTitle("channel (0.5 ns each)");
  // hleadingtime->GetYaxis()->SetTitle("counts");
  // sprintf(ctitle,"TDC pulse width ");
  // htdcpulsewidth = new TH1F("htdcpulsewidth",ctitle,1000,0,1000.);
  // htdcpulsewidth->GetXaxis()->SetTitle("channel (0.5 ns each)");
  // htdcpulsewidth->GetYaxis()->SetTitle("counts");

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
    Int_t mymode = 1;
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




  // // seek the number of events
  // int handle;
  // int b_handle;
  // //  int status = evOpen((char*)filename.Data(),"r",&handle); 
  // int status;

  // status = evOpen(dtitle,"r",&handle); 
  // if (status<0) {
  //   //    printf("Unable to open file %s, status = 0x%x\n",filename.Data(),status);
  //     printf("Unable to open file %s, status = 0x%x\n",dtitle,status);
  //     exit(-1);
  // }

  // int nevents;
  // nevents = evOpenSearch(handle,&b_handle);
  // evCloseSearch(b_handle);
  // evClose(handle);

  // //  cout << "Number of events found in the file " << filename.Data() << ": " << nevents << endl;
  // cout << "Number of events found in the file " << dtitle << ": " << nevents << endl;


  cout <<"Ready to go! starting event loop" << endl;
  // Loop over events

  for (Int_t iev = 0; iev < maxevent; iev++)  {//the loop over the events
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
  // for (Int_t i=0;i<NUMADCSLOTS;i++){
  //   printf("In Slot %d, Number of events with bad adc = %d \n",slotindADC[i],nbad_adc[i]);
  // }

  // TDC
  // for (Int_t i=0;i<NUMTDCSLOTS;i++){
  //   printf("In Slot %d, Number of events with bad tdc = %d \n",slotindTDC[i],nbad_tdc[i]);
  // }
  cout << endl;

  coda->codaClose();

  hfile.Write();
  hfile.Close();

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
  memset (adcrawdat, 0, NUMADCSLOTS*NUMCHANA*sizeof(Int_t));
  memset (adcmped, 0, NUMADCSLOTS*NUMCHANA*sizeof(Int_t));
  memset (adcchan, 0, NUMADCSLOTS*NUMCHANA*sizeof(Int_t));
  memset (numslothitsADC, 0, NUMADCSLOTS*sizeof(Int_t));

  //TDC
  memset (numslothitsTDC, 0, NUMTDCSLOTS*sizeof(Int_t));
  memset (tdctdat, 0, NUMTDCSLOTS*NUMCHANT*sizeof(Int_t));
  memset (tdcldat, 0, NUMTDCSLOTS*NUMCHANT*sizeof(Int_t));
  memset (tdcchan, 0, NUMTDCSLOTS*NUMCHANT*sizeof(Int_t));

  //MAPPING
  // memset(adcmapped, -1, MAXPMTIDS*NPIXELS*sizeof(Int_t));
  // memset(adcrawmapped, -1, MAXPMTIDS*NPIXELS*sizeof(Int_t));
  // memset(tdctmapped, -1, MAXPMTIDS*NPIXELS*sizeof(Int_t));
  // memset(tdclmapped, -1, MAXPMTIDS*NPIXELS*sizeof(Int_t));
}

//void decode (Int_t* data, TTree *tree_adc, TTree *tree_tdc) {
void decode (Int_t* data, TTree *tree) {
  // ----- init -----
  Int_t ichan = 0, rdata = 0;
  evlen = data[0] + 1;
  evtype = data[1]>>16;
  evnum = data[4];
  static Int_t dodump = 0;  // dump the raw data

  Int_t edgetype = 0;
  Int_t slotnew=0;

  if (evtype > 10) return;

  // dump event?
  if (dodump) {
    cout << "\n\n Event number " << dec << evnum;
    cout << " length " << evlen << " type " << evtype << endl;
    Int_t ipt = 0;
    for (Int_t j = 0; j < (evlen/5); j++) {
      cout << dec << "\n evbuffer[" << ipt << "] = ";
      for (Int_t k=j; k<j+5; k++) {
        cout << hex << data[ipt++] << " ";
      }
      cout << endl;
    }
    if (ipt < evlen) {
      cout << dec << "\n evbuffer[" << ipt << "] = ";
      for (Int_t k=ipt; k<evlen; k++) {
        cout << hex << data[ipt++] << " ";
      }
      cout << endl;
    }
  }

  //printf("#localTrig = %d \n",numLocal);

  // ----- loop to find header of ADC block read -----
  Int_t index=0;
  Int_t indexlast=0;
  if (debug) cout << "\n\n Event number " << dec << evnum << endl;
  if (debug) cout << " START LOOKING AT  at index = " << index << " last = "<< indexlast << " " << roclen[myroc]<< endl;
  while (((data[index]&0xffffffff)!=0x0da000011)) {
    index++;
  }
  while (((data[indexlast]&0xffffffff)!=0x0da000022)) {
    indexlast++;
  }
  if (debug) cout << " ADC header at index = " << index << " last = "<< indexlast << endl;

  // ----- loop to find header of TDC block read. Start at indexlast+1, end at indexend -----
  Int_t indexend=indexlast;
  while (((data[indexend]&0xffffffff)!=0x0da000033)) {
    indexend++;
  }
  if (debug) cout << " TDC header at index = " << indexlast << " last = "<< indexend << endl;
  // ----- main ADC loop -----

  if (debug) cout << "*** Starting ADC loop"<<endl;  

  slotnew = 0;
  for (Int_t j = index+1; j < indexlast; j++) {
    if (debug) printf("data[%d] = 0x%x = (dec)  \n",j,data[j]);

    Int_t slot = (data[j]&0xf8000000)>>27; 
    Int_t slot_ndata = 0;
    Int_t slotindex=0;

    if (slot!=slotnew) {
      slotindex=0;
      slotnew=slot;
      slot_ndata= (data[j]&0x7f) - 1;
      if (debug) cout << slot_ndata << " " << slot << endl;
      for (Int_t jj=0;jj<NUMADCSLOTS;jj++) {
        if (slot==slotindADC[jj]) slotindex=jj;
      }
      if (slot_ndata > NUMCHANA) {
        printf("*** ADC bad evnum= %d,slot_ndata = %d, data[%d] = 0x%x = (dec)  \n",evnum,slot_ndata,j,data[j]);
	//   nbad_adc[slotindex]++;
        return;
      }
    }
    if (debug) cout << "ADC slot "<<slot<< " slot index = " << slotindex << " slot indata = " << slot_ndata << endl;

    //    tree_adc_slot = slotindex;
    //
    for (Int_t ii=0;ii<slot_ndata;ii++){
      j++;
      if (debug==2) printf("data[%d] = 0x%x = (dec) %d %d  \n",j,data[j],(data[j]&0xfe0000)>>17,data[j]&0x3fff);
      ichan = (data[j]&0xfe0000)>>17; // 1881
      rdata = data[j]&0x3fff;  // 1881
      //	   ichan = (data[j]&0xfe0000)>>17; // 1877
      //rdata = data[j]&0xffff;  // 1877
      if (ichan >= 0 && ichan < NUMCHANA && rdata >= 0) {
        adcrawdat[slotindex][ii] = rdata; 
        adcchan[slotindex][ii] = ichan; 

        //tree_adc_chan = ichan;
	Int_t ped = ped_threshold[slotindex][ichan];

        //tree_adc_hits[tree_adc_slot][tree_adc_chan] = rdata;

	//tree_adc_ped_chan[tree_adc_slot][tree_adc_chan] = ped;
	adcmped[slotindex][ichan] = rdata - ped;
     
        //	cout << "************* CHECK PEDESTAL slot "<< slot<<" slotindex " << slotindex << " chan "<< tree_adc_chan <<" tree_adc_hit " << tree_adc_hits <<" pedestal_chan " << 	tree_adc_pedestal_chan << endl;

	// mapping
	// Int_t pmt = map_adc_pmt[ tree_adc_slot ][ tree_adc_chan ];
	// Int_t pixel = map_adc_pixel[ tree_adc_slot ][ tree_adc_chan ];
	// adcrawmapped[ pmt ][ pixel ] = rdata;
	// adcmapped[ pmt ][ pixel ] = rdata - ped;
      }
    }
    numslothitsADC[slotindex]=slot_ndata;
  } // end of ADC loop


  // ----- main TDC loop -----
  // adapted from ADC loop
  Int_t notdc = 0;

  if (debug) cout << "*** Starting TDC loop"<<endl;  

  slotnew = 0;
  for (Int_t j = indexlast+1; j < indexend; j++) {
    if (debug) printf("data[%d] = 0x%x = (dec)  \n",j,data[j]);

    Int_t slot = (data[j]&0xf8000000)>>27; 
    Int_t slot_ndata = 0;
    Int_t slotindex=0;

    if (slot!=slotnew) {
      slotindex=0;
      slotnew=slot;
      slot_ndata= (data[j]&0x7ff) - 1;
      if (debug) cout << slot_ndata << " " << slot << endl;
      for (Int_t jj=0;jj<NUMTDCSLOTS;jj++) {
        if (slot==slotindTDC[jj]) slotindex=jj;
      }
      if (slot_ndata > NUMCHANT) {
        printf("*** TDC bad evnum= %d,slot_ndata = %d, data[%d] = 0x%x = (dec)  \n",evnum,slot_ndata,j,data[j]);
	//   nbad_tdc[slotindex]++;
        return;
      }
    }
    if (debug) cout << "TDC slot "<<slot<< " slot index = " << slotindex << " slot indata = " << slot_ndata << endl;  

    //tree_tdc_slot = slotindex;

    for (Int_t ii=0;ii<slot_ndata;ii++){
      j++;
      if (debug==2) printf("data[%d] = 0x%x = (dec) %d %d  \n",j,data[j],(data[j]&0xfe0000)>>17,data[j]&0xffff);
      ichan = (data[j]&0xfe0000)>>17; // 1881
      rdata = data[j]&0xffff;  // 1881
      //	   ichan = (data[j]&0xfe0000)>>17; // 1877
      //rdata = data[j]&0xffff;  // 1877
      if (ichan >= 0 && ichan < NUMCHANT && rdata > 0) {
        tdcchan[slotindex][ii] = ichan; 
	edgetype = (data[j]&0x10000)>>16; // the type of the edge is in bit 16

        //tree_tdc_chan = ichan;

	if (edgetype == TRAILING_EDGE)
	  tdctdat[slotindex][ii] = rdata; 
	//tree_tdc_hits[tree_tdc_slot][tree_tdc_chan] = rdata;

	if (edgetype == LEADING_EDGE)
	  tdcldat[slotindex][ii] = rdata; 
	//tree_tdc_le_hits[tree_tdc_slot][tree_tdc_chan] = rdata;
     
	// mapping
	// Int_t pmt = map_tdc_pmt[ tree_tdc_slot ][ tree_tdc_chan ];
	// Int_t pixel = map_tdc_pixel[ tree_tdc_slot ][ tree_tdc_chan ];
	// if (edgetype == TRAILING_EDGE) tdctmapped[ pmt ][ pixel ] = rdata;
	// if (edgetype == LEADING_EDGE) tdclmapped[ pmt ][ pixel ] = rdata;

      }
      if (rdata == 0) {
	notdc++;
	cout<<"***WARNING: no TDC data for event "<< evnum<< endl; 
      }
    }
    numslothitsTDC[slotindex]=slot_ndata;
  } // end of TDC loop
  
    // ----- fill -----
  if (notdc > 0) {
    cout<<" Event skipped (no TDC data)"<<endl;
  } else {
    tree->Fill(); 
  }
} // end of decode()

void analysis() {
  Int_t islot, ichan,ihit, rdata;
  // Int_t rawtimes[2];
  // if (debug) cout << " analysis " << endl;

  // ADC
  // //if (numLocal >= 1){
  for (islot = 0; islot < NUMADCSLOTS; islot++) {
    //hadc1->Fill(slotindADC[islot] ,numslothitsADC[islot]);

    if (debug) cout << " slot " << slotindADC[islot] << " " << numslothitsADC[islot] << endl;
    for (ihit = 0; ihit < numslothitsADC[islot] ; ihit++) {
      rdata = adcrawdat[islot][ihit];
      ichan = adcchan[islot][ihit];
      hbranch->Fill(branchnum);

      // if (ichan == mychan[0] &&  slotindADC[islot]== myslot[0]) hadc2[0]->Fill(rdata);
      // if (ichan == mychan[1] &&  slotindADC[islot]== myslot[1]) hadc2[1]->Fill(rdata);
      // if (ichan == mychan[2] &&  slotindADC[islot]== myslot[2]) hadc2[2]->Fill(rdata);
      // if (ichan == mychan[3] &&  slotindADC[islot]== myslot[3]) hadc2[3]->Fill(rdata);

      // hadc3[islot]->Fill(rdata,float(ichan)); // all channels of this slot
      hadcped[islot][ichan]->Fill(rdata); // for reading peadestals
    }
  } // end of ADC

  //   // TDC
  //   //   if (numLocal >= 1){
  // Int_t ngoodpulses = 0; // stores the number of hit couples, that is: a leading edge and a trailing edge in successive hits
  // for (islot = 0; islot < NUMTDCSLOTS; islot++) {
  //   // if (debug) cout << "TDC  analysis: islot " << islot << " slotindTDC[islot] " << slotindTDC[islot] << " numslothitsTDC[islot] " << numslothitsTDC[islot] << endl;
  //   htdc1->Fill(slotindTDC[islot] ,numslothitsTDC[islot]);
  //   ihit=0;
  //   while (ihit < numslothitsTDC[islot] ) { 
  //     if (tdcdat[LEADING_EDGE][islot][ihit] > 0 && tdcdat[TRAILING_EDGE][islot][ihit+1] > 0) {
  // 	ngoodpulses++;
  // 	rawtimes[TRAILING_EDGE] = tdcdat[TRAILING_EDGE][islot][ihit+1];
  // 	rawtimes[LEADING_EDGE]  = tdcdat[LEADING_EDGE][islot][ihit];
  // 	ichan = tdcchan[islot][ihit];
  // 	if (debug) cout << " TDC (before filling) slot " << slotindTDC[islot] << " " << numslothitsTDC[islot] << " " << ichan << " " << rawtimes[0]<< " " << rawtimes[1]<< " " << rawtimes[TRAILING_EDGE] - rawtimes[LEADING_EDGE] << endl;
  // 	hleadingtime->Fill(float(rawtimes[LEADING_EDGE]));
  // 	htdcpulsewidth->Fill(float(rawtimes[TRAILING_EDGE] - rawtimes[LEADING_EDGE]));
  // 	htdc[islot][ichan]->Fill(rawtimes[LEADING_EDGE]); // for reading pedestals
  //     }
  //     ihit++;
  //   }

  //   // if (debug) cout << "    TDC analysis: event = "<< ievent <<" slot = " << islot << " numslothitsTDC = " << numslothitsTDC[islot] <<", found " << ngoodpulses << " good pulses " << endl;
  // } // end of TDC
  //   //	}
}

// calculate and save pedestal channels
void pedsup(){
  ofstream outfile("ped_test.dat");
  Double_t maxbin, pedx;
  //int maxx, minx, mean = 0, sigma=0;
  Int_t pedend, content;
  Int_t thres = 0;
  //TF1* gausfit=new TF1("gausfit","gaus",0,1000);
  for (Int_t i = 0; i < NUMADCSLOTS;i++) {
    outfile << "slot=" << slotindADC[i] << endl;
    for (Int_t j = 0; j < NUMCHANA;j++) {
      maxbin =  hadcped[i][j]->GetMaximumBin();
      pedx =  hadcped[i][j]->GetBinCenter(maxbin);
      // maxx = pedx + 50;
      // minx = pedx - 30; 
      // hadcped[i][j]->Fit("gausfit","","",minx,maxx);
      // mean  = gausfit->GetParameter(1);
      // sigma = gausfit->GetParameter(2);
      // thres = int(mean);
      thres = pedx;
      pedend = 0;
      for (Int_t ichan=1000; ichan>thres; ichan--){
	content = hadcped[i][j]->GetBinContent(ichan);
	if (content > 0 && pedend == 0) pedend = ichan;
      }
      outfile <<  thres;
      outfile << " " << pedend;
      outfile << endl;
    }
  }
  outfile.close();
}

// void pedsup_old(){
//   ofstream outfile("ped_test.dat");
//   Double_t maxbin, pedx, maxx, minx, mean = 0, sigma=0; 
//   Int_t thres = 0;
//   TF1* gausfit=new TF1("gausfit","gaus",0,2000);
//   for (Int_t i = 0; i < NUMADCSLOTS;i++) {
//     outfile << "slot=" << slotindADC[i] << endl;
//     for (Int_t j = 0; j < NUMCHANA;j++) {
//       maxbin =  hadcped[i][j]->GetMaximumBin();
//       pedx =  hadcped[i][j]->GetBinCenter(maxbin);
//       maxx = pedx + 50;
//       minx = pedx - 50; 
//       hadcped[i][j]->Fit("gausfit","","",minx,maxx);
//       mean  = gausfit->GetParameter(1);
//       sigma = gausfit->GetParameter(2);
//       thres = int(mean + 5*sigma);
//       outfile <<  thres << endl;

//       cout<<"slot "<<slotindADC[i]<<" index "<< j <<" mean,sigma,thres "<< mean <<" "<< sigma <<" "<< thres <<endl;

//     }
//   }
//   outfile.close();
// }

// read out pedestal channels from file
void readpedsup(){
  ifstream infile("ped_test.dat");
  char thres[4];
  Int_t slot = 0;
  char pedend[4];
  char header[8]; // will contain "slot=.." where .. are the slot digits
  char slotstr[2]; // helps parsing header

  for (Int_t i = 0; i < NUMADCSLOTS;i++) {
    infile >> header; // slot=..
    //cout << "readpedsup(): read line: " << header << endl;
    slotstr[0] = header[5]; // first slot digit
    slotstr[1] = header[6]; // second slot digit
    slot = atoi(slotstr);   // convert header substring (saved in slotstr) to slot#
    cout << "readpedsup(): reading pedestals of ADC slot " <<slot << "..." << endl;
    for (Int_t j = 0; j < NUMCHANA;j++) {
      infile >>  thres >> pedend;
      ped_threshold[i][j] = atoi(thres);
      ped_ends[i][j] = atoi(pedend);
      if (debug) cout << "readpedsup(): chan = "<< j<<" thres = '" << thres << "' value = " << ped_threshold[i][j] <<" pedend = " << ped_ends[i][j] <<endl;
    }
  }

  infile.close();
  cout << "readpedsup(): done! "<<endl;
}


// Int_t mapping_read(const char* mapfilein){

//   // pre-processing, let's remove all non-data lines; the result is
//   // saved to a temp file
//   char file[100];
//   sprintf(file,"%s.tmp", mapfilein);

//   gROOT->ProcessLine(Form(".! awk '/^[ \t]*[0123456789]/' %s > %s",mapfilein,file));
  
//   // now we read the temp file
//   Int_t pmtID, ADCslot, ADCchan, TDCslot, TDCchan, NINOID;

//   ifstream in(file);

//   // resetting maps
//   for (Int_t ii = 0; ii < MAXPMTS; ii++){
//       map_pmt_id[ ii ] = -1;
//       map_nino_id[ ii ] = -1;
//       for (Int_t jj = 0; jj < NPIXELS; jj++){
// 	map_adc_slot[ ii ][ jj ] = -1;
// 	map_tdc_slot[ ii ][ jj ] = -1;
// 	map_adc_chan[ ii ][ jj ] = -1;
// 	map_tdc_chan[ ii ][ jj ] = -1;
//       }
//   }
//   for (Int_t ii = 0; ii < NUMADCSLOTS; ii++){
//     for (Int_t jj = 0; jj < NUMCHANA; jj++){
//       map_adc_pmt[ ii ][ jj ] = -1;
//       map_adc_pixel[ ii ][ jj ] = -1;
//     }
//   }
//   for (Int_t ii = 0; ii < NUMTDCSLOTS; ii++){
//     for (Int_t jj = 0; jj < NUMCHANT; jj++){
//       map_tdc_pmt[ ii ][ jj ] = -1;
//       map_tdc_pixel[ ii ][ jj ] = -1;
//     }
//   }

//   Int_t nlines = 0;
//   while ( !in.eof() ){
//     nlines++; // we start at 1

//     in >> pmtID;
//     in >> ADCslot;
//     in >> ADCchan;
//     in >> TDCslot;
//     in >> TDCchan;
//     in >> NINOID;

//     //cout<< "read: " << pmtID <<" "<< ADCslot <<" "<<  ADCchan <<" "<<  TDCslot <<" "<<  TDCchan <<" "<<  NINOID <<endl;

//     map_nino_id[ nlines ] = NINOID;
//     map_pmt_id[ pmtID ] = nlines;
 
//     // pmt, pixel -> slot, chan
//     //-------------------------
//     // unroll channel mappings
//     for (Int_t i = 1; i < NPIXELS; i++){
//       // slots do not change with pixel
//       map_adc_slot[ nlines ][ i ] = ADCslot;
//       map_tdc_slot[ nlines ][ i ] = TDCslot;
//       // channels do
//       map_adc_chan[ nlines ][ i ] = ADCchan + i;
//       map_tdc_chan[ nlines ][ i ] = TDCchan + i;
//     }

//     // slot, chan -> pmt, pixel
//     //-------------------------
//     for (Int_t pixel = 1; pixel < NPIXELS; pixel++){
//       map_adc_pmt  [ ADCslot ][ ADCchan + pixel - 1] = pmtID;
//       map_adc_pixel[ ADCslot ][ ADCchan + pixel - 1] = pixel;  

//       map_tdc_pmt  [ TDCslot ][ TDCchan + pixel - 1] = pmtID;
//       //for (Int_t j=0; j<6; j++) map_tdc_pixel[ TDCslot ][ TDCchan + pixel + j ]= pixel;
//       map_tdc_pixel[ TDCslot ][ TDCchan + pixel - 1]= pixel;
//     }
//   }
  
//   in.close();
 
//   if (nlines > MAXPMTS) cout << "WARNING: number of pmt config lines are > MAXPMTS = " << MAXPMTS << endl;

//   if (debug){
//     cout << "Test pmt, pixel -> slot, chan" <<endl;
//     for (Int_t i=0; i<nlines+1; i++){
//       for (Int_t j=0; j<16; j++){
// 	cout << "PMT "<< i <<" pixel "<< j 
// 	     <<" adc slot = " << map_adc_slot[ i ][ j ]
// 	     <<" tdc slot = " << map_tdc_slot[ i ][ j ]
// 	     <<" adc chan = " << map_adc_chan[ i ][ j ]
// 	     <<" tdc chan = " << map_tdc_chan[ i ][ j ]
// 	     << endl;
//       }
//     }
    
//     cout << "Test slot, chan -> pmt, pixel" <<endl;
//     cout << "Now ADC:"<<endl;
//     for (Int_t i=0; i<NUMADCSLOTS; i++){
//       for (Int_t j=0; j<NUMCHANA; j++){
// 	cout << "adc slot "<< i <<" chan "<< j 
// 	     <<" pmt = " << map_adc_pmt[ i ][ j ]
// 	     <<" pixel = " << map_adc_pixel[ i ][ j ]
// 	     << endl;
//       }
//     }
//     cout << "Now TDC:"<<endl;
//     for (Int_t i=0; i<NUMTDCSLOTS; i++){
//       for (Int_t j=0; j<NUMCHANT; j++){
// 	cout << "tdc slot "<< i <<" chan "<< j 
// 	     <<" pmt = " << map_tdc_pmt[ i ][ j ]
// 	     <<" pixel = " << map_tdc_pixel[ i ][ j ]
// 	     << endl;
//       }
//     }
//   }
  
//   map_len = nlines;
//   return nlines;
// }

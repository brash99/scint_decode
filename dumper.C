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

#define NUMADCSLOTS 4   // ADC # of slots. =4 for test crate, =4 ? for the other ones
#define NUMTDCSLOTS 3   // TDC # of slots. =3 for test crate, =9 for the other ones

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "THaCodaFile.h"
#include "THaEtClient.h"


using namespace std;

void dump(int *data);
void usage();
void clear();
void readpedsup();

// Global data 
// -----------
Int_t fWordSeen;
Int_t evlen, evtype, evnum;
Int_t *roclen = new Int_t[MAXROC];
Int_t myroc = 14;
Int_t runno = -1;

//ADC
//Int_t numslothitsADC[NUMADCSLOTS];
Int_t adcdat[NUMADCSLOTS][NUMCHANA];
Int_t adcchan[NUMADCSLOTS][NUMCHANA];
Int_t slotindADC[NUMADCSLOTS];

Int_t pedestal_threshold[NUMADCSLOTS][NUMCHANA];
char dumpfile[20];
ofstream fdumpfile;


//TDC
//Int_t numslothitsTDC[NUMTDCSLOTS];
Int_t tdctdat[NUMTDCSLOTS][NUMCHANT];
Int_t tdcldat[NUMTDCSLOTS][NUMCHANT];
Int_t tdcchan[NUMTDCSLOTS][NUMCHANT];
Int_t slotindTDC[NUMTDCSLOTS]={14,15,16};

int main(int argc, char* argv[])
{
  //============ Initializations ============
  THaCodaData *coda;      
  Int_t istatus;
  Int_t lslot;

  Int_t ievent = 0;
  // Int_t fHasHeader=1;   // =1 for 1877 and 1881, but =0 for 1875

  // Int_t numLocal = 0;
  // Int_t branchnum = 999;

  // Int_t myslot[4], mychan[4];
  char dtitle[100];

  // Int_t *irn = new Int_t[MAXROC];
  // Int_t *rocpos = new Int_t[MAXROC];

  Int_t choice1 = 1; // =1 CODA file, else ET connection

  Int_t maxevent = MAXEVENTS;


  // ----- argc, argv -----
  if (argc >= 1) {
    runno = atoi(argv[1]); // required arg: run#
  }
  if (argc >= 2){
    myroc = atoi(argv[2]); // optional arg: roc#
  }
  if (argc >= 3){ // optional arg: maxevent
    maxevent = atoi(argv[3]);
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

  cout << "my roc is " << myroc << endl;

  for (int i = 0; i < NUMADCSLOTS; i++){
    slotindADC[i] = lslot;
    lslot++;
  }

  sprintf(dtitle,"../data/scint_%d.dat",runno);

  cout << "Fastbus analysis "<<endl;
  cout << "Events to process "<<maxevent<<endl;

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

  sprintf(dumpfile,"sbs_%d_dump.txt", runno);

  fdumpfile.open(dumpfile);
  fdumpfile << "Ev.num A_T slot / chan ped_edge value / ..."<<endl;

  readpedsup();
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
      dump( coda->getEvBuffer() );
    }
  }

 end1:

  fdumpfile.close();
  cout << endl;

  printf("\n");
  printf("Number of events analyzed = %d \n",ievent);

  coda->codaClose();

  return 0;
}; //end of main function

void clear() {
  fWordSeen = 0;
  //ADC
  memset (adcdat, 0, NUMADCSLOTS*NUMCHANA*sizeof(Int_t));
  memset (adcchan, 0, NUMADCSLOTS*NUMCHANA*sizeof(Int_t));
  //  memset (numslothitsADC, 0, NUMADCSLOTS*sizeof(Int_t));

  //TDC
  //  memset (numslothitsTDC, 0, NUMTDCSLOTS*sizeof(Int_t));
  memset (tdctdat, 0, NUMTDCSLOTS*NUMCHANT*sizeof(Int_t));
  memset (tdcldat, 0, NUMTDCSLOTS*NUMCHANT*sizeof(Int_t));
  memset (tdcchan, 0, NUMTDCSLOTS*NUMCHANT*sizeof(Int_t));
}

void dump (int* data) {
  // ----- init -----
  Int_t ichan = 0, rdata = 0;
  evlen = data[0] + 1;
  evtype = data[1]>>16;
  evnum = data[4];

  Int_t edgetype = 0;
  int slotnew=0;

  char buffer[30];

  if (evtype > 10) return;
  //if (evnum < 10 ) return; 
// Sometimes, the first few events are bad. No idea why. Let's skip them

  // dump event?
  
  if (0){
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
  cout << " START LOOKING AT  at index = " << index << " last = "<< indexlast << " " << roclen[myroc]<< endl;
  while (((data[index]&0xffffffff)!=0x0da000011)) {
    index++;
  }
  while (((data[indexlast]&0xffffffff)!=0x0da000022)) {
    indexlast++;
  }
  cout << " ADC header at index = " << index << " last = "<< indexlast << endl;

  // ----- loop to find header of TDC block read. Start at indexlast, end at indexend -----
  int indexend=indexlast;
  while (((data[indexend]&0xffffffff)!=0x0da000033)) {
    indexend++;
  }
  cout << " TDC header at index = " << indexlast << " last = "<< indexend << endl;

  cout << endl;
  cout << "Ev.num A_T slot / chan ped_edge value / ..."<<endl;


  // ----- main ADC loop -----
  slotnew = 0;
  for (int j = index+1; j < indexlast; j++) {
    //printf("data[%d] = 0x%x = (dec)  \n",j,data[j]);

    int slot = (data[j]&0xf8000000)>>27; 
    int slot_ndata = 0;
    Int_t slotindex=0;

    if (slot!=slotnew) {
      slotindex=0;
      slotnew=slot;
      slot_ndata= (data[j]&0x7f) - 1;
      //cout << slot_ndata << " " << slot << endl;
      for (int jj=0;jj<NUMADCSLOTS;jj++) {
        if (slot==slotindADC[jj]) slotindex=jj;
      }
      if (slot_ndata > NUMCHANA) {
        printf("*** ADC bad evnum= %d,slot_ndata = %d, data[%d] = 0x%x = (dec)  \n",evnum,slot_ndata,j,data[j]);
        return;
      }
    }
    //cout << "ADC slot "<<slot<< " slot index = " << slotindex << " slot indata = " << slot_ndata << endl;


    sprintf(buffer,"%6d adc %d", evnum, slotindex);
    cout << buffer;
    fdumpfile << buffer;


    //
    for (int ii=0;ii<slot_ndata;ii++){
      j++;
      //printf("data[%d] = 0x%x = (dec) %d %d  \n",j,data[j],(data[j]&0xfe0000)>>17,data[j]&0x3fff);
      ichan = (data[j]&0xfe0000)>>17; // 1881
      rdata = data[j]&0x3fff;  // 1881
      //	   ichan = (data[j]&0xfe0000)>>17; // 1877
      //rdata = data[j]&0xffff;  // 1877
      if (ichan >= 0 && ichan < NUMCHANA && rdata >= 0) {
        adcdat[slotindex][ii] = rdata; 
        adcchan[slotindex][ii] = ichan; 

	sprintf(buffer," / %2d %4d %4d", ichan, pedestal_threshold[slotindex][ichan], rdata);
	cout << buffer;
	fdumpfile << buffer;

      }
    }
    //numslothitsADC[slotindex]=slot_ndata;

    cout<< endl;
    fdumpfile << endl;

  } // end of ADC loop

  cout<< endl;
  fdumpfile << endl;

  // ----- main TDC loop -----

  //adapted from ADC loop
  slotnew=0;

  for (int j = indexlast+1; j < indexend; j++) {
    //printf("data[%d] = 0x%x = (dec)  \n",j,data[j]);

    int slot = (data[j]&0xf8000000)>>27; 
    int slot_ndata = 0;
    Int_t slotindex=0;

    if (slot!=slotnew) {
      slotindex=0;
      slotnew=slot;
      //      slot_ndata= (data[j]&0x7f) - 1;
      slot_ndata= (data[j]&0x7ff) - 1;
      //cout << slot_ndata << " " << slot << endl;
      for (int jj=0;jj<NUMTDCSLOTS;jj++) {
        if (slot==slotindTDC[jj]) slotindex=jj;
      }
      if (slot_ndata > NUMCHANT) {
        printf("*** TDC bad evnum= %d,slot_ndata = %d, data[%d] = 0x%x = (dec)  \n",evnum,slot_ndata,j,data[j]);
        return;
      }
    }
    //cout << "TDC slot "<<slot<< " slot index = " << slotindex << " slot indata = " << slot_ndata << endl;

    bool has_leading = false;
    bool has_trailing = false;
    //bool more_than_1 = false;

    //numslothitsTDC[slotindex]=slot_ndata;

    sprintf(buffer, "%6d tdc %d", evnum, slotindex);
    cout << buffer;
    fdumpfile << buffer;

    if (slot_ndata > 0){ // we have TDC data!
      for (int ii=0;ii<slot_ndata;ii++){
	j++;
	//printf("data[%d] = 0x%x = (dec) %d %d  \n",j,data[j],(data[j]&0xfe0000)>>17,data[j]&0xffff);
      ichan = (data[j]&0xfe0000)>>17; // 1877
      rdata = data[j]&0xffff;  // 1877
      if (ichan >= 0 && ichan < NUMCHANT && rdata >= 0) {
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
  	    //more_than_1 = true;
  	  }
  	}
  	if (edgetype == TRAILING_EDGE) {
  	  if (!has_trailing) {
  	    // first leading edge encountered for this event
  	    has_trailing = true;
  	  } else {
  	    // what?!? more than one trailing edge?
  	    //more_than_1 = true;
  	  }
  	}

	if (edgetype == LEADING_EDGE) {
	  sprintf(buffer, " / %2d LEAD %4d", ichan, rdata);
	} else {
	  sprintf(buffer, " / %2d TRAIL %4d", ichan, rdata);
	}

	cout << buffer;
	fdumpfile << buffer;


	//is_good = (is_good) && (has_leading) && (has_trailing);// && (!more_than_1);
	
	// tdcdat[edgetype][slotindex][ii] = rdata; 
      } // matches if (ichan >= 0 && ichan < NUMCHANT && rdata >= 0) {


      } // matches for (int ii=0;ii<slot_ndata;ii++){

      cout<< endl;
      fdumpfile << endl;


    }// matches if (slot_ndata > 0){


    cout<< endl;
    fdumpfile << endl;

  } // end of TDC loop, matches for (int j = indexlast+1; j < indexend; j++) 


} // end of dump

void readpedsup(){
  ifstream infile("ped_test.dat");
  char thres[4];//, endp[4];
  Int_t slot = 0;
  char header[8]; // will contain "slot=.." where .. are the slot digits
  char slotstr[2]; // helps parsing header

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
      //if (debug) cout << "readpedsup(): chan = "<< j<<" thres = '" << thres << "' value = " << pedestal_threshold[i][j] <<endl;

      // if ( debug && (pedestal_threshold[i][j] < 1 || pedestal_threshold[i][j] > MAXPED) ) cout << "WARNING: pedestal of physical slot "<< slot <<" channel "<< j <<" is "<< pedestal_threshold[i][j] <<" which is out of boundary!"<<endl;

    }
  }
  //debug=0;
  infile.close();
  cout << "readpedsup(): done! "<<endl;
}

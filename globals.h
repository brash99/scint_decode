// DEFINES

#define MAXROC     32
#define MAXHITS    100
#define MAXEVENTS  1000000
// constants used to distinguish which case the bit 16 of TDC data word is referring to:
#define TRAILING_EDGE 0
#define LEADING_EDGE  1

#define NUMCHANA   64 // ADC chans per slot
#define NUMCHANT   96*2*3 // (TDC chans per slot) * (leading+trailing edge) * (number of pulses recorded)

// If the following values change, (un)comment as well the relevant initial values of
//     nbad_adc[NUMADCSLOTS], nbad_tdc[NUMTDCSLOTS], slotindTDC[NUMTDCSLOTS]
// which are below
#define NUMADCSLOTS 4   // ADC # of slots. =4 for test crate, =4 ? for the other ones
#define NUMTDCSLOTS 3   // TDC # of slots. =3 for test crate, =9 for the other ones

// INCLUDES
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

// GLOBAL VARS
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
Int_t a[15][17];
Int_t t[15][17];

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

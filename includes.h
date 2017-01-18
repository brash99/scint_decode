#include <stdio.h>

#include <fstream>

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

// This macro, to be called as
// .x delta1.C( PARAMS )
// where the order of the PARAMS is defined below, takes as arguments
// a NINO ID (range 1-14, correspondent to known ADC slot-chans) and
// the "pedestal" and "signal from pulser" runs.

// The parameter ADCstrip allow user to specify to which ADC strip
// (the black cables labeled "card1" up to "card14") the NINO is
// connected to.

// It outputs a file, named after the NINO ID, containing an header
// and a row with the deltas "peak position in signal run" - "peak
// position in pedestal run" for each NINO chan

static const int NUMADCSLOTS = 4;
static const int NUMTDCSLOTS = 3;
static const int NUMCHANA = 64;
static const int NUMCHANT = 96;

int run; // run number. Used in the titles
float positions[2][16], rms[2][16]; // NINO has 16 chans

#include <fstream.h>
#include <TCanvas.h>

void init(int runno);
void extract(int slot, int StartChan, int shift, int type);
TString LSdelay(int slot, int chan);

void delta1(int NINOID=1, int ADCstrip, int pedestalRun1_8, int signalRun1_8, int pedestalRun9_16, int signalRun9_16){
  //TCanvas *c1 = new TCanvas();
  TString filename;
  TFile *file0;
  TTree *t;
  
  int deltac;
  int slot;

  TString delay;

  if ( ADCstrip >= 1 && ADCstrip <= 4 ){
    slot = 0;
    deltac = 16 * (ADCstrip - 1);
  };
  if (ADCstrip >=5 && ADCstrip <= 8){
    slot = 1;
    deltac = 16 * (ADCstrip - 5);
  };
  if (ADCstrip >=9 && ADCstrip <= 12){
    slot = 2;
    deltac = 16 * (ADCstrip - 9);
  };
  if (ADCstrip >=13 && ADCstrip <= 14){
    slot = 3;
    deltac = 16 * (ADCstrip - 13);
  };

  delay = LSdelay(slot, deltac);

  run=pedestalRun1_8;
  filename.Form("sbs_%d_14.root",run);
  file0 = TFile::Open(filename);
  t = (TTree *)file0->Get("t"); 
  extract(slot, deltac, 0, 0);
  file0->Close();

  run=pedestalRun9_16;
  filename.Form("sbs_%d_14.root",run);
  file0 = TFile::Open(filename);
  t = (TTree *)file0->Get("t"); 
  extract(slot, deltac, 8, 0);
  file0->Close();

  run=signalRun1_8;
  filename.Form("sbs_%d_14.root",run);
  file0 = TFile::Open(filename);
  t = (TTree *)file0->Get("t");
  extract(slot, deltac, 0, 1); 
  file0->Close();

  run=signalRun9_16;
  filename.Form("sbs_%d_14.root",run);
  file0 = TFile::Open(filename);
  t = (TTree *)file0->Get("t");
  extract(slot, deltac, 8, 1); 
  file0->Close();

  char outfile[100];
  if (NINOID < 10) {
    sprintf(outfile,"delta_NINO_0%d.txt", NINOID);
  } else {
    sprintf(outfile,"delta_NINO_%d.txt", NINOID);
  }
  
  ofstream file(outfile);
  cout<<"Opened output file "<< outfile <<endl;

  TString format;
  float pos1, pos2, rms1, rms2;
  bool error = false;
  int delta;

  // header
  format = "#NINOID d1   d2   d3   d4   d5   d6   d7   d8   d9   d10  d11  d12  d13  d14  d15  d16 ADCstrip delay mean";

  if (NINOID == 1) file << format <<endl;

  cout << "NINO ID " << NINOID <<" Cable "<< ADCstrip << endl;
  if (NINOID < 10){
    file << NINOID << "     ";
  } else {
    file << NINOID << "    ";
  }

  float mean = 0.;
  for (int chan = 0; chan < 16; chan++){
    pos1 = positions[0][chan];
    pos2 = positions[1][chan];
    //rms1 = rms[0][chan];
    //rms2 = rms[1][chan];
    delta = int(pos2 - pos1 + 0.5); // approx to nearest integer

    //    format.Form(" %.2f", delta);
    format.Form("  %d", delta);

    file << format;
    cout << "chan: "<< chan << " ped: "<< pos1 << " sig: "<< pos2 << " delta: "<< delta;
    if ( delta < 50 ) {
      cout << " ***** WARNING DELTA TOO LOW ****";
      file << "<*";
      error = true;
    };
    cout << endl;

    mean += delta;
  }

  file <<"    "<< ADCstrip <<"     "<< delay <<"   "<< mean/16;

  if (error) file << " !CHANNEL(S) TO CHECK!";
  
  file << endl;
  file.close();

  delete c1;
  cout<<"Done"<<endl;
}

void init(int runno){
	TString filename;
	filename.Form("sbs_%d_14.root",runno);
	TFile *_file0 = TFile::Open(filename);

	run=runno;

	// read tree
	TTree *t = (TTree *)_file0->Get("t"); 
}

// slot = ADC slot
// StartChan = first ADC chan of whole NINO (0, 16, 32, 48);
// shift = 0 or 8, depending on the group of 8 NINO chans (1-8 or 9-16) we want to look at;
// type = 0 for "pedestal" run, 1 for "input signal from pulser" run
void extract(int slot, int StartChan, int shift, int type){
  int slot;
  TString draw;
  int nbin=4100; // 100
  int min=-100, max=4000; // max 500 is too low
  TH1D *htmp=new TH1D("htmp","htmp",nbin,min,max);

   for (int chan = StartChan + shift; chan < StartChan + shift + 8; chan++){
    draw.Form("adc[%d][%d]>>htmp",slot,chan);
    t->Draw(draw);
    positions[ type ][ chan - StartChan ] = htmp->GetMean();
    rms [ type ][ chan - StartChan ] = htmp->GetRMS();

    // cout << run <<" "<< slot <<" "<< chan <<" "<< type <<" "<< htmp->GetMean() << " "<< htmp->GetRMS();
    // if ( htmp->GetMean() < 100 ) cout << " ***** WARNING MEAN TOO LOW ****";

    cout << endl;
  }
}

TString LSdelay_aux(int slot, int chan){
  TString delay = "L";

  switch (slot){
  case 0:
    if (chan <= 15) delay = "U"; // U for Unkwnown, not yet measured
    if (chan >= 48 && chan <= 55) delay = "S";
    break;
  case 1:
    if (chan >= 0 && chan <= 7) delay = "S";
    break;
  case 2:
    if (chan >= 8 && chan <= 23) delay = "S";
    if (chan >= 32 && chan <= 39) delay = "S";
    break;
  case 3:
    if (chan >= 16 && chan <= 23) delay = "S";
    break;
  default:
    delay = "U"; // wrong slot number
  }

  return delay;
}

TString LSdelay(int slot, int chan){
  return LSdelay_aux(slot, chan) + " " + LSdelay_aux(slot, chan+8);
}

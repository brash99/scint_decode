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

Int_t baread[NBARS]; // ADC of bar as read from database
Int_t ba[NBARS], bt[NBARS]; // ADC of bar corrected for gain, TDC of bar
Double_t gainfactor[NBARS];
Int_t aslotbar[NBARS], achanbar[NBARS], tslotbar[NBARS], tchanbar[NBARS];

// read-from-database mapping (work in progress):
TString dbname = "database.txt";
void readdb(TString fname = dbname){
  ifstream in;
  in.open(fname);
  TString header;
  in >> header;
  int i;
  int bar, placeholder;

  for (i=0; i<NBARS; i++){
    in >> bar 
       >> placeholder >> placeholder >> placeholder 
       >> aslotbar[ bar] >> achanbar[ bar ] 
       >> tslotbar[ bar ] >> tchanbar [bar ]
       >> placeholder >> placeholder
       >> gainfactor[ bar ];
  }
  in.close();
}

// to do: create a tree friend of t, then try and fill it with ba[],
// baread[], bt[], calculated from the arrays read from database via
// TString.Form(); this should trick root, whose TFormula class, used
// in t->Draw(), cannot call a function which return a string, it
// accept just integers, (thanks, root dev team, mapping would be too
// simple otherwise)


// manual mapping (USE THIS!):

// pixel (counting from 0) -> ADC slot & chan
Int_t aslot(Int_t pixel){
  Int_t slot = -1;
  if ( (pixel >= 0 && pixel <= 15) || 
       (pixel >= 32 && pixel <= 47) || 
       (pixel >= 64 && pixel <= 79) ||
       (pixel >= 96 && pixel <= 111) ) slot = 0;
  if ( (pixel >= 16 && pixel <= 31) ||
       (pixel >= 128 && pixel <= 143) ) slot = 1;
  if ( (pixel >= 48 && pixel <= 63) ||
       (pixel >= 80 && pixel <= 95) ||
       (pixel >= 112 && pixel <= 127) ||
       (pixel >= 144 && pixel <= 159) ) slot = 2;
  if ( (pixel >= 160 && pixel <= 223) ) slot = 3;
  return slot;
}

Int_t achan(Int_t pixel){
  Int_t chan = -1;
  if (pixel >= 0 && pixel <= 15) chan = pixel;
  if (pixel >= 16 && pixel <= 31) chan = pixel + 32;
  if (pixel >= 32 && pixel <= 47) chan = pixel - 16;
  if (pixel >= 48 && pixel <= 63) chan = pixel - 48;
  if (pixel >= 64 && pixel <= 79) chan = pixel - 32;
  if (pixel >= 80 && pixel <= 95) chan = pixel - 64;
  if (pixel >= 96 && pixel <= 111) chan = pixel - 48;
  if (pixel >= 112 && pixel <= 127) chan = pixel - 80;
  if (pixel >= 128 && pixel <= 143) chan = pixel - 128;
  if (pixel >= 144 && pixel <= 159) chan = pixel - 96;
  if (pixel >= 160 && pixel <= 175) chan = pixel - 128;
  if (pixel >= 176 && pixel <= 191) chan = pixel - 176;
  if (pixel >= 192 && pixel <= 207) chan = pixel - 144;
  if (pixel >= 208 && pixel <= 223) chan = pixel - 192;
  return chan;
}

// global pixel number 0..NBAR
Int_t globalpixel(Int_t pmt, Int_t pixel){
  return (pmt-1)*NPIXEL + (pixel-1);
}

// adc slot of NINO ID
int aslotNINO(int NINOID){
  int value=-1;
  switch (NINOID){
  case 2:
  case 4:
  case 5: 
  case 18:
    value=0; break;
  case 3:
  case 7:
    value=1; break;
  case 8:
  case 11:
  case 17:
  case 19:
    value=2; break;
  case 6:
  case 12:
  case 13:
  case 14:
    value=3; break;
  default:
    cout <<" Unknown NINO ID"<<endl;
  }
  return value;
} 
// first ADC chan of NINO ID (the 15 channels which follow are
// supposed to be in the same order in both NINO card and ADC module)
int achanNINO(int NINOID){
  int value=-1;
  switch(NINOID){
  case 7:
  case 12:
  case 17:
  case 18:
    value=0; break;
  case 2:
  case 14:
  case 19:
    value=16; break;
  case 4:
  case 8:
  case 13:
    value=32; break;
  case 3:
  case 5:
  case 6:
  case 11:
    value=48; break;
  default:
    cout<<" Unknown NINO ID"<<endl;
  }
  return value; 
}

// used by Fastbus_main1.C:
//-------------------------
// pmt (>=1, <=14) -> ADC slot
Int_t handmapping_adc_slot(Int_t pmt){
  switch (pmt){
  case 1:
  case 3:
  case 5:
  case 7:
    return 0;
    break;
  case 2:
  case 9:
    return 1;
    break;
  case 4:
  case 6:
  case 8:
  case 10:
    return 2;
    break;
  case 11:
  case 12:
  case 13:
  case 14:
    return 3;
    break;
  default:
    return -1;
  }
}

// pmt (>=1,<=14), pixel (>=1, <=16) -> ADC channel
Int_t handmapping_adc_chan(Int_t pmt, Int_t pixel){
  pixel--;
  switch(pmt){
  case 1:
  case 4:
  case 9:
  case 12:
    return pixel;
    break;
  case 3:
  case 6:
  case 14:
    return 16+pixel;
    break;
  case 5:
  case 8:
  case 11:
    return 32+pixel;
    break;
  case 2:
  case 7:
  case 10:
  case 13:
    return 48+pixel;
    break;
  default:
    return -1;
  }
}

// pmt (>=1, <=14) -> TDC slot
Int_t handmapping_tdc_slot(Int_t pmt){
  switch (pmt){
  case 1:
  case 3:
  case 5:
  case 7:
  case 9:
  case 11:
    return 0;
    break;
  case 2:
  case 4:
  case 6:
  case 8:
  case 10:
  case 13:
    return 1;
    break;
  case 12:
  case 14:
    return 2;
    break;
  default:
    return -1;
  }
}

// pmt (>=1,<=14), pixel (>=1, <=16) -> TDC channel
Int_t handmapping_tdc_chan(Int_t pmt, Int_t pixel){
  pixel--;
  switch(pmt){
  case 1:
  case 12:
  case 13:
    return pixel;
    break;
  case 2:
  case 3:
  case 14:
    return 16+pixel;
    break;
  case 4:
  case 5:
    return 32+pixel;
    break;
  case 6:
  case 7:
    return 48+pixel;
    break;
  case 8:
  case 9:
    return 64+pixel;
    break;
  case 10:
  case 11:
    return 80+pixel;
  default:
    return -1;
  }
}
//-------------------------

// static const int NUMADCSLOTS = 4;
// static const int NUMTDCSLOTS = 3;
// static const int NUMCHANA = 64;
// static const int NUMCHANT = 96;
// static const float adc_charge = 50*1e-15; // 50 fC, corresponding to an ADC chan
// static const float e = 1.6e-19; // C, electron charge
// static const int NPMT = 14;
// static const int NPIXEL = 16;
// static const int NBARS = NPMT*NPIXEL;

// WARNING: global pixel starts from 0, but it is displayed as if counting from 1
TCanvas * checkpixels(Int_t threshold = 40){
  Int_t slot, chan, bar;
  TH1D htmp ("htmp","htmp",2200,-200,2000);
  TString buffer, cut, drawme;

  Double_t counts, x[NBARS], y[NBARS];
  Int_t subtot[NPMT], notconn[NPMT][NPIXEL];

  TCanvas * ccheckpixels = new TCanvas("ccheckpixels", "ccheckpixels");

  memset(x, 0, NBARS * sizeof(Int_t));
  memset(y, 0, NBARS * sizeof(Int_t));

  Int_t total = 0;
  Int_t subtotal = 0;
  Int_t subpixel;
  Int_t pmt;
  Int_t totalev;

  for (bar = 0; bar < NBARS; bar++){
    slot = aslot(bar);
    chan = achan(bar);
    drawme.Form("adc[%d][%d]>>htmp", slot, chan);
    cut.Form("adc[%d][%d] > %d", slot, chan, threshold);

    t->Draw(drawme);
    totalev = htmp.GetEntries();

    // apply cuts, store resulting histogram in htmp
    t->Draw(drawme, cut);

    counts = htmp.GetEntries();
    float rms = htmp.GetRMS();

    subpixel = bar % NPIXEL;

    pmt = int(bar/NPIXEL);

    if ( subpixel == 0){    
      subtotal = 0;
    }

    cout << "'pixel' "<< bar+1;
    cout << " (pmt "<< pmt+1 <<" pixel "<< subpixel+1 <<" )";
    cout << " total events: "<< totalev;
    cout << " RMS: "<< rms ;
    cout <<" has "<< counts <<" entries > bin "<< threshold;

    if ( rms < 8 ) {
      cout << " CHECK CONNECTION";
      notconn[pmt][subtotal] = subpixel+1;
      subtotal++;
      total++;
    } else {
      if ((1.*counts)/totalev < 0.0005) {
	cout << " CHECK SIGNAL";
	notconn[pmt][subtotal] = subpixel+1;
	subtotal++;
	total++;
      } else {
	cout << " SEEMS OK";
      }
    }


    // if ((1.*counts)/totalev < 0.0005) {
    //   if (rms > 8) {
    // 	cout << " NO BAR";
    //   } else {
    // 	cout << " NO CONNECTION";
    //   }
    //   notconn[pmt][subtotal] = subpixel+1;
    //   subtotal++;
    //   total++;
    // }

    cout << endl;

    if ( subpixel == NPIXEL - 1){
      cout << "-----------------" << endl;
      subtot[pmt] = subtotal;
    }

    // save stuff in arrays
    x[bar] = bar;
    y[bar] = counts;
  }


  for (int pmt = 0; pmt < NPMT; pmt++){
    cout << "PMT "<< pmt+1 << " has " << subtot[pmt] << " not good pixels";
    if (subtot[pmt] > 0) {
      cout << ":";
      for(int i=0; i<subtot[pmt]; i++){
	cout << " " << notconn[pmt][i];
      }
    }
    cout << endl;
  }
  
  cout << "Total = "<< total <<" for threshold = " << threshold <<endl;

  // time to plot
  buffer.Form("Check of good pixels vs pixel number, run %d; pixel; counts > threshold",run);
  TGraph * gcheckpixels = new TGraph( NBARS, x, y);
  gcheckpixels->SetTitle(buffer);
  gcheckpixels->SetMarkerStyle(21);
  gcheckpixels->SetMarkerSize(0.5);

  gcheckpixels->Draw("AP");
  return ccheckpixels;
}


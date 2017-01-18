// static const int NUMADCSLOTS = 4;
// static const int NUMTDCSLOTS = 3;
// static const int NUMCHANA = 64;
// static const int NUMCHANT = 96;
// static const float adc_charge = 50*1e-15; // 50 fC, corresponding to an ADC chan
// static const float e = 1.6e-19; // C, electron charge
// static const int NPMT = 14;
// static const int NPIXEL = 16;
// static const int NBARS = NPMT*NPIXEL;

// static const int NNEIGH = 4; // number of neighbour bars to consider in cuts.
// to set unconnected bar, see init()

Double_t fitf(Double_t *x,Double_t *par) {
  Double_t arg = 0;
  Double_t fitval = 0;  
  if (par[2]!=0) {
    arg = (x[0] - par[1])/par[2];
  }
  fitval = par[0]*exp(-0.5*arg*arg);
  
  return fitval;
}


TCanvas * Avert(Int_t thresh_current_bar, Int_t thres_neigh_bars, Int_t mean_thresh = 60){
  Int_t slot, chan, bar, n;
  TH1D htmp ("htmp","htmp",2200,-200,2000);
  TString buffer, cut, drawme;

  Double_t average, x[NBARS], y[NBARS], rms;

  TCanvas * cAvert = new TCanvas("cAvert", "cAvert");

  memset(x, 0, NBARS * sizeof(Int_t));
  memset(y, 0, NBARS * sizeof(Int_t));

  cout << "# of neighbours considered: "<<  NNEIGH << endl;

  TF1 *fit = new TF1("fit",fitf,0,2000,3);

  for (bar = 0; bar < NBARS; bar++){
    if ( connected[bar] ){
      slot = aslot(bar);
      chan = achan(bar);
      drawme.Form("adc[%d][%d]>>htmp", slot, chan);
      cut.Form("adc[%d][%d] > %d && adc[%d][%d] < 250", slot, chan, thresh_current_bar, slot, chan); // cut on current

      // let's build the cut on neighbours
      for( int obar = bar - NNEIGH; obar < bar + NNEIGH + 1; obar++){
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
      rms = htmp.GetRMS();
      average = htmp.GetMean();
      n = htmp.GetEntries();

      if (bar % NPIXEL == 0) cout <<"     PMT # " << 1 + bar / NPIXEL << endl;
    
      //     if (rms > 8 && n > 20) { // bar should be connected, and hopefully enough statistics
      fit->SetParameter(1, average);
      fit->SetParLimits(1, average - rms/2, average + rms/2);
      fit->SetParameter(2, rms);
      fit->SetParLimits(2, 0, rms * 2);
      
      htmp.Fit("fit","Q"); // quiet mode. do not print results
      
      average = fit->GetParameter(1);
      rms = fit->GetParameter(2);
      // } else {
      // 	cout << "so not fitted: ";
      // }

    } else { // current bar not connected
      average = 0;
      rms = htmp.GetRMS();
      n = 0;
    }

    printf("global pixel %3d histo mean %6.2f rms %6.2f FIT mean %6.2f rms %6.2f chi2 %6.2f / %d",
	   bar+1, htmp.GetMean(), htmp.GetRMS(), average, rms, fit->GetChisquare(), n-2);

    //    if (rms < 8) cout <<"         CHECK CONNECTION";
    cout << endl;

    // save stuff in arrays
    x[bar] = bar+1;
    //y[bar] = average; // from fit
    y[bar] = htmp.GetMean(); // from histogram
  }

  // time to plot
  buffer.Form("<A-vert> vs pixel number, run %d; pixel; <ADC>",run);
  TGraph * gAvert = new TGraph( NBARS, x, y);
  gAvert->SetTitle(buffer);
  gAvert->SetMarkerStyle(21);
  gAvert->SetMarkerSize(0.5);

  // gAvert->GetXaxis()->SetNoExponent();
  // gAvert->GetXaxis()->SetLimits(1,NBARS);
  // gAvert->GetXaxis()->SetNdivisions(-500 -NPMT);

  gAvert->Draw("AP");
  gPad->SetGrid();

  // for each pmt, calculate average of the mean values
  cout << endl;
  for (int pmt = 0; pmt < NPMT; pmt++){
    double mean = 0.;
    int npoints = 0;
    for (int pixel = 0; pixel < NPIXEL; pixel++){
      bar = pmt * NPIXEL + pixel;
      if (y[bar] >= mean_thresh){
	npoints++;
	mean += y[bar];
      }
    }
    printf("Group %2d mean = %6.2f, computed on %d points", pmt+1, mean/npoints, npoints);
    cout<<endl;
  }
  cout<<endl;
  
  return cAvert;
}


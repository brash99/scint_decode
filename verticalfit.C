void verticalfit(int pmt, int pixel, int threshold, int curthres=90, int neigh=NNEIGH){
  TString drawme, drawme1, cut, tempcut;
  Int_t slot, chan, gpixel; 
  Int_t slotn, chann;

  int lastbin = 300;
  int nbins = 100;

  // An histogram on the stack would be preferable, but then root would not allow the access to the fit function 
  TH1D *htemp = new TH1D("htemp","htemp",nbins, 0, lastbin);
  TH1D *htemp1 = new TH1D("htemp1","htemp1",nbins, 0, lastbin);

  gpixel = globalpixel(pmt, pixel);
  // conversion to adc slot, chan
  slot = aslot(gpixel);
  chan = achan(gpixel);

  drawme.Form("adc[%d][%d]>>htemp", slot, chan);
  drawme1.Form("adc[%d][%d]>>htemp1", slot, chan);

  cut.Form("adc[%d][%d] > %d && adc[%d][%d] < %d", slot, chan, curthres, slot, chan, lastbin);

  for (int curpixel = gpixel - neigh; curpixel < gpixel + neigh+1; curpixel++){
    if (curpixel != gpixel){
      tempcut = "";
      slotn = aslot(curpixel);
      chann = achan(curpixel);
      if (slotn != -1 && chann != -1)
	tempcut.Form(" && adc[%d][%d] < %d", slotn, chann, threshold);
      cut = cut + tempcut;
    }
  }

  t->Draw(drawme1); // whole spectrum
  t->Draw(drawme, cut,"goff");

  //cout << drawme << endl;
  //cout << cut << endl;

  Int_t n = htemp->GetEntries();
  cout << n << " entries after cut" << endl;

  // fit on subrange
  TF1 *fit = new TF1("fit","gaus",curthres, lastbin);
  htemp->Fit(fit,"R","goff");

  gPad->Update();
  TPaveStats *ps = (TPaveStats*)htemp->GetListOfFunctions()->FindObject("stats");
  ps->SetOptStat(10);
  gPad->Modified();
  gPad->Update();

  TString title;
  title.Form("ADC: %d neighbors, with ADC < %d; ADC; counts", neigh, threshold);
  htemp1->SetTitle(title);

  double mean =  fit->GetParameter(1);
  double rms =  fit->GetParameter(2);
  cout << "mean: " << mean << endl;
  cout << "rms: " << rms << endl;
  cout << "num phe: " << (mean/rms)*(mean/rms) << endl;

  htemp->SetLineColor(kBlue);
  htemp1->SetLineColor(kRed);

  htemp1->Draw();
  htemp->Draw("same");
  fit->Draw("same");
}

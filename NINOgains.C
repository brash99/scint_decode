void NINOgains(int NINOID, double normalizeto, int pedestal_run, int group1_run, int group2_run, int threshold=40){
  int i;
  int slot = aslotNINO(NINOID);
  int chan = achanNINO(NINOID);

  int pedposition[NPIXEL]={0};
  int peakposition[NPIXEL]={0};

  TH1D *histo = new TH1D("histo","histo",2100,-100,2000);

  TString drawme;

  double gain;

  // DB#4: NINO ID, NINO channel, gain correction factor
  TString db4name = "db4.TMP";
  ofstream output;
  output.open(db4name);

  // pedestal
  read(pedestal_run);
  cout<<"Reading pedestal run "<< pedestal_run<<endl;
  for (i=0; i<NPIXEL; i++){
    drawme.Form("adcraw[%d][%d]>>histo",slot, i + chan);
    t->Draw(drawme,"","goff");
    pedposition[i] = histo->GetMean();
  }

  // group 1 of NINO channels
  read(group1_run);
  cout<<"Reading group 1 calibration run "<< group1_run<<endl;
  for (i=0; i<NPIXEL/2; i++){
    drawme.Form("adcraw[%d][%d]>>histo",slot, i + chan);
    t->Draw(drawme,"","goff");
    peakposition[i] = histo->GetMean();
  }

  // group 2 of NINO channels
  read(group2_run);
  cout<<"Reading group 2 calibration run "<< group2_run<<endl;
  for (i=NPIXEL/2; i<NPIXEL; i++){
    drawme.Form("adcraw[%d][%d]>>histo",slot, i + chan);
    t->Draw(drawme,"","goff");
    peakposition[i] = histo->GetMean();
  }

  // output
  output << NINOID;
  cout << NINOID <<endl;

  for (i=0; i<NPIXEL; i++){
    gain = 0;
    if (peakposition[i]-pedposition[i] > threshold)
      gain = normalizeto / (peakposition[i]-pedposition[i] + 0.);

    output << " "<< gain;
    cout << i <<" ped = "<< pedposition[i] <<" peak = "<< peakposition[i] <<" gain = "<< gain <<endl;

  }

  output << endl;
  output.close();

  delete histo;
}

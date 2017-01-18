// Bogdan, feel free to add any command you like in this script. To run it in root, just type
// .x bogdan.C
// and press the enter key (in a terminal with a running session of root, of course)

// I expect you to already know what follows, but just in case:
// ============================================================================
// || Do not forget to add a ; at the end of each line, and to save the file ||
// ============================================================================
// if you don't the former, root will give you a warning, and nothing else

// Variable names:
// adc[ slot ][ channel ] = adc with pedestal at 0
// adcraw[ slot ][ channel ] = adc as read
// tdct[ slot ][ channel ] = tdc trailing edge
// tdcl[ slot ][ channel ] = tdc leading edge

void bogdan(){
//sample command
// 200
t->Draw("adc[0][62]","tdct[1][8]>1200");

//t->Draw("tdct[1][8]","adc[3][56]>30&&adc[3][56]<800&&adc[3][55]<30&&adc[3][57]<30&&adc[3][54]<30&&adc[3][58]<30&&adc[3][53]<30&&adc[3][59]<30");


//152
//t->Draw("adc[2][56]","adc[2][56]>30&&adc[2][56]<800&&adc[2][55]<30&&adc[2][57]<30&&adc[2][54]<30&&adc[2][58]<30&&adc[2][53]<30&&adc[2][59]<30");

//150
//t->Draw("adc[2][54]","adc[2][54]>80&&adc[2][54]<800&&adc[2][53]<30&&adc[2][55]<30&&adc[2][52]<30&&adc[2][56]<30&&adc[2][51]<30&&adc[2][57]<30");


}







/*

// for August 24 presentation: run 1042
// begin ////////////////////////////////////////////////////
// selected pmt 3 pixel 5, or adc[0][20], because it is in the middle of a group of 8 channels (so neighbour pixels are easy to identify) with no "bad" pixels, and also because its tdc indexes are [0][20] as well
//----------------------------

// tdc

t->Draw("tdct[0][20]", "tdct[0][20]>-1 && tdct[0][20]<1500");
htemp->SetTitle("TDC; TDC bin; counts");
c1->SetLogy();
c1->SaveAs("run_1042_tdc_0_20.pdf");

//----------------------------

TH1D * htemp1 = new TH1D("htemp1","htemp1",400, -100, 300);
TH1D * htemp2 = new TH1D("htemp2","htemp2",400, -100, 300);
TH1D * htemplow1 = new TH1D("htemplow1","htemplow1", 100, 0, 300);
TH1D * htemplow2 = new TH1D("htemplow2","htemplow2", 100, 0, 300);
TH1D * htemptdc = new TH1D("htemptdc", "htemptdc", 2100, -100, 2000);

// ------------------------------------

// Whole ADC plot

// no cuts
t->Draw("adc[0][20]>>htemp1", "adc[0][20]<300");
htemp1->SetTitle("ADC: no cuts on neighbours; ADC bin; counts");
htemp1->SetLineColor(kRed);
htemp1->Draw();

// now apply tdc cut
t->Draw("adc[0][20]>>htemp2", "adc[0][20]<300 && tdct[0][20]>1200","same");
htemp2->SetLineColor(kBlue);
htemp2->Draw("same");
 c1->SaveAs("run_1042_adc_0_20.pdf");

//low res version

// no cuts
t->Draw("adc[0][20]>>htemplow1", "adc[0][20]<300");
htemplow1->SetTitle("ADC: no cuts on neighbours; ADC bin; counts");
htemplow1->SetLineColor(kRed);
htemplow1->Draw();

// now apply tdc cut
t->Draw("adc[0][20]>>htemplow2", "adc[0][20]<300 && tdct[0][20]>1200","same");
htemplow2->SetLineColor(kBlue);
htemplow2->Draw("same");
c1->SaveAs("run_1042_adc_0_20_low.pdf");

//-------------------------

// horizontal track

t->Draw("adc[0][20]>>htemplow1", "adc[0][20]<300 && adc[0][20]>0 && adc[0][21]>40 && adc[0][21]<80 && adc[0][19]>40 && adc[0][19]<80");
htemplow1->SetTitle("ADC: 1 neighbour per side, with ADC > 40; ADC bin; counts");
htemplow1->SetLineColor(kBlue);
htemplow1->Draw();
c1->SaveAs("run_1042_adc_0_20_1neigh_gt_40_lt_80.pdf");

//---------------------------

// vertical

t->Draw("adc[0][20]>>htemplow1", "adc[0][20]<300 && adc[0][20]>0 && adc[0][21]<30 && adc[0][19]<30 && adc[0][22]<30 && adc[0][18]<30 && adc[0][23]<30 && adc[0][17]<30");
htemplow1->SetLineColor(kRed);
htemplow1->SetTitle("ADC: 3 neighbours per side, with ADC < 30; ADC bin; counts");

// plus cut on current adc:

t->Draw("adc[0][20]>>htemplow2", "adc[0][20]<300 && adc[0][20]>90 && adc[0][21]<30 && adc[0][19]<30 && adc[0][22]<30 && adc[0][18]<30 && adc[0][23]<30 && adc[0][17]<30","same");
htemplow2->SetLineColor(kBlue);
htemplow2->Draw("same");
c1->SaveAs("run_1042_adc_0_20_3neigh_lt_30_currentgt_90.pdf");

//----------------------------

// fit of vertical track

verticalfit(3, 5, 30, 90, 3);
// end ////////////////////////////////////////////////////////

*/

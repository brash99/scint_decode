TFile *f1 = new TFile("sbs_1052_14.root","READ");
TTree *tree1 = (TTree*)f1->Get("T");
TCanvas *c1 = new TCanvas("c1","c1",100,100,500,270);
c1->Divide(2,2)
c1->cd(1)
ADC_chan_slot_20->Draw()
c1->cd(2)
ADC_chan_slot_21->Draw()
c1->cd(3)
ADC_chan_slot_22->Draw()
c1->cd(4)
ADC_chan_slot_23->Draw()

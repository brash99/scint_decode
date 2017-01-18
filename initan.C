
#include <iostream>
#include <fstream>
using namespace std;
// this  script finds the maximum of ADC data and fits a line to an extended region
// after the peak. This is for all ADC from [0][0] to [3][63]. The slopes are then 
// compiled into a single array, and the reciprocal is taken.
// note that there are 4 loops instead of one nested loop. This is because root hates strings
// Parker Reid - Jan 1
void initan() {

 // open file
  int i;
  int n; 
 ofstream myfile;
  myfile.open ("slopean.txt");
    double q[100000000][10000000];
    double result[100000000][100000000];
	// start the loop of ADC channels
	// for i = 0 - 3  {
	// for n = 0 - 63 {
	//   LOGIC HERE 
		//  Get maximum of desired ADC this is 0.
	 // fit a gaussian to the adc
	 // use range 0 to +1 sigma to determine slope
      
    TF1 *myfit = new TF1("myfit", "[0]+[1]*x"); // line fit
    myfit->SetParName(0,"po");
    myfit->SetParName(1,"m"); // this is the slope parameter which we need
    
    int i = 0;
    for ( n = 0; n<64; n++){


    
   TString adcname;
   adcname.Form("adc[0][%i]",n);
   TString adclim;
   adclim.Form("adc[0][%i]>0&&adc[0][%i]<20",n,n);
   
    t->Fit("myfit",adcname, adclim);
    
    q[i][n] = myfit -> GetParameter(1);
    result[i][n] = 1/(q[i][n]);
    myfile <<  i<<" "<< n <<" " << result[i][n] << endl;
}
//loop 2
i++;
 for ( n = 0; n<64; n++){
 

    
   TString adcname;
   adcname.Form("adc[1][%i]",n);
   TString adclim;
   adclim.Form("adc[1][%i]>0&&adc[1][%i]<20",n,n);
   
    t->Fit("myfit",adcname, adclim);
    
    q[i][n] = myfit -> GetParameter(1);
    result[i][n] = 1/(q[i][n]);
    myfile <<  i<<" "<< n <<" " << result[i][n] << endl;
}
// loop 3
i++;
 for ( n = 0; n<64; n++){


    
   TString adcname;
   adcname.Form("adc[2][%i]",n);
   TString adclim;
   adclim.Form("adc[2][%i]>0&&adc[2][%i]<20",n,n);
   
    t->Fit("myfit",adcname, adclim);
    
    q[i][n] = myfit -> GetParameter(1);
    result[i][n] = 1/(q[i][n]);
    myfile <<  i<<" "<< n <<" " << result[i][n] << endl;
} 
// loop 4
i++;
for ( n = 0; n<64; n++){
    

    
   TString adcname;
   adcname.Form("adc[3][%i]",n);
   TString adclim;
   adclim.Form("adc[3][%i]>0&&adc[3][%i]<20",n,n);
   
    t->Fit("myfit",adcname, adclim);
    
    q[i][n] = myfit -> GetParameter(1);
    result[i][n] = 1/(q[i][n]);
    myfile <<  i<<" "<< n <<" " << result[i][n] << endl;
}

myfile.close();
}

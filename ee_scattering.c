#include <cstdio>
#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TFile.h"

void angles()
{
	//Define the max and min values of cosθ
	Int_t xmax = 1;
	Int_t xmin = -1;
	
	// Create a new canvas
	TCanvas *c1 = new TCanvas("c1", " ", 800, 600);
	
	// Create a histogram
	TH1F *h1 = new TH1F("s1", "PDF of polar angle", 100, -1, 1);
	
	// Create the pdf function
	TF1 *f = new TF1("f", "(3./8)*(1 + x*x)", -1, 1);
	
	// Max value of pdf will be for x = -1 or x = 1 since x = cosθ
	Double_t fmax = f->GetMaximum();

	// Create a second histogram
	TH1F *h2 = new TH1F("s2", "PDF of azimuthial angle", 100, -10, 190);	
	
	// Create a 2D histogram
    	TH2F *h2D = new TH2F("h2D", "2D Histogram", 100, -1, 1, 100, -10, 190);

	// Declare the a and b parameters in the pdf which are the max and min azimuthial angles
	Int_t a = 0; //Degrees
	Int_t b = 180; //Degrees
	
	// Create the pdf of the uniform distribution
	TF1 *f2 = new TF1("f2", "[0]", a, b);
	f2->SetParameter(0, 1.0/(b - a));
	
	// Max value of pdf will be for x = -1 or x = 1 since x = cosθ
	Double_t fmax2 = f2->GetMaximum();
	
	Int_t i;
	
	Int_t N = 1000000;
	
	// Generate random numbers from uniform distribution, r1, r2, r3,..
	for(i=0 ; i < N ; i++){
		
		// For the polar angle
		Double_t r_arr = gRandom->Uniform(0, 1);
		Double_t s_arr = gRandom->Uniform(0, 1);
		
		Double_t xi = (xmax - xmin)*r_arr + xmin;
		
		Double_t u = fmax*s_arr;
		
		Double_t fi = f->Eval(xi);
		
		if(fi > u){
		
			h1->Fill(xi);
		
		}
		
		// For the azimuthial angle
		Double_t r_arr2 = gRandom->Uniform(0, 1);
		Double_t s_arr2 = gRandom->Uniform(0, 1);
		
		Double_t xj = (b - a)*r_arr2 + a;
		
		Double_t u2 = fmax2*s_arr2;
		
		Double_t fj = f2->Eval(xj);
		
		if(fj > u2){
		
			h2->Fill(xj);
		
		}
		
		// For the 2D histogram
		if(fi > u && fj > u2){
			
			h2D->Fill(xi, xj);
		}
			
	}
	// Normalize histogram by the total number of events and the bin width
	// v = (counts in each bin)/(total data)*(bin width)
	// This normalization will give the pdf in y-axis
    	Int_t bin; //Declare a variable in order to access the content of each bin
    	
    	Int_t Nbins = h1->GetNbinsX(); //Get number of bins
    	
    	Double_t n = h1->GetEntries(); // Total data in the bins
    	
    	for(bin=0 ; bin <= Nbins ; bin++){
    	
    		// Get counts of each bin
    		Double_t counts = h1->GetBinContent(bin);
    		
    		// Get the bin width
    		Double_t width = h1->GetBinWidth(bin);
    		
    		// Normalize the counts
    		h1->SetBinContent(bin, counts/(n*width));
    	
    	}
    	
    	// Fit the known pdf 
    	h1->Fit(f, "R");
	
	h1->Draw("hist");
	
	h1->GetXaxis()->SetTitle("cos(theta)");
	
	f->Draw("same");
	
	// Normalize histogram by the total number of events and the bin width
	// v = (counts in each bin)/(total data)*(bin width)
	// This normalization will give the pdf in y-axis
    	Int_t Nbins2 = h2->GetNbinsX(); //Get number of bins
    	
    	Double_t n2 = h2->GetEntries(); // Total data in the bins
    	
    	for(bin=0 ; bin <= Nbins2 ; bin++){
    	
    		// Get counts of each bin
    		Double_t counts2 = h2->GetBinContent(bin);
    		
    		// Get the bin width
    		Double_t width2 = h2->GetBinWidth(bin);
    		
    		// Normalize the counts
    		h2->SetBinContent(bin, counts2/(n2*width2));
    	
    	}
    	
    	// Create a second canvas for the azimuthial angle
	TCanvas *c2 = new TCanvas("c2", " ", 800, 600);	
	
    	// Fit the uniform distribution 
    	h2->Fit(f2, "R");
	
	h2->Draw("hist");
	
	h2->GetXaxis()->SetTitle("azimuthial [degrees]");
	
	f2->Draw("same");
	
	// Create a third canvas for the 2D histogram
	TCanvas *c3 = new TCanvas("c3", " ", 800, 600);	

	Double_t n2D = h2D->GetEntries();
	
	Double_t NbinsX = h2D->GetNbinsX(); // Number of bins along x-axis
	Double_t NbinsY = h2D->GetNbinsY(); // Number of bins along y-axis
	
	// Normalize the hist by entering one bin at a time
    	for(i = 0 ; i <= NbinsX ; i++){ //For bins along x-axis
   
		for(Int_t j = 0 ; j <= NbinsY ; j++){// For bins along y-axis
		
			// Get the bin content 
			Double_t counts2D = h2D->GetBinContent(i, j);
			
			// Get the bin width in X and Y
			Double_t widthX = h2D->GetXaxis()->GetBinWidth(i);
			Double_t widthY = h2D->GetYaxis()->GetBinWidth(j);
		
			// Normalize the counts
    			h2D->SetBinContent(i, j, counts2D/(n2D*widthX*widthY));
    			
	    }
    	}
	
	h2D->Draw("colz");
	
	h2D->GetXaxis()->SetTitle("polar");
    	h2D->GetYaxis()->SetTitle("azimuthial");
    	
    	// Save everything in a .root file
    	TFile *myFile = new TFile("Elec_mion.root", "RECREATE");
    	
    	// Save the canvases
  	
  	c1->Write("Polar_canvas");
  	c2->Write("Azimuthial_canvas");
  	c3->Write("2D_canvas");
  	
  	myFile->Close();
	
}

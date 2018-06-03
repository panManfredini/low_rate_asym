double compute_extended_term( double mu, double *fs_fb_array, int n_obs, double Nb){
	
	double ex_term =0.;
	for(int k=0; k<n_obs; k++){
		ex_term += (fs_fb_array[k] -1.) * Nb / ( fs_fb_array[k] * mu + Nb ) ;
	}

	return ex_term;

}

void ex_term_approx_attemp(){

 // ATTEMP to compute the exended term distribution for a fixed mu value
 
 //=================== OPTIONS ============================//
   	int type = 0;
 //=======================================================//

 //==================== PART 1 ===========================//
 
	// likelihood
	pdfLikelihood* likeHood = getTheLikelihoodByType(50, type);
	likeHood->initialize();

	double Nb = plotHelpers::giveNb(likeHood);
	double sig_scale = likeHood->getSignalMultiplier();

	TH2F bkg_pdf = likeHood->getOverallBkgHisto();
	bkg_pdf.Scale(1. / Nb);   // scaling to be integral 1.
	TH2F signal_pdf = likeHood->signal_component->getInterpolatedHisto();
	signal_pdf.Scale( sig_scale );// scaling to be integral 1.
	
	double x_bins[120]; 
	for(int	i = 0; i < 120; i++){
		x_bins[i] = ( TMath::Exp(((double)i)*0.1)  -1. ) /500.;
	}
	TH1F *pdf_of_ratio = new TH1F("pdf_of_ratio" , "", 119, x_bins);

	// loop over histos and fill PDF
	for(int x = 1; x <= bkg_pdf.GetNbinsX(); x++){
		for(int y = 1; y <= bkg_pdf.GetNbinsY(); y++){
			if(bkg_pdf.GetBinContent(x,y) > 1E-7  ) {
				double content = signal_pdf.GetBinContent(x,y) / bkg_pdf.GetBinContent(x,y) ;
				pdf_of_ratio->Fill(content, bkg_pdf.GetBinContent(x,y));
			}
			else { cout << "MISSED!" << bkg_pdf.GetBinContent(x,y) << "    " << signal_pdf.GetBinContent(x,y) << endl; }
		}		
	}

	pdf_of_ratio->Draw();
	cout << "Integral : " << pdf_of_ratio->Integral() << endl;
	pdf_of_ratio->Scale(1./pdf_of_ratio->Integral());  // rescale in order to have integral 1, we miss some entries...

 
 //==================== PART 2 ===========================//
 	double fs_fb[1000] = {0.};
	TRandom3 rambo;
	TH1F *ex_term_distro =  new TH1F("ex_term_distro","real",200, -1000.,1000.);
	TH1F *ex_term_apx1 =  new TH1F("ex_term_apx1","apx 1nd",200, -1000.,1000.);
	TH1F *ex_term_apx2 =  new TH1F("ex_term_apx2","apx 2nd",200, -1000.,1000.);
	TH1F *ex_term_apx3 =  new TH1F("ex_term_apx3","apx 3rd",200, -1000.,1000.);
	TH1F *ex_term_apx4 =  new TH1F("ex_term_apx4","apx 4th",200, -1000.,1000.);
	
	double mu = 3;
	
	//loop over datasets
	for(int i=0; i< 10000; i++){

		double N_obs = 200.;//rambo.Poisson(Nb);	
 		double fs_fb_sum = 0.;
 		double fs_fb_sum_2 = 0.;
 		double fs_fb_sum_3 = 0.;
 		double fs_fb_sum_4 = 0.;
 		double fs_fb_sum_5 = 0.;
	
		// compute N_obs array of fs/fb
		for(int j=0; j < N_obs; j++){
			fs_fb[j] = pdf_of_ratio->GetRandom();
			//cout << "random fs/fb " << fs_fb[j] << endl;
			fs_fb_sum += fs_fb[j];
			fs_fb_sum_2 += fs_fb[j] * fs_fb[j] ;
			fs_fb_sum_3 += fs_fb[j] * fs_fb[j] * fs_fb[j];
			fs_fb_sum_4 += fs_fb[j] * fs_fb[j] * fs_fb[j] *fs_fb[j];
			fs_fb_sum_5 += fs_fb[j] * fs_fb[j] * fs_fb[j] * fs_fb[j] * fs_fb[j];

		}


		TGraph zero_finder(1);


		ex_term_distro->Fill(compute_extended_term( mu, fs_fb, N_obs, Nb));
			
		double ex_term_approx_1  =  (fs_fb_sum - N_obs) - (fs_fb_sum_2 - fs_fb_sum ) * mu / Nb ; 
		double ex_term_approx_2  =  (fs_fb_sum - N_obs) - (fs_fb_sum_2 - fs_fb_sum ) * mu / Nb + ( fs_fb_sum_3 - fs_fb_sum_2 ) * mu*mu / Nb / Nb;
		double ex_term_approx_3  =  ex_term_approx_2 - ( fs_fb_sum_4 - fs_fb_sum_3) *mu*mu*mu / Nb/Nb/Nb;
		double ex_term_approx_4  =  ex_term_approx_3 + ( fs_fb_sum_5 - fs_fb_sum_4) *mu*mu*mu*mu / Nb/Nb/Nb/Nb;
		ex_term_apx1->Fill(ex_term_approx_1);
		ex_term_apx2->Fill(ex_term_approx_2);
		ex_term_apx3->Fill(ex_term_approx_3);
		ex_term_apx4->Fill(ex_term_approx_4);


	}
	new TCanvas();
	ex_term_distro->SetLineWidth(3);
	ex_term_apx1->SetLineWidth(3);
	ex_term_apx2->SetLineWidth(3);
	ex_term_apx3->SetLineWidth(3);
	ex_term_apx4->SetLineWidth(3);

	ex_term_distro->Draw("hist PLC");
//	ex_term_apx1->Draw("hist same PLC");
	ex_term_apx2->Draw("hist same PLC");
	ex_term_apx3->Draw("hist same PLC");
	ex_term_apx4->Draw("hist same PLC");
     gPad->BuildLegend(); 
 //=======================================================//


}

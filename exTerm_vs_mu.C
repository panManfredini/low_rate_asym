double compute_extended_term( double mu, double *fs_fb_array, int n_obs, double Nb){
	
	double ex_term =0.;
	for(int k=0; k<n_obs; k++){
		ex_term += (fs_fb_array[k] -1.) * Nb / ( fs_fb_array[k] * mu + Nb ) ;
	}

	return ex_term;

}

void exTerm_vs_mu(){

 // ATTEMP to compute mu_hat distribution from distro of fs/fb
 // 1) compute distro of fs/fb given bkg only hypo
 // 2) compute distro of sum in "the formula" viam MC (for now)
 // 3) compute mu_hat on many toys (for now) and store in histo
 
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
	cout << "Integral X > 100 = " << pdf_of_ratio->Integral(pdf_of_ratio->FindBin(100), -1) << endl;

 
 //==================== PART 2 ===========================//
 	double fs_fb[1000] = {0.};
	TRandom3 rambo;
	TH1F *mu_distro =  new TH1F("mu_distro","",100, -100.,20.);
	vector <TH1F*> ex_term_distro;
	const unsigned int n_mu_val = 20;
	double mu = 10.;
	double step = 1;
	double mu_values[n_mu_val];

        for(unsigned int k=0 ; k < n_mu_val; k++) { 
		mu_values[k] = mu;
		ex_term_distro.push_back( new TH1F("ex_term_distro"+TString::Itoa(k,10),TString::Format("mu = %1.1f",mu),200, -1000.,1000.) );
		ex_term_distro[k]->SetLineWidth(4);
		mu = mu - step;
	}
	
	//loop over datasets
	for(int i=0; i< 10000; i++){

		int N_obs = 200.; //rambo.Poisson(Nb);	
	
		// compute N_obs array of fs/fb
		for(int j=0; j < N_obs; j++){
			fs_fb[j] = pdf_of_ratio->GetRandom();
		}


		// loop on mu

 		for(unsigned int j=0; j < n_mu_val; j++){
	         	double ex_term = compute_extended_term( mu_values[j], fs_fb, N_obs, Nb);
			ex_term_distro[j]->Fill(ex_term);

		}

	}

	new TCanvas();
	ex_term_distro[0]->Draw("hist PLC");
	ex_term_distro[5]->Draw("hist SAME PLC");
	ex_term_distro[10]->Draw("hist same PLC ");
	ex_term_distro[15]->Draw("hist same PLC");

	   gPad->BuildLegend();	
 //=======================================================//


}


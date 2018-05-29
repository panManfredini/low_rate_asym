double compute_extended_term( double mu, double *fs_fb_array, int n_obs, double Nb){
	
	double ex_term =0.;
	for(int k=0; k<n_obs; k++){
		ex_term += (fs_fb_array[k] -1.) * Nb / ( fs_fb_array[k] * mu + Nb ) ;
	}

	return ex_term;

}

void distro_mu_hat_attemp(){

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

 
 //==================== PART 2 ===========================//
 	double fs_fb[1000] = {0.};
	TRandom3 rambo;
	TH1F *mu_distro =  new TH1F("mu_distro","",100, -100.,20.);
	
	//loop over datasets
	for(int i=0; i< 1000; i++){

		double N_obs = rambo.Poisson(Nb);	
	
		// compute N_obs array of fs/fb
		for(int j=0; j < N_obs; j++){
			fs_fb[j] = pdf_of_ratio->GetRandom();
			//cout << "random fs/fb " << fs_fb[j] << endl;
		}


		TGraph zero_finder(1);

		// loop on mu
		double mu = 10.;
		double step = 0.1;

		bool isPositive = (N_obs -mu -Nb  + compute_extended_term( mu, fs_fb, N_obs, Nb) ) > 0 ;

		for(int j=0; j < 2000; j++){
			
			double ex_term = compute_extended_term( mu, fs_fb, N_obs, Nb);
			double current_f = N_obs -mu -Nb  + ex_term ;
			
			//cout << "mu " << mu << "  ex_term " << ex_term << " current_f " << current_f << endl;

			zero_finder.SetPoint(j, current_f, mu);
			// find when flip sign
			if(isPositive && current_f < 0. ) break;
			else if(!isPositive && current_f > 0.) break;

			mu = mu - step;
		}
	
		double mu_hat = zero_finder.Eval(0.);
		cout << "mu_hat " << mu_hat << endl;
	//	zero_finder.Draw("A*L");
	//	new TCanvas();


		mu_distro->Fill(mu_hat);
	}
	new TCanvas();
	mu_distro->Draw();

       cout << "integral " <<	mu_distro->Integral(mu_distro->FindBin(0.), -1) << endl;
 //=======================================================//


}

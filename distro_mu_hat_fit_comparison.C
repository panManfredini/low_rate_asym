void distro_mu_hat_fit_comparison(){

	// formula analitica per trovare il massimo della likelihood
	// con questo script voglio testare su un mucchio di toys quando la formula fallusce. 
	// N_obs -mu -N_b +  SUM_i_Nobs[ (F_s -F_b)N_B / (mu*F_s +N_b*F_b) ] = 0
	
	// 1) loop over toy datasets and compute minimum with formula
	// 2) Save the computed value in tree of postfit previously done 
	

	
	// general Options:
	int type = 0;
	int Gen = 1;
	const int scan_points = 4000;
	double step = 0.01;
	bool isPositive = false;

	// likelihood
	pdfLikelihood* likeHood = getTheLikelihoodByType(50, type);
	likeHood->initialize();
	likeHood->getPOI()->setCurrentValue(1.);

	
	double Nb = plotHelpers::giveNb(likeHood);
	double sig_scale = likeHood->getSignalMultiplier();

	TH2F bkg_pdf = likeHood->getOverallBkgHisto();
	bkg_pdf.Scale(1. / Nb);   // scaling to be integral 1.
	TH2F signal_pdf = likeHood->signal_component->getInterpolatedHisto();
	signal_pdf.Scale( sig_scale );// scaling to be integral 1.


	// simple data looper using datahandler
	TString toy_name = "low_rate_asy_M50_mu0_G0_V0";
	dataHandler *toyRunner = new dataHandler("toyRunner","../build/RESULTS/GENtrees/low_rate_asy_M50_mu0_G0_V0.root", toy_name + "_0");
	toyRunner->setPrefixTree(toy_name);
        //toyRunner->setPrintLevel(DEBUG);
	
	// Tree of postfits, to be modified adding a branch 'mu_hat_analitic'
	TFile *post_fit_file = TFile::Open("../build/RESULTS/FITtrees/mod_post_fit_low_rate_asy_M50_muTrue0_muFit0_G0.root", "update"); // update a copy of the original file 
	TTree *post_fit_tree = (TTree*)post_fit_file->Get("post_fit_tree");
	float mu_analitic =0 ;
	TBranch *b_mu_analitic = post_fit_tree->Branch("mu_analitic",&mu_analitic,"mu_analitic/F"); 
	
	// file that will contain modified tree
	// doenst work :(
    	//TFile *mod_fit_file = new TFile("../build/RESULTS/FITtrees/mod_post_fit_low_rate_asy_M50_muTrue0_muFit0_G0.root","RECREATE");

 
	float cs1 =0;
	float cs2 =0;

    /// LOOP on many datasets
    for(int k=0; k< 5000; k++){
	
	toyRunner->setTreeIndex(k);
	Long64_t N = toyRunner->getEntries();
	
	post_fit_tree->GetEntry(k); 

	// FIND mu_hat	
	double mu = 10.;
	TGraph zero_finder(1);
	for(int j=0; j<scan_points; j++){


		// COMPUTING EXTENDED TERM 
		double 	extended_term = 0;
		for ( Long64_t i =0; i < N ; i++){

			cs1 = toyRunner->getS1(i);	
			cs2 = toyRunner->getS2(i);	

			// this is fs normalized to 1 
			double fs = signal_pdf.GetBinContent(signal_pdf.FindBin(cs1, cs2));
			// normalized to 1 fb
			double fb = bkg_pdf.GetBinContent( bkg_pdf.FindBin(cs1, cs2));

			extended_term += ( fs - fb ) * Nb / ( mu*fs + Nb*fb);

		}	
		
		double formula =  N -mu -Nb + extended_term ;
		
		if(j == 0) isPositive = ( formula > 0 );
		
		zero_finder.SetPoint(j, formula , mu);

		// find when flip sign
		if(isPositive && formula < 0. ) break;
		else if(!isPositive && formula > 0.) break;

		mu = mu - step;
		

	}
	
	// filling the new branch	
	mu_analitic = zero_finder.Eval(0.);
	if( mu_analitic < -20.) mu_analitic = -20;
	cout << "mu hat " << mu_analitic << endl;
	
	b_mu_analitic->Fill();  


    }

    post_fit_tree->Write();

}








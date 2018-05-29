{

	// formula analitica per trovare il massimo della likelihood
	// con questo script voglio testare se ho fatto la derivata correttamente ;)
	// N_obs -mu -N_b +  SUM_i_Nobs[ (F_s -F_b)N_B / (mu*F_s +N_b*F_b) ] = 0
	
	// APPROXIMATED:: N_obs -mu - N_b -K_b + N_b/mu*K_s = 0

	// 1) plot della formula per alcuni toys, find zeros by eye.
	// 2) find real max with minuit for some toys and compare with formula.
	// 3) find out how far I get using the approximated 
	

	//=================== PART 1 =========================//
	
	// general Options:
	int tree_index = 0;
	int type = 0;
	int Gen = 1;

	// simple data
	TFile *toy_file = TFile::Open("../build/RESULTS/GENtrees/apr_ER2_RG2_W2_M50_mu5_G"+TString::Itoa(Gen,10)+"_V"+TString::Itoa(type,10)+".root");

	TString toy_name = "apr_ER2_RG2_W2_M50_mu5_G1_V0_";
	TTree *toy_tree = (TTree*)toy_file->Get( toy_name + TString::Itoa(tree_index,10) );

	// set tree branches
	Float_t cs1 =0;
	Float_t cs2 =0;
	toy_tree->SetBranchAddress("cs1", &cs1);
	toy_tree->SetBranchAddress("cs2", &cs2);
	Long64_t N = toy_tree->GetEntries();

	// likelihood
	pdfLikelihood* likeHood = getTheLikelihoodByType(50, type);
	likeHood->initialize();
	likeHood->getPOI()->setCurrentValue(1.);

	const int scan_points = 20;
	TGraph mu_scan(scan_points);
	
	double mu = -1.;
	double step = 0.1;
	

	double Nb = plotHelpers::giveNb(likeHood);
	double sig_scale = likeHood->getSignalMultiplier();

	TH2F bkg_pdf = likeHood->getOverallBkgHisto();
	bkg_pdf.Scale(1. / Nb);   // scaling to be integral 1.
	TH2F signal_pdf = likeHood->signal_component->getInterpolatedHisto();
	signal_pdf.Scale( sig_scale );// scaling to be integral 1.
 
	// moshe check 
	// cout << bkg_pdf.Integral() << "    " << signal_pdf.Integral() << endl;

	for(int j=0; j<scan_points; j++){

		double 	extended_term = 0;

		// loop over events
		for ( Long64_t i =0; i < N ; i++){

			toy_tree->GetEntry(i);

			// this is fs normalized to 1 
			double fs = signal_pdf.GetBinContent(signal_pdf.FindBin(cs1, cs2));
			// normalized to 1 fb
			double fb = bkg_pdf.GetBinContent( bkg_pdf.FindBin(cs1, cs2));

			extended_term += ( fs - fb) * Nb / ( mu*fs + Nb*fb);

		}	
		
		double formula =  N -mu -Nb + extended_term ;

		mu_scan.SetPoint(j, mu, formula);

		mu += step;

	}

	mu_scan.GetXaxis()->SetTitle("#mu");
	mu_scan.GetYaxis()->SetTitle("formula_scan");
	mu_scan.Draw("A*C");



	//=============== PART 2 =======================//
	pdfLikelihood* like_to_fit = getTheLikelihoodToFit("apr_ER2_RG2_W2", 5, 50, type, Gen);
	
	// setting all parameter to be fixed
	map <int, LKParameter*> *params = like_to_fit->getParameters();
	int itr = 0;
	for(ParameterIterator ip=params->begin(); ip!=params->end(); ip++){
		itr++;
        	if(itr == 1 ) continue;
		LKParameter *par = ip->second;
		par->setType(FIXED_PARAMETER);
	}

	// no safeguard 
	like_to_fit->setWithSafeGuard(false);
	like_to_fit->setTreeIndex(tree_index);
	
	like_to_fit->getPOI()->setMinimum(-0.5);
	like_to_fit->getPOI()->setMaximum(5.);


	TGraph* ll_scan = like_to_fit->getGraphOfLogLikelihood( 20 );
	new TCanvas();
	ll_scan->Draw("A*L");

	//like_to_fit->maximize(false);
        //cout << "MU best fit:: " << like_to_fit->getSigmaHat() << endl;   

 		
}








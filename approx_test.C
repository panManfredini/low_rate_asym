{


	double Nb = 600.;
	double Ns = 2.;


	TH1F mu_distro("mu_distro","",100, -Nb - 10.,10.);

	TRandom3 rambo;

	//loop over datasets
	for(int i=0; i< 1000; i++){

		double obs_bkg = rambo.Poisson(Nb);	
		double obs_sig = rambo.Poisson(Ns);
		double N_obs = obs_bkg + obs_sig;
		
		TGraph zero_finder(1);

		// loop on mu
		double mu = 10.;
		double step = 0.1;
		bool isPositive = (N_obs -mu -Nb -obs_bkg + Nb/mu*obs_sig) > 0 ;

		cout << N_obs -mu -Nb -obs_bkg + Nb/mu*obs_sig << " mu " << mu << " isPositive " << isPositive << endl;

		for(int j=0; j < 2000; j++){

			double current_f = N_obs -mu -Nb -obs_bkg + Nb/mu*obs_sig;

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


		mu_distro.Fill(mu_hat);
	}


	mu_distro.Draw();
}

{

	// 1) compute the pdf of:  [ fs / (mu*fs + Nb*fb) ]
	// 2) compute the integral as a function of mu:
	//   Nb * integral [ y * pdf(y) dy ] ;   where y = fs / (mu*fs + Nb*fb)
	
		

	// likelihood
	pdfLikelihood* likeHood = getTheLikelihoodByType(50, 0);
	likeHood->initialize();
	likeHood->getPOI()->setCurrentValue(1.);

	
	double Nb = plotHelpers::giveNb(likeHood);
	double sig_scale = likeHood->getSignalMultiplier();

	TH2F bkg_pdf = likeHood->getOverallBkgHisto();
	bkg_pdf.Scale(1. / Nb);   // scaling to be integral 1.
	TH2F signal_pdf = likeHood->signal_component->getInterpolatedHisto();
	signal_pdf.Scale( sig_scale );// scaling to be integral 1.



	double mu = 0.;

	
	// Exponential binning between [xmin, xmax]
	double xmin = 0.;
	double xmax = 2.;
	const int tau = 10000;
	double x_bins[tau+1]; 

	for(int	i = 0; i <= tau ; i++){
		double xbin_temp = ( TMath::Exp( ((double)i)/((double)tau))  -1. ) * ( xmax - xmin ) / (TMath::Exp( 1. ) -1.) + xmin ;
		x_bins[i] = xbin_temp; 
	}
	TH1F *pdf_of_y = new TH1F("pdf_of_y" , "", 1000000, -10., 3000.);
	//TH1D *pdf_of_y = new TH1D("pdf_of_y" , "", tau, x_bins);

	
	double mu_min = -1000.;


	double fs_max = 0.;
	double fb_atmax = 0.;
	
	// loop over histos and fill PDF
	for(int x = 1; x <= bkg_pdf.GetNbinsX(); x++){

		for(int y = 1; y <= bkg_pdf.GetNbinsY(); y++){
		
			double fs = signal_pdf.GetBinContent(x,y);
			double fb = bkg_pdf.GetBinContent(x,y);
			

			if(bkg_pdf.GetBinContent(x,y) > 1E-9  ) {
				if(mu_min < ((-1.)*Nb * fb / fs ) && fs > 1E-9 ) mu_min = (-1.)*Nb * fb / fs; 
				if(fs > fs_max) { fs_max = fs; fb_atmax = fb; }

				double content = fs / fb ;
				//double content = (fs -fb)*Nb/ (mu * fs + Nb * fb) ;
				//double content = fs / (mu * fs + Nb * fb) ;
			       // cout << " content " << content << endl;	
				pdf_of_y->Fill(content, fb );
			}
			else { cout << "MISSED!" << bkg_pdf.GetBinContent(x,y) << "    " << signal_pdf.GetBinContent(x,y) << endl; }
		}		
	}

	pdf_of_y->Draw();
	cout << "Integral : " << pdf_of_y->Integral() << endl;
	double norm_pdf = 0;
	double mean =0;
	for(int j = 1; j <= pdf_of_y->GetNbinsX(); j++){
		
		double f_x = pdf_of_y->GetBinContent(j);
		double width = pdf_of_y->GetXaxis()->GetBinWidth(j);

		double x_bin = pdf_of_y->GetXaxis()->GetBinCenter(j);
		mean += x_bin * f_x * width;
		
		norm_pdf += f_x * width;
		//cout << "x_bin " << x_bin << "  f_x " << f_x << endl;
	}
	cout << "Integral X  = " << pdf_of_y->Integral() << "   norm = " << norm_pdf  << endl;
	
	pdf_of_y->Scale(1./norm_pdf);  // rescale in order to have integral 1, we miss some entries...

	cout << "mu_min " << mu_min << endl;
	cout << "fs_max  " << fs_max << endl;
	cout << "fb_atmax  " << fb_atmax << endl;
	cout << "Nb  " << Nb << endl; 
	cout << "average " << mean / norm_pdf << endl; 
	cout << "Real average " << Nb * mean / norm_pdf << endl; 
       cout << "smallest bin size "<< 	pdf_of_y->GetXaxis()->GetBinWidth(1) << endl;


/*
	//============== PART 2 (integral computation)
	
	double mu = 0.;
	double step = 0.001;
	const int N_steps = 5000;
	TGraph *average =  new TGraph(1);

	//loop on mu
	for(int k=0; k <= N_steps; k++){

		double integral = 0;
		
		// loop over bins
		for(int i = 1; i <= pdf_of_y->GetNbinsX(); i++){

			double x = pdf_of_y->GetXaxis()->GetBinCenter(i);
			double pdf_x = pdf_of_y->GetBinContent(i);
			double width = pdf_of_y->GetXaxis()->GetBinWidth(i);
			
			integral += (x -1.) * Nb / ( mu*x +Nb) * pdf_x * width ; //  * width  // <---- I think the width is not needed here because pdf_of_y has sum_of_bins =1
										      //  but I might be wrong ;)
		}

		// cout << "Integral for mu " << mu << " is  " << integral << endl;
		
		average->SetPoint(k, mu, integral * Nb -mu);

		mu = mu - step;

	}

	new TCanvas();
		
	average->Draw("A*L");
*/
	
}





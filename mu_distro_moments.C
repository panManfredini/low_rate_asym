{

	// 1) compute the pdf of fs/fb = x
	// 2) compute the integral as a function of mu:
	// 	Nb * Integral[ (x-1) *Nb / (mu*x + Nb) Pdf(x) dx ]
	
	//================= PART 1 =============================//
		
	// likelihood
	pdfLikelihood* likeHood = getTheLikelihoodByType(50, 0);
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
	double norm_pdf = 0;
	for(int j = 1; j <= pdf_of_ratio->GetNbinsX(); j++){
		
		double f_x = pdf_of_ratio->GetBinContent(j);
		double width = pdf_of_ratio->GetXaxis()->GetBinWidth(j);

		norm_pdf += f_x * width;
	}
	cout << "Integral X  = " << pdf_of_ratio->Integral() << "   norm = " << norm_pdf  << endl;
	
	pdf_of_ratio->Scale(1./norm_pdf);  // rescale in order to have integral 1, we miss some entries...




	//============== PART 2 (integral computation)
	
	double mu = 0.;
	double step = 0.001;
	const int N_steps = 5000;
	TGraph *average =  new TGraph(1);

	//loop on mu
	for(int k=0; k <= N_steps; k++){

		double integral = 0;
		
		// loop over bins
		for(int i = 1; i <= pdf_of_ratio->GetNbinsX(); i++){

			double x = pdf_of_ratio->GetXaxis()->GetBinCenter(i);
			double pdf_x = pdf_of_ratio->GetBinContent(i);
			double width = pdf_of_ratio->GetXaxis()->GetBinWidth(i);
			
			integral += (x -1.) * Nb / ( mu*x +Nb) * pdf_x * width ; //  * width  // <---- I think the width is not needed here because pdf_of_ratio has sum_of_bins =1
										      //  but I might be wrong ;)
		}

		// cout << "Integral for mu " << mu << " is  " << integral << endl;
		
		average->SetPoint(k, mu, integral * Nb -mu);

		mu = mu - step;

	}

	new TCanvas();
		
	average->Draw("A*L");
}





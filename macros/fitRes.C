void fitRes()
{
  TFile* inFile_Na22 = TFile::Open("plots/run074-084_Ch2.root","READ");
  TH1F* h_charge_Na22 = (TH1F*)( inFile_Na22->Get("h_charge") );
  TF1* fit_511 = (TF1*)( h_charge_Na22->GetFunction("fit_511") );
  TF1* fit_1275 = (TF1*)( h_charge_Na22->GetFunction("fit_1275") );
  
	TFile* out = TFile::Open("plots/fitRes.root","RECREATE");

  TCanvas* c1 = new TCanvas();
  TH1F* hPad1 = (TH1F*)( gPad->DrawFrame(0.,0.,100.,10000.) );
  hPad1 -> SetTitle(Form(";charge [V#upoints];entries"));
  hPad1 -> Draw();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
	gStyle -> SetOptFit(11111);

  h_charge_Na22 -> Draw("same");



	//511keV fit

  TF1* fit_gaus_pol2 = new TF1("fit_SB","gaus(0)+pol2(3)",fit_511->GetParameter(1)-3.5*fit_511->GetParameter(2),fit_511->GetParameter(1)+5*fit_511->GetParameter(2));
  fit_gaus_pol2 -> SetParameter(0, fit_511->GetParameter(0));
	fit_gaus_pol2 -> SetParameter(1, fit_511->GetParameter(1));
	fit_gaus_pol2 -> SetParameter(2, fit_511->GetParameter(2));
	fit_gaus_pol2 -> SetParameter(3, 2.48397e+03);
	fit_gaus_pol2 -> SetParameter(4, -5.70163e+01);
	fit_gaus_pol2 -> SetParameter(5, -1.09253e-10);
	fit_gaus_pol2 -> SetParameter(4, -2*fit_gaus_pol2 -> GetParameter(1)*fit_gaus_pol2 -> GetParameter(5));
	fit_gaus_pol2 -> SetParLimits(5, -1000,0);
 	h_charge_Na22 -> Fit(fit_gaus_pol2,"NQRS+");
	fit_gaus_pol2 -> FixParameter(4, -2*fit_gaus_pol2 -> GetParameter(1)*fit_gaus_pol2 -> GetParameter(5));
	fit_gaus_pol2 -> SetLineColor(kYellow);
	h_charge_Na22 -> Fit(fit_gaus_pol2,"RS+");
	out -> cd();
	fit_gaus_pol2 -> Write("fit_gaus_pol2");
  //fit_gaus_pol2 -> Draw("same");


	TF1* fit_pol2 = new TF1("pol2","pol2",fit_511->GetParameter(1)-3.5*fit_511->GetParameter(2),fit_511->GetParameter(1)+5*fit_511->GetParameter(2));
  fit_pol2 -> SetParameter(0, fit_gaus_pol2->GetParameter(3));
	fit_pol2 -> SetParameter(1, fit_gaus_pol2->GetParameter(4));
	fit_pol2 -> SetParameter(2, fit_gaus_pol2->GetParameter(5));
	fit_pol2 -> SetLineColor(kGreen);
	fit_pol2 -> Draw("same");
	out -> cd();
	fit_pol2 -> Write("fit_pol2");

	TF1* fit_gaus = new TF1("Gaus","gaus",fit_511->GetParameter(1)-3.5*fit_511->GetParameter(2),fit_511->GetParameter(1)+5*fit_511->GetParameter(2));
  fit_gaus -> SetParameter(0, fit_gaus_pol2->GetParameter(0));
	fit_gaus -> SetParameter(1, fit_gaus_pol2->GetParameter(1));
	fit_gaus -> SetParameter(2, fit_gaus_pol2->GetParameter(2));
	fit_gaus -> SetLineColor(kBlue);
	fit_gaus -> Draw("same");
	out -> cd();
	fit_gaus -> Write("fit_gaus");

	//1275keV fit
	//fit_1275->GetParameter(1)-2.*fit_1275->GetParameter(2)
	TF1* fit_gaus_pol2_2 = new TF1("fit_SB_2","gaus(0)+pol2(3)",60,76);
  fit_gaus_pol2_2 -> SetParameter(0, fit_1275->GetParameter(0));
	fit_gaus_pol2_2 -> SetParameter(1, fit_1275->GetParameter(1));
	fit_gaus_pol2_2 -> SetParameter(2, fit_1275->GetParameter(2));
	/*fit_gaus_pol2_2 -> SetParameter(3, 2.48397e+03);
	fit_gaus_pol2_2 -> SetParameter(4, -5.70163e+01);
	fit_gaus_pol2_2 -> SetParameter(5, -1.09253e-10);*/
	fit_gaus_pol2_2 -> SetParameter(4, -2*fit_gaus_pol2_2 -> GetParameter(1)*fit_gaus_pol2_2 -> GetParameter(5));
	fit_gaus_pol2_2 -> SetParLimits(5, -1000,0);
 	h_charge_Na22 -> Fit(fit_gaus_pol2_2,"NQRS+");
	fit_gaus_pol2_2 -> FixParameter(4, -2*fit_gaus_pol2_2 -> GetParameter(1)*fit_gaus_pol2_2 -> GetParameter(5));
	fit_gaus_pol2_2 -> SetLineColor(kYellow);
	h_charge_Na22 -> Fit(fit_gaus_pol2_2,"RS+");
	out -> cd();
	fit_gaus_pol2_2 -> Write("fit_gaus_pol2_2");
  //fit_gaus_pol2_2 -> Draw("same");


	TF1* fit_pol2_2 = new TF1("pol2_2","pol2",60,76);
  fit_pol2_2 -> SetParameter(0, fit_gaus_pol2_2 -> GetParameter(3));
	fit_pol2_2 -> SetParameter(1, fit_gaus_pol2_2 -> GetParameter(4));
	fit_pol2_2 -> SetParameter(2, fit_gaus_pol2_2 -> GetParameter(5));
	fit_pol2_2 -> SetLineColor(kGreen);
	fit_pol2_2 -> Draw("same");
	out -> cd();
	fit_pol2_2 -> Write("fit_pol2_2");

	TF1* fit_gaus_2 = new TF1("Gaus_2","gaus",60,76);
  fit_gaus_2 -> SetParameter(0, fit_gaus_pol2_2 -> GetParameter(0));
	fit_gaus_2 -> SetParameter(1, fit_gaus_pol2_2 -> GetParameter(1));
	fit_gaus_2 -> SetParameter(2, fit_gaus_pol2_2 -> GetParameter(2));
	fit_gaus_2 -> SetLineColor(kBlue);
	fit_gaus_2 -> Draw("same");
	out -> cd();
	fit_gaus_2 -> Write("fit_gaus_2");
	
  
}


   

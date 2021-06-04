void tau()
{
  TFile* inFile_Na22 = TFile::Open("plots/run027-028_Ch2.root","READ");
  TProfile* pAvg_511 = (TProfile*)( inFile_Na22->Get("p_avgWf_511"));
  TProfile* pAvg_1275 = (TProfile*)( inFile_Na22->Get("p_avgWf_1275"));
	TFile* out = TFile::Open("plots/tau.root","RECREATE");

  TCanvas* c1 = new TCanvas();
  gPad -> SetGridx();
  gPad -> SetGridy();
	pAvg_511 -> SetMarkerColor(kBlue+2);
	pAvg_511 -> SetMarkerStyle(22);
	pAvg_511 -> SetMarkerSize(1);
	gStyle -> SetOptStat(0);
	gStyle -> SetOptFit(0);
	pAvg_511 -> Draw();

  TF1* fit_expo = new TF1("fit_expo", "expo",10,200);
 	pAvg_511 -> Fit(fit_expo,"RS+");

	/*TF1* fit_expo_rise = new TF1("fit_expo_rise", "[0]*(1-exp(-1*(x-[1])/[2]))",-40,-57);
	fit_expo_rise -> SetParameter(0, 0.8);
	fit_expo_rise -> SetParameter(1, -50);
	fit_expo_rise -> SetParameter(2, 100);
	fit_expo_rise -> SetLineWidth(5);
	fit_expo_rise -> SetLineColor(kBlue+2);
 	pAvg_511 -> Fit(fit_expo_rise,"RS+");
	//fit_expo_rise -> Draw("same");*/

	pAvg_511 -> GetXaxis() -> SetRangeUser(-100, 250);
	
	std::cout<<-1/fit_expo->GetParameter(1)<<std::endl;

	TLegend* legend = new TLegend(0.53,0.55,0.85,0.75);
  //legend->AddEntry(fit_expo_rise,Form("#tau_{R}=%f", fit_expo_rise->GetParameter(2)),"l");
	legend->SetBorderSize(0);
  legend->AddEntry(fit_expo,Form("#tau_{Decay} = %.0f ns", -1/fit_expo->GetParameter(1)),"l");
  legend->Draw("same");


	/*float max = pAvg_511 -> GetMaximum();
	std::cout<<"max "<<max<<std::endl;
	float trenta = max * 0.3;
	float settanta = max * 0.7;
	std::cout<<"trenta "<<pAvg_511 -> GetX(trenta)<<std::endl;
	float timeTrenta = pAvg_511 -> GetBinCenter(pAvg_511 -> FindBin(trenta));
	float timeSettanta = pAvg_511 -> GetBinCenter(pAvg_511 -> FindBin(settanta));
	std::cout<<"t30 "<<timeTrenta<<std::endl;
	std::cout<<"t70 "<<timeSettanta<<std::endl;*/

	//TFile* inFile_LED = TFile::Open("plots/run039_Ch2.root","READ");
	TFile* inFile_LED = TFile::Open("plots/analyzeNpe_07Apr_lowLed_6.0Vov_1stSiPM_new_Ch2.root","READ");
  TProfile* pAvg_1pe = (TProfile*)( inFile_LED->Get("p_avgWf_1pe"));
  TProfile* pAvg_3p3 = (TProfile*)( inFile_LED->Get("p_avgWf_3pe"));


  TCanvas* c2 = new TCanvas();
  gPad -> SetGridx();
  gPad -> SetGridy();
	pAvg_1pe -> SetMarkerColor(kBlue+2);
	pAvg_1pe -> SetMarkerStyle(22);
	pAvg_1pe -> SetMarkerSize(0.8);
	gStyle -> SetOptFit(0);
	pAvg_1pe -> Draw();

	gPad -> SetLeftMargin(-200);
	pAvg_1pe -> GetYaxis()->SetTitleOffset(1.3);

  //TF1* fit_expo_LED = new TF1("fit_expo_LED", "expo",45,120);
	TF1* fit_expo_LED = new TF1("fit_expo_LED", "expo",-115,-50);
 	pAvg_1pe -> Fit(fit_expo_LED,"RS+");

	TLegend* legend2 = new TLegend(0.53,0.55,0.85,0.75);
  //legend->AddEntry(fit_expo_rise,Form("#tau_{R}=%f", fit_expo_rise->GetParameter(2)),"l");
	legend2->SetBorderSize(0);
  legend2->AddEntry(fit_expo,Form("#tau_{Decay} = %.0f ns", -1/fit_expo_LED->GetParameter(1)),"l");
  legend2->Draw("same");

	//pAvg_1pe -> GetXaxis() -> SetRangeUser(20, 150);
	pAvg_1pe -> GetXaxis() -> SetRangeUser(-170, 0);


  
}


   

void drawSigmaE_laser()
{
  TFile* inFile_Na22 = TFile::Open("plots/analyzeNpe_07Apr_Na22_2.0Vov_dryCoupling_2ndSiPM_new_Ch2.root","READ");
  TH1F* h_charge_Na22 = (TH1F*)( inFile_Na22->Get("h_charge") );
  TF1* fit_511 = (TF1*)( h_charge_Na22->GetFunction("fit_511") );
  TF1* fit_1275 = (TF1*)( h_charge_Na22->GetFunction("fit_1275") );
  
  TGraphErrors* g_calib = new TGraphErrors();
  g_calib -> SetPoint(0,fit_511->GetParameter(1),511.);
  g_calib -> SetPointError(0,fit_511->GetParError(1),0.);
  g_calib -> SetPoint(1,fit_1275->GetParameter(1),1275.);
  g_calib -> SetPointError(1,fit_1275->GetParError(1),0.);

  TCanvas* c1 = new TCanvas();
  TH1F* hPad1 = (TH1F*)( gPad->DrawFrame(0.,0.,70.,1800.) );
  hPad1 -> SetTitle(Form(";charge [V#upoints];energy [keV]"));
  hPad1 -> Draw();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  g_calib -> SetMarkerStyle(20);
  g_calib -> SetMarkerSize(0.7);
  g_calib -> Draw("P,same");
  TF1* fit_calib = new TF1("fit_calib","pol1",0.,100.);
  fit_calib -> SetParameters(0,1275./fit_1275->GetParameter(1));
  fit_calib -> FixParameter(0,0.);
  g_calib -> Fit(fit_calib,"NRS");
  fit_calib -> Draw("same");
  
  
  
  std::vector<int> laserTunes;
  
  laserTunes.push_back(56);
  laserTunes.push_back(57);
  laserTunes.push_back(58);
  laserTunes.push_back(60);
  laserTunes.push_back(65);
  laserTunes.push_back(70);
  laserTunes.push_back(75);
  laserTunes.push_back(80);
  laserTunes.push_back(85);
  laserTunes.push_back(90);
  laserTunes.push_back(94);
  laserTunes.push_back(95);
  laserTunes.push_back(97);
  laserTunes.push_back(98);
  laserTunes.push_back(100);
  
  /*
  laserTunes.push_back(60);
  laserTunes.push_back(65);
  laserTunes.push_back(67);
  laserTunes.push_back(68);
  laserTunes.push_back(70);
  laserTunes.push_back(72);
  laserTunes.push_back(75);
  laserTunes.push_back(77);
  laserTunes.push_back(80);
  laserTunes.push_back(85);
  laserTunes.push_back(90);
  laserTunes.push_back(92);
  laserTunes.push_back(95);
  laserTunes.push_back(97);
  laserTunes.push_back(98);
  laserTunes.push_back(99);
  laserTunes.push_back(100);*/
  
  
  TGraphErrors* g = new TGraphErrors();
  TGraphErrors* g_new = new TGraphErrors();
  
  for(auto laserTune : laserTunes)
  {
    TFile* inFile = TFile::Open(Form("plots/analyzeNpe_07Apr_laser%03d_2.0Vov_dryCoupling_2ndSiPM_new_Ch2.root",laserTune),"READ");
    
    TH1F* h_charge = (TH1F*)( inFile->Get("h_charge") );
    TF1* fit = (TF1*)( h_charge->GetFunction("fit") );
    
    TF1* fit_new = new TF1("fit_new","gaus(0)",0.,100.);
    h_charge -> Fit(fit_new,"QNRLS+");
    
    float Q = fit -> GetParameter(1);
    float sigmaQ = fit -> GetParameter(2);
    float QErr = fit -> GetParError(1);
    float sigmaQErr = fit -> GetParError(2);
    
    g -> SetPoint(g->GetN(),Q,sigmaQ/Q);
    g -> SetPointError(g->GetN()-1,QErr,sigmaQ/Q*sqrt(pow(sigmaQErr/sigmaQ,2)+pow(QErr/Q,2)));
    
    Q = fit_new -> GetParameter(1);
    sigmaQ = fit_new -> GetParameter(2);
    QErr = fit_new -> GetParError(1);
    sigmaQErr = fit_new -> GetParError(2);
    
    g_new -> SetPoint(g_new->GetN(),Q,sigmaQ/Q);
    g_new -> SetPointError(g_new->GetN()-1,QErr,sigmaQ/Q*sqrt(pow(sigmaQErr/sigmaQ,2)+pow(QErr/Q,2)));
  }
  
  
  TCanvas* c2 = new TCanvas();
  TH1F* hPad2 = (TH1F*)( gPad->DrawFrame(0.,0.,70.,0.1) );
  hPad2 -> SetTitle(Form(";charge [V#upoints];#sigma_{Q}/Q"));
  hPad2 -> Draw();
  gPad -> SetGridx();
  gPad -> SetGridy();  
  
  // g -> SetMarkerStyle(20);
  // g -> SetMarkerSize(0.7);
  // g -> Draw("AP");

  g_new -> SetMarkerStyle(20);
  g_new -> SetMarkerSize(0.7);
  g_new -> Draw("P,same");
  
  TF1* fit = new TF1("fit","sqrt(pow([0]/x/x/x,2)+pow([1]/sqrt(x),2)+pow([2],2))",0.,100.);
  fit -> SetParameters(1.,0.1,0.005);
  fit -> FixParameter(2,0.);
  g_new -> Fit(fit,"N","",9.,100.);
  fit -> Draw("same");
  
  TF1* fit_noise = new TF1("fit","[0]/x/x/x",0.,100.);
  fit_noise -> SetParameter(0,fit->GetParameter(0));
  fit_noise -> SetLineStyle(7);
  fit_noise -> SetLineWidth(1);
  fit_noise -> SetLineColor(kRed+1);
  fit_noise -> Draw("same");
  
  TF1* fit_stoc = new TF1("fit","[0]/sqrt(x)",0.,100.);
  fit_stoc -> SetParameter(0,fit->GetParameter(1));
  fit_stoc -> SetLineStyle(7);
  fit_stoc -> SetLineWidth(1);
  fit_stoc -> SetLineColor(kRed+1);
  fit_stoc -> Draw("same");
  
  TF1* fit_const = new TF1("fit","[0]",0.,100.);
  fit_const -> SetParameter(0,fit->GetParameter(2));
  fit_const -> SetLineStyle(7);
  fit_const -> SetLineWidth(1);
  fit_const -> SetLineColor(kRed+1);
  fit_const -> Draw("same");
  
  
  TGraph* g_nPE = new TGraph();
  TGraph* g_nPE_sub = new TGraph();
  for(int point = 0; point < g_new->GetN(); ++point)
  {
    double Q = g_new -> GetPointX(point);
    double sigmaQoQ = g_new -> GetPointY(point);
    double sigmaQoQ_sub = sqrt( pow(g_new->GetPointY(point),2) - pow(fit_noise->Eval(Q),2) );
    
    g_nPE -> SetPoint(g_nPE->GetN(),fit_calib->Eval(Q)/1000.,1.03732*1.03732/sigmaQoQ/sigmaQoQ/(fit_calib->Eval(Q)/1000.));
    g_nPE_sub -> SetPoint(g_nPE_sub->GetN(),fit_calib->Eval(Q)/1000.,1.03732*1.03732/sigmaQoQ_sub/sigmaQoQ_sub/(fit_calib->Eval(Q)/1000.));
		
		
  }
  
  TCanvas* c3 = new TCanvas();
  TH1F* hPad3 = (TH1F*)( gPad->DrawFrame(0.,0.,1.8,3000) );
  hPad3 -> SetTitle(Form(";energy [MeV];N_{pe} / MeV"));
  hPad3 -> Draw();
  gPad -> SetGridx();
  gPad -> SetGridy();  
  
  g_nPE_sub -> SetMarkerStyle(20);
  g_nPE_sub -> SetMarkerSize(0.7);
  g_nPE_sub -> SetMarkerColor(kRed);
  g_nPE_sub -> Draw("P,same");
  
  g_nPE -> SetMarkerStyle(20);
  g_nPE -> SetMarkerSize(0.7);
  g_nPE -> Draw("P,same");

}


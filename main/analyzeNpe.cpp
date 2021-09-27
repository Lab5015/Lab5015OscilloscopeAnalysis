#include "interface/WFMReader.h"
#include "interface/histoFunc.h"
#include "interface/FFTAnalyzer.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"



double MultiGaussian(double* x, double* par)
{
  double xx = x[0];
  int nGaus = par[0];
  double val = 0.;
  for(int ii = 0; ii < nGaus; ++ii)
  {
    val += par[1+3*ii+0] * TMath::Gaus(xx,par[1+3*ii+1],par[1+3*ii+2]);
  }
  val += par[1+3*nGaus+0] * TMath::Gaus(xx,par[1+3*nGaus+1],par[1+3*nGaus+2]);
  
  return val;
}

double MultiGaussian_dFixed(double* x, double* par)
{
  double xx = x[0];
  int nGaus = par[0];
  double val = 0.;
  val = par[1] * TMath::Gaus(xx,par[2],par[3]);
  for(int ii = 0; ii < nGaus; ++ii)
  {
    val += par[6+2*ii] * TMath::Gaus(xx,par[4] + ii* par[5] ,par[7 +2*ii]);
  }
  return val;
}

double ExpoSum(double* x, double* par)
{
  double xx = x[0];
  double x0 = par[0];
  double C = par[1];
  double N = par[2];
  double tRise = par[3];
  double tDec = par[4];
  if( xx < x0 ) return C;
  else return C + N * ( exp(-1.*(xx-x0)/tDec) - exp(-1.*(xx-x0)/tRise) );
}




int main(int argc, char** argv)
{
  //--- parse parameters
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::vector<std::string> inFileNames = opts.GetOpt<std::vector<std::string> >("Input.inFileNames");
  std::string runType = opts.GetOpt<std::string>("Input.runType");
  std::string fitType = opts.GetOpt<std::string>("Input.fitType");
  
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  
  
  
  
  //--- define WFM reader
  std::map<std::string,int> nFrames;
  int nPoints;
  double tMin;
  double tMax;
  double tUnit;
  double VUnit;
  
	
  std::map<std::string,WFMReader*> map_reader;
  for(auto inFileName : inFileNames)
  {
    map_reader[inFileName] = new WFMReader(inFileName,true);
    nFrames[inFileName] = map_reader[inFileName]->GetNFrames();
    nPoints = map_reader[inFileName]->GetNPoints();
    tMin = map_reader[inFileName]->GetTMin();
    tMax = map_reader[inFileName]->GetTMax();
    tUnit = map_reader[inFileName]->GetTUnit();
    VUnit = map_reader[inFileName]->GetVUnit();
  }
  
  
  
  
  //--- define integration ranges
  int dynamicPedSub = opts.GetOpt<int>("Params.dynamicPedSub");
  
  double pedIntMin = opts.GetOpt<double>("Params.pedIntMin");
  double pedIntMax = opts.GetOpt<double>("Params.pedIntMax");
  double pedIntMin2 = opts.GetOpt<double>("Params.pedIntMin2");
  double pedIntMax2 = opts.GetOpt<double>("Params.pedIntMax2");
  double sigIntMin = opts.GetOpt<double>("Params.sigIntMin");
  double sigIntMax = opts.GetOpt<double>("Params.sigIntMax");
  
  int minSamplePed  = int( (pedIntMin-tMin) / tUnit);
  int maxSamplePed  = int( (pedIntMax-tMin) / tUnit);
  int minSamplePed2 = int( (pedIntMin2-tMin) / tUnit);
  int maxSamplePed2 = int( (pedIntMax2-tMin) / tUnit);
  int minSampleSig  = int( (sigIntMin-tMin) / tUnit);
  int maxSampleSig  = int( (sigIntMax-tMin) / tUnit);
  
  if( (minSamplePed < 0 || minSamplePed > nPoints) ||
      (maxSamplePed < 0 || maxSamplePed > nPoints) ||
      (minSampleSig < 0 || minSampleSig > nPoints) ||
      (maxSampleSig < 0 || maxSampleSig > nPoints) )
  {
    std::cout << "analyzeNpe::Error: integration limit(s) out of allowed range" << std::endl;
    return(-1);
  }
  
  
  
  
  //--- define outfile and histograms
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
  outFile -> cd();
  
  TH1F* h_ped;
  if( runType == "lowLED" )                                              h_ped = new TH1F("h_ped","", 512, -10.*VUnit, 5.*VUnit);
  if( runType == "Na22" || runType == "pedestal" || runType == "laser" ) h_ped = new TH1F("h_ped","",2000,-10.*VUnit,20.*VUnit);
  if( runType == "PMT" )                                                 h_ped = new TH1F("h_ped","",1000, -5.*VUnit, 0.5*VUnit);
  h_ped -> SetTitle(";pedestal charge [V#upointns];entries");
  TH1F* h_ped_sub;
  if( runType == "lowLED" )                                              h_ped_sub = new TH1F("h_ped_sub","", 512, -10.*VUnit, 10.*VUnit);
  if( runType == "Na22" || runType == "pedestal" || runType == "laser" ) h_ped_sub = new TH1F("h_ped_sub","",2000,-10.*VUnit,20.*VUnit);
  if( runType == "PMT" )                                                 h_ped_sub = new TH1F("h_ped_sub","",1000, -5.*VUnit, 0.5*VUnit);
  h_ped -> SetTitle(";pedestal charge [V#upointns];entries");  
  TH1F* h_precharge;
  if( runType == "lowLED" )                    h_precharge = new TH1F("h_precharge","", 1024, -1000.*VUnit, 5000.*VUnit);
  if( runType == "Na22" || runType == "laser") h_precharge = new TH1F("h_precharge","", 1024, 0., 40000.*VUnit);
  if( runType == "PMT")                        h_precharge = new TH1F("h_precharge","", 1024, 0., 4000.*VUnit);
  if( runType == "pedestal" )                  h_precharge = new TH1F("h_precharge","",  512, -1000*VUnit, 1000.*VUnit);
  h_precharge -> SetTitle(";charge [V#upointns];entries"); 
  TH1F* h_postcharge;
  if( runType == "lowLED" )                    h_postcharge = new TH1F("h_postcharge","", 1024, -1000.*VUnit, 5000.*VUnit);
  if( runType == "Na22" || runType == "laser") h_postcharge = new TH1F("h_postcharge","", 1024, 0., 40000.*VUnit);
  if( runType == "PMT")                        h_postcharge = new TH1F("h_postcharge","", 1024, 0., 4000.*VUnit);
  if( runType == "pedestal" )                  h_postcharge = new TH1F("h_postcharge","",  512, -1000*VUnit, 1000.*VUnit);
  h_postcharge -> SetTitle(";charge [V#upointns];entries"); 
  TH1F* h_charge;
  if( runType == "lowLED" )                    h_charge = new TH1F("h_charge","",  256*2, -1000.*VUnit, 5000.*VUnit);
  if( runType == "Na22" || runType == "laser") h_charge = new TH1F("h_charge","", 1024, 0., 50000.*VUnit);
  if( runType == "PMT")                        h_charge = new TH1F("h_charge","", 1024, 0., 4000.*VUnit);
  if( runType == "pedestal" )                  h_charge = new TH1F("h_charge","",  512, -1000*VUnit, 1000.*VUnit);
  h_charge -> SetTitle(";charge [V#upointns];entries"); 
  
  TProfile* p_avgWf_1pe = new TProfile("p_avgWf_1pe","",nPoints,tMin-0.5*tUnit,tMax+0.5*tUnit);
  p_avgWf_1pe -> SetTitle(";sample time [ns];sample value [V]");
  TProfile* p_avgWf_3pe = new TProfile("p_avgWf_3pe","",nPoints,tMin-0.5*tUnit,tMax+0.5*tUnit);
  p_avgWf_3pe -> SetTitle(";sample time [ns];sample value [V]");
  TProfile* p_avgWf_all = new TProfile("p_avgWf_all","",nPoints,tMin-0.5*tUnit,tMax+0.5*tUnit);
  p_avgWf_all -> SetTitle(";sample time [ns];sample value [V]");
  TProfile* p_avgWf_1275 = new TProfile("p_avgWf_1275","",nPoints,tMin-0.5*tUnit,tMax+0.5*tUnit);
  p_avgWf_1275 -> SetTitle(";sample time [ns];sample value [V]");
  TProfile* p_avgWf_511 = new TProfile("p_avgWf_511","",nPoints,tMin-0.5*tUnit,tMax+0.5*tUnit);
  p_avgWf_511 -> SetTitle(";sample time [ns];sample value [V]");  

  TProfile* p_avgWf_1275_norm = new TProfile("p_avgWf_1275_norm","",nPoints,tMin-0.5*tUnit,tMax+0.5*tUnit);
  p_avgWf_1275_norm -> SetTitle(";sample time [ns];sample value [V]");
  TProfile* p_avgWf_511_norm = new TProfile("p_avgWf_511_norm","",nPoints,tMin-0.5*tUnit,tMax+0.5*tUnit);
  p_avgWf_511_norm -> SetTitle(";sample time [ns];sample value [V]");  

  
  TProfile* p_avgWf_1275_norm_aligned = new TProfile("p_avgWf_1275_norm_aligned","",nPoints,tMin-0.5*tUnit,tMax+0.5*tUnit);
  p_avgWf_1275_norm_aligned -> SetTitle(";sample time [ns];sample value [V]");
  TProfile* p_avgWf_511_norm_aligned = new TProfile("p_avgWf_511_norm_aligned","",nPoints,tMin-0.5*tUnit,tMax+0.5*tUnit);
  p_avgWf_511_norm_aligned -> SetTitle(";sample time [ns];sample value [V]");  
  
  
  
  //--- 1st loop over frames
  std::vector<double> ped_values;
  
  for(auto mapIt : map_reader)
  {
    WFMReader reader = *(mapIt.second);
    for(int iFrame = 0; iFrame < nFrames[mapIt.first]; ++iFrame)
    {
      std::cout << ">>> reading frame " << iFrame << " / " << nFrames[mapIt.first] << " \r" << std::flush;
      std::pair<std::vector<double>,std::vector<double> > wf = reader.GetFrame(iFrame,0.,1.);
      //%if( iFrame%8 != 0 ) continue;
      
      double charge = 0.;
      for(int iSample = minSamplePed; iSample < maxSamplePed; ++iSample)
        charge += wf.second.at(iSample)*tUnit;
      
      ped_values.push_back( charge / double(maxSamplePed-minSamplePed) );
      
      h_ped -> Fill( charge / double(maxSamplePed-minSamplePed) );
    }
    reader.CloseFile();
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  
  //--- find pedestal
  float ped_mean = h_ped -> GetBinCenter( h_ped->GetMaximumBin() );
  float ped_rms = h_ped->GetRMS();
  int nPeaks_ped = 3;
  TSpectrum* spectrum_ped = new TSpectrum(nPeaks_ped);
  int nFound_ped = spectrum_ped -> Search(h_ped,5,"",0.01);
  double* peaks_ped;
  peaks_ped = spectrum_ped -> GetPositionX(); 
  std::vector<double> vec_peaks_ped;
  for(int ii = 0; ii < nFound_ped; ++ii)
    vec_peaks_ped.push_back(peaks_ped[ii]);
  std::sort(vec_peaks_ped.begin(),vec_peaks_ped.end()); 
  TF1* f_ped = new TF1("f_ped","gaus(0)",vec_peaks_ped.back()-0.5*ped_rms,vec_peaks_ped.back()+5.*ped_rms);
  h_ped -> Fit(f_ped,"QNRS");
  ped_mean = f_ped -> GetParameter(1);
  ped_rms = f_ped -> GetParameter(2);
  TF1* f_ped2 = new TF1("f_ped2","gaus(0)",ped_mean-1*ped_rms,ped_mean+5*ped_rms);
  h_ped -> Fit(f_ped2,"RS+");
  ped_mean = f_ped2 -> GetParameter(1);
  ped_rms = f_ped2 -> GetParameter(2);
  
  
  
  
  //--- 2nd loop over frames
  int event = 0;
  for(auto mapIt : map_reader)
  {
    WFMReader reader = *(mapIt.second);
    for(int iFrame = 0; iFrame < nFrames[mapIt.first]; ++iFrame)
    {
      std::cout << ">>> reading frame " << iFrame << " / " << nFrames[mapIt.first] << " \r" << std::flush;
      
      std::pair<std::vector<double>,std::vector<double> > wf;
      if( dynamicPedSub ) wf  = reader.GetFrame(iFrame,ped_values[event]/tUnit,-1.);
      else                wf  = reader.GetFrame(iFrame,ped_mean/tUnit,-1);
      
      //if( iFrame%8 != 0 ) continue;
      
      double charge_ped = 0.;
      for(int iSample = minSamplePed2; iSample < maxSamplePed2; ++iSample)
        charge_ped += wf.second.at(iSample)*tUnit;
			
      h_ped_sub -> Fill( charge_ped / double(maxSamplePed2-minSamplePed2) );
    
      ++event;
    }
    reader.CloseFile();
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  
  TF1* fitPedSub = new TF1( "fitPedSub", "gaus", h_ped_sub->GetMean() - 2*h_ped_sub->GetRMS(), h_ped_sub->GetMean() + 2*h_ped_sub->GetRMS());
  fitPedSub -> SetParameters(h_ped_sub->GetMaximum(), h_ped_sub->GetMean(), h_ped_sub->GetRMS()); 
  h_ped_sub -> Fit( fitPedSub, "NQR");
  TF1* fitPedSub2 = new TF1( "fitPedSub2", "gaus", fitPedSub->GetParameter(1) - 1*fitPedSub->GetParameter(2), fitPedSub->GetParameter(1) + 1*fitPedSub->GetParameter(2));
  fitPedSub2 -> SetParameters(fitPedSub->GetParameter(0), fitPedSub->GetParameter(1), fitPedSub->GetParameter(2)); 
  h_ped_sub -> Fit( fitPedSub2, "QRS+");
  
  std::cout << "Ped inf limit    " << fitPedSub2->GetParameter(1) - 3 * fitPedSub2->GetParameter(2) << std::endl;
  std::cout << "Ped  sup limit   " << fitPedSub2->GetParameter(1) + 3 * fitPedSub2->GetParameter(2) << std::endl;
  std::cout << "N events  " << h_ped_sub->Integral( h_ped_sub->FindBin(fitPedSub2->GetParameter(1) - 3 * fitPedSub2->GetParameter(2)), h_ped_sub->FindBin(fitPedSub2->GetParameter(1) + 3 * fitPedSub2->GetParameter(2))) << std::endl;
  
  
  
  
  //--- 3rd loop over frames
  std::vector<double> integral_values;
  std::vector<double> t_trigger_values;

  event = 0;
  for(auto mapIt : map_reader)
  {
    WFMReader reader = *(mapIt.second);
    for(int iFrame = 0; iFrame < nFrames[mapIt.first]; ++iFrame)
    {
      std::cout << ">>> reading frame " << iFrame << " / " << nFrames[mapIt.first] << " \r" << std::flush;
      
      std::pair<std::vector<double>,std::vector<double> > wf;
      if( dynamicPedSub ) wf  = reader.GetFrame(iFrame,ped_values[event]/tUnit,-1.);
      else                wf  = reader.GetFrame(iFrame,ped_mean/tUnit,-1);
      
      //if( iFrame%8 != 0 ) continue;
      ++event;
      
      double charge_ped = 0.;
      for(int iSample = minSamplePed2; iSample < maxSamplePed2; ++iSample)
      	charge_ped += wf.second.at(iSample)*tUnit;
      
      if( (charge_ped / double(maxSamplePed2-minSamplePed2)) < fitPedSub2->GetParameter(1) - 3 * fitPedSub2->GetParameter(2) ||
          (charge_ped / double(maxSamplePed2-minSamplePed2)) > fitPedSub2->GetParameter(1) + 3 * fitPedSub2->GetParameter(2) )
        continue;
      
      double charge = 0.;
      for(int iSample = minSampleSig; iSample < maxSampleSig; ++iSample)
          charge += wf.second.at(iSample)*tUnit; 
      h_charge -> Fill(1.*charge);

      /*
      float precharge = 0.;
      for(int iSample = minSample; iSample < maxSample; ++iSample)
        precharge += wf.second.at(iSample)*tUnit; 
      h_precharge -> Fill(precharge);
      
      float postcharge = 0.;
      minSample = std::max(0,int((sigIntMax-tMin)/tUnit));
      maxSample = std::min(nPoints,int((sigIntMin-tMin)/tUnit));
      for(int iSample = minSample; iSample < maxSample; ++iSample)
        precharge += wf.second.at(iSample)*tUnit; 
      h_precharge -> Fill(precharge);
      */
      
      //if( iFrame < 2000 && iFrame%100 == 0 )
      if( iFrame < 200 && iFrame%20 == 0 )
      {
        TGraph* g_wf = new TGraph();
        for(unsigned int iSample = 0; iSample < wf.first.size(); ++iSample)
          g_wf -> SetPoint(g_wf->GetN(),wf.first.at(iSample),wf.second.at(iSample));
        
        outFile -> cd();
        g_wf -> Write(Form("g_wf_frame%02d",iFrame));
      }
      
    }
    reader.CloseFile();
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  
  TSpectrum* spectrum;
  int nFound;
  int nPeaks;
  double* peaks;
  std::vector<double> vec_peaks;
  
  TF1* fit_singlePe;
  TF1* fit_511;
  TF1* fit_1275;
  if( runType == "lowLED" )
  {
    nPeaks = 10;
    spectrum = new TSpectrum(nPeaks);
    nFound = spectrum -> Search(h_charge,5,"",0.05);
    peaks = spectrum -> GetPositionX();
    
    
    for(int ii = 0; ii < nFound; ++ii)
      vec_peaks.push_back(peaks[ii]);
    std::sort(vec_peaks.begin(),vec_peaks.end());
    
    TF1* fit_0pe = new TF1("fit_0pe","gaus(0)",-200.*VUnit,100.*VUnit);
    fit_0pe -> SetNpx(10000);
    fit_0pe -> SetLineColor(kGreen);
    fit_0pe -> SetParameter(0,h_charge->GetBinContent(h_charge->FindBin(0.)));
    fit_0pe -> SetParameter(1,0.);
    fit_0pe -> SetParameter(2,100.*VUnit);
    h_charge -> Fit(fit_0pe,"QRLS+");
    float sigma = fit_0pe -> GetParameter(2);
    
    if( fitType == "dFree")
    {
      fit_singlePe = new TF1("fit_singlePe",MultiGaussian,-1000.*VUnit,vec_peaks.back()+2.*sigma,1+3*nFound+3);
      fit_singlePe -> SetNpx(10000);
      fit_singlePe -> FixParameter(0,nFound);
      for(int ii = 0; ii < nFound; ++ii)
      {
        fit_singlePe -> SetParameter(1+ii*3+0,h_charge->GetBinContent(h_charge->FindBin(vec_peaks[ii])));
        fit_singlePe -> SetParameter(1+ii*3+1,vec_peaks[ii]);
        fit_singlePe -> SetParameter(1+ii*3+2,sigma);
      }
    }
    
    if(fitType == "dFixed")
    {
      fit_singlePe = new TF1("fit_singlePe",MultiGaussian_dFixed,-1000.*VUnit,vec_peaks.back()+2.*sigma,3+2*nFound+3);
      fit_singlePe -> SetNpx(10000);
      fit_singlePe -> FixParameter(0,nFound);
      
      for(int l = 0; l < nFound; ++l)
      {
        TF1* fitGaus = new TF1("fit_gaus","gaus",vec_peaks[l]-0.015,vec_peaks[l]+0.015);
        fitGaus -> SetLineColor(kBlue);
        h_charge -> Fit(fitGaus, "QRS");
        fit_singlePe -> SetParameter(6+2*l, fitGaus->GetParameter(0));
        fit_singlePe -> SetParameter(7+2*l, fitGaus->GetParameter(2));
        if( l == 0 ) fit_singlePe -> SetParameter(4, fitGaus->GetParameter(1));
        if( l == 1 ) fit_singlePe -> SetParameter(5, fitGaus->GetParameter(1)-fit_singlePe->GetParameter(4));
      }
    }
    
    TH1F* h_charge_bkg = new TH1F("h_charge_bkg","",
                                  h_charge->GetNbinsX(),
                                  h_charge->GetXaxis()->GetBinLowEdge(1),
                                  h_charge->GetXaxis()->GetBinLowEdge(h_charge->GetNbinsX())+h_charge->GetXaxis()->GetBinWidth(1));
    
    for(int ii = 0; ii < nFound-1; ++ii)
    {
      float x = 0.5 * (vec_peaks[ii] + vec_peaks[ii+1]);
      float delta = vec_peaks[ii+1] - vec_peaks[ii];
      for(int bin = 1; bin <= h_charge->GetNbinsX(); ++bin)
      {
        if( h_charge->GetXaxis()->GetBinCenter(bin) > (x - 0.1*delta) &&
            h_charge->GetXaxis()->GetBinCenter(bin) < (x + 0.1*delta) )
        {
          h_charge_bkg -> Fill(h_charge->GetXaxis()->GetBinCenter(bin),h_charge->GetBinContent(bin));
        }
      }
    }
    TF1* fit_charge_bkg = new TF1("fit_charge_bkg","gaus(0)",-1000.*VUnit,vec_peaks.back()+2.*sigma);
    h_charge_bkg -> Fit(fit_charge_bkg,"QRS+");
    
    if( fitType == "dFree")
    {
      fit_singlePe -> SetParameter(1+nFound*3+0,fit_charge_bkg->GetParameter(0));
      fit_singlePe -> SetParameter(1+nFound*3+1,fit_charge_bkg->GetParameter(1));
      fit_singlePe -> SetParameter(1+nFound*3+2,fit_charge_bkg->GetParameter(2));
    }
    
    if( fitType == "dFixed")
      {
	//fit_singlePe -> SetParameter(1,fit_charge_bkg->GetParameter(0));
	//fit_singlePe -> SetParameter(2,fit_charge_bkg->GetParameter(1));
	//fit_singlePe -> SetParameter(3,fit_charge_bkg->GetParameter(2));
	fit_singlePe -> FixParameter(1,0.);
	fit_singlePe -> FixParameter(2,0.);
	fit_singlePe -> FixParameter(3,0.);
    }
    
    for(int l = 0; l < nFound; ++l)
      {
        TF1* fitGaus = new TF1("fit_gaus","gaus",vec_peaks[l]-0.015,vec_peaks[l]+0.015);
	h_charge -> Fit(fitGaus, "QNRS+");
        fitGaus -> SetLineColor(kBlue);
        fit_singlePe -> SetParameter(6+2*l, fitGaus->GetParameter(0)-fit_charge_bkg->Eval(fitGaus->GetParameter(1)));
        if( l == 0 ) fit_singlePe -> SetParameter(4, fitGaus->GetParameter(1));
        if( l == 1 ) fit_singlePe -> SetParameter(5, fitGaus->GetParameter(1)-fit_singlePe->GetParameter(4));
      }
    
    // for(int iPar = 0; iPar < 3+2*nFound+3; ++iPar)
    //   fit_singlePe -> FixParameter(iPar,fit_singlePe->GetParameter(iPar));
    
    h_charge -> Fit(fit_singlePe,"QRS");
    
    if( fitType == "dFree")
    {
      TGraphErrors* g = new TGraphErrors();
      for(int ii = 0; ii < nFound; ++ii)
      {
        g -> SetPoint(g->GetN(),ii,fit_singlePe->GetParameter(1+ii*3+1)); 
        g -> SetPointError(g->GetN()-1,0.,fit_singlePe->GetParError(1+ii*3+1));
      }
      
      TF1* fit_charge_vs_Npe = new TF1("fit_charge_vs_Npe","pol1",-1.,10.);
      g -> Fit(fit_charge_vs_Npe,"QRS+");
      
      g -> SetTitle(";N_{pe};charge [V#upointns]");
      g -> SetName("g_charge_vs_Npe");
      g -> Write();
    }
    
    if( fitType == "dFixed")
    {
      TGraphErrors* g = new TGraphErrors();
      g -> SetPoint(g->GetN(),1,fit_singlePe->GetParameter(5)); 
      g -> SetPointError(g->GetN()-1,0.,fit_singlePe->GetParError(5));
      
      g -> SetTitle(";N_{pe};charge [V#upointns]");
      g -> SetName("g_charge_vs_Npe");
      g -> Write();
    }
  }
  
  
  if( runType == "Na22" )
  {
    nPeaks = 6;
    spectrum = new TSpectrum(nPeaks);
    nFound = spectrum -> Search(h_charge, 20.0, "", 0.005);
    peaks = spectrum -> GetPositionX();
    for(int ii = 0; ii < nFound; ++ii)
      vec_peaks.push_back(peaks[ii]);
    std::sort(vec_peaks.begin(),vec_peaks.end());
    
    if( nFound > 0 )
      {
	// fit 511
	int bin_511 = h_charge->FindBin(peaks[0]);
	float sigma_511 = 0.;
	for(int bin = bin_511; bin <= h_charge->GetNbinsX(); ++bin)
	  {
	    if( h_charge->GetBinContent(bin) < 0.5*h_charge->GetBinContent(bin_511) )
	      {
		sigma_511 = (h_charge->GetBinCenter(bin) - h_charge->GetBinCenter(bin_511))/(2.35/2.);
		break;
	      }
	  }
	
	fit_511 = new TF1 ("fit_511", "gaus", peaks[0]-2.*sigma_511,peaks[0]+2.*sigma_511);
	fit_511 -> SetParameters(h_charge->GetBinContent(bin_511),peaks[0],sigma_511);
	h_charge -> Fit(fit_511, "QNRS");
	h_charge -> Fit(fit_511, "RS+", "", fit_511->GetParameter(1)-fit_511->GetParameter(2), fit_511->GetParameter(1)+fit_511->GetParameter(2));
	
	TF1* fit_511_gausSum = new TF1("fit_511_gausSum","gaus(0)+gaus(3)",peaks[0]-2.*sigma_511,peaks[0]+2.*sigma_511);
	fit_511_gausSum -> SetParameters(h_charge->GetBinContent(bin_511)*2./3.,peaks[0],sigma_511/2.,h_charge->GetBinContent(bin_511)/3.,peaks[0],3.*sigma_511);
	h_charge -> Fit(fit_511_gausSum,"RS+");
	
	fit_1275 = new TF1 ("fit_1275", "gaus", vec_peaks[nFound-1] - vec_peaks[nFound-1]*0.08, vec_peaks[nFound-1] + vec_peaks[nFound-1]*0.2);
	h_charge -> Fit(fit_1275, "QNR");
	h_charge -> Fit(fit_1275, "RS+", "", fit_1275->GetParameter(1)-1.*fit_1275->GetParameter(2), fit_1275->GetParameter(1)+1.5*fit_1275->GetParameter(2));
      }
  }
  
  
  if( runType == "laser" || runType == "pedestal" )
  {
    double min = -0.5;
    nPeaks = 1;
    spectrum = new TSpectrum(nPeaks);
    if( runType == "laser") min = 20.0;
    nFound = spectrum -> Search(h_charge, min, "", 0.005);
    peaks = spectrum -> GetPositionX();
    for(int ii = 0; ii < nFound; ++ii)
      vec_peaks.push_back(peaks[ii]);
    std::sort(vec_peaks.begin(),vec_peaks.end());
    
    // fit
    int bin_max = h_charge->FindBin(peaks[0]);
    float sigma = 0.;
    for(int bin = bin_max; bin <= h_charge->GetNbinsX(); ++bin)
    {
      if( h_charge->GetBinContent(bin) < 0.5*h_charge->GetBinContent(bin_max) )
      {
        sigma = (h_charge->GetBinCenter(bin) - h_charge->GetBinCenter(bin_max))/(2.35/2.);
        break;
      }
    }
    
    TF1* fit = new TF1 ("fit", "gaus", peaks[0]-2.*sigma,peaks[0]+2.*sigma);
    fit -> SetParameters(h_charge->GetBinContent(bin_max),peaks[0],sigma);
    h_charge -> Fit(fit, "QNRS");
    h_charge -> Fit(fit, "RS+", "", fit->GetParameter(1)-fit->GetParameter(2), fit->GetParameter(1)+fit->GetParameter(2));
  }
  
  
  
  //--- 4th loop over frames
  event = 0;
  for(auto mapIt : map_reader)
  {
    WFMReader reader = *(mapIt.second);
    
    for(int iFrame = 0; iFrame < nFrames[mapIt.first]; ++iFrame)
    {
      std::cout << ">>> reading frame " << iFrame << " / " << nFrames[mapIt.first] << " \r" << std::flush;
      
      std::pair<std::vector<double>,std::vector<double> > wf;
      if( dynamicPedSub ) wf  = reader.GetFrame(iFrame,ped_values[event]/tUnit,-1.);
      else                wf  = reader.GetFrame(iFrame,ped_mean/tUnit,-1);

      //normalize and align wf
      std::pair<std::vector<double>,std::vector<double> > wf_norm;
      std::pair<std::vector<double>,std::vector<double> > wf_norm_aligned;

      //if( iFrame%8 != 0 ) continue;
      ++event;
      double charge_ped = 0.;
      for(int iSample = minSamplePed2; iSample < maxSamplePed2; ++iSample)
      	charge_ped += wf.second.at(iSample)*tUnit;
      
      if ((charge_ped / double(maxSamplePed-minSamplePed)) < fitPedSub2->GetParameter(1) - 3 * fitPedSub2->GetParameter(2)) continue;
      
      //integrate charge and find maximum amplitude of waveform
      double charge = 0.;
      double max_amp = 0.;
      double t_trigger = 0;
      for(int iSample = minSampleSig; iSample < maxSampleSig; ++iSample)
      {
          charge += wf.second.at(iSample)*tUnit; 
	  if (wf.second.at(iSample)>max_amp) max_amp = wf.second.at(iSample);
      }
      h_charge -> Fill(1.*charge);

      //find trigger time
      for(int iSample = minSampleSig; iSample < maxSampleSig; ++iSample)
      {
	  if (wf.second.at(iSample)>0.5*max_amp) 
	    { 
	      t_trigger = wf.first.at(iSample);
	      break;
	    }
      }
      //fill normalized and aligned waveforms
      for(int iSample = 0; iSample < nPoints; ++iSample)
      {
	  wf_norm.first.push_back(wf.first.at(iSample));
	  wf_norm.second.push_back(wf.second.at(iSample)/charge);

	  wf_norm_aligned.first.push_back(wf.first.at(iSample)-t_trigger);
	  wf_norm_aligned.second.push_back(wf.second.at(iSample)/charge);
      }
   
      if( runType == "lowLED" && fabs(charge-(fit_singlePe->GetParameter(4)+1.*fit_singlePe->GetParameter(5))) < fit_singlePe->GetParameter(9) )
      {
        for(unsigned int iSample = 0; iSample < wf.first.size(); ++iSample)
          p_avgWf_1pe -> Fill(wf.first.at(iSample),wf.second.at(iSample));
      }
      if( runType == "lowLED" && fabs(charge-(fit_singlePe->GetParameter(4)+3.*fit_singlePe->GetParameter(5))) < fit_singlePe->GetParameter(13) )
      {
        for(unsigned int iSample = 0; iSample < wf.first.size(); ++iSample)
          p_avgWf_3pe -> Fill(wf.first.at(iSample),wf.second.at(iSample));
      }
      
      if( runType == "Na22" )
      {
        if( fabs(charge-fit_511->GetParameter(1)) < fit_511->GetParameter(2) )
        {
          for(unsigned int iSample = 0; iSample < wf.first.size(); ++iSample)
	    {
	      p_avgWf_511 -> Fill(wf.first.at(iSample),wf.second.at(iSample));
	      p_avgWf_511_norm -> Fill(wf_norm.first.at(iSample),wf_norm.second.at(iSample));
	      p_avgWf_511_norm_aligned -> Fill(wf_norm_aligned.first.at(iSample),wf_norm_aligned.second.at(iSample));
	    }
        }
        
        if( fabs(charge-fit_1275->GetParameter(1)) < fit_1275->GetParameter(2) )
        {
          for(unsigned int iSample = 0; iSample < wf.first.size(); ++iSample)
	    {
	      p_avgWf_1275 -> Fill(wf.first.at(iSample),wf.second.at(iSample));
	      p_avgWf_1275_norm -> Fill(wf_norm.first.at(iSample),wf_norm.second.at(iSample));
	      p_avgWf_1275_norm_aligned -> Fill(wf_norm_aligned.first.at(iSample),wf_norm_aligned.second.at(iSample));
	    }
        }
      }
      
      if( runType == "laser" )
      {
        for(unsigned int iSample = 0; iSample < wf.first.size(); ++iSample)
          p_avgWf_all -> Fill(wf.first.at(iSample),wf.second.at(iSample));
      }
    }
    std::cout << std::endl;
    
  }
  std::cout << std::endl;
  
  // TF1* fit_avgWf = new TF1("fit_avgWf",ExpoSum,tMin,tMax,5);
  // fit_avgWf -> SetNpx(10000);
  // fit_avgWf -> SetParameters(-120.,0.,0.7,10.,40);
  // p_avgWf_all -> Fit(fit_avgWf, "RS+","",-120.,-95.);
  // fit_avgWf -> Write();
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}

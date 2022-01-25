#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"


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


//remove spikes
std::pair<std::vector<double>,std::vector<double> > removeSpike(std::pair<std::vector<double>,std::vector<double> > wf)
{
  std::pair<std::vector<double>,std::vector<double> > wf_cleaned;
  
  std::vector<double> time = wf.first;
  std::vector<double> amp_cleaned;
  
  float grad = 0.012;
  bool spike_found = false;

  for(int iSample = 0; iSample < int(wf.first.size()); ++iSample)
  {
    double amp = wf.second.at(iSample);
    if (iSample > 1 && iSample < int(wf.first.size())-1)
    {
      double prev_amp = wf.second.at(iSample-1);
      double post_amp = wf.second.at(iSample+1);

      if (spike_found)//remove also second sample after spike
      {
	amp = post_amp;
	spike_found = false;
      }
    
      if (amp-prev_amp>grad && amp>0.013)
      {
	amp = prev_amp;
	spike_found = true;
	//	std::cout << "spike found" << std::endl;
      }
    }
    amp_cleaned.push_back(amp);
  }

  wf_cleaned = std::make_pair(time, amp_cleaned);

  return wf_cleaned;
  
}

//find trigger time sample
int trigTimeSample (std::pair<std::vector<double>,std::vector<double> > wf, float max_amp)
{
  int trig_t_sample = 100;
  //  std::cout << "searching for trigger time..." << std::endl;

  for(int iSample = 1; iSample < int(wf.first.size())-1; iSample++)
  {
    double amp = wf.second.at(iSample);
    if (amp>0.5*max_amp)
    {
      trig_t_sample = iSample;
      break;
    }
  }
  //std::cout << "found trig_time = " << trig_t_sample << std::endl;
  return trig_t_sample;
}


TF1 * GetFitFunc (TGraph *g, int iMaxSample)
{
  //
  TGraph * gFit =  new TGraph();

  //  --- fit  with gaus around maximum
  TF1 * fitPulse = new TF1 ("fitPulse", "gaus", 0, 1024);
  for (int iS = 0; iS <5; iS++)
  {
    Double_t t,s;
    if (iMaxSample-2+iS>0  && iMaxSample-2+iS<1024) 
    {
  	g->GetPoint(iMaxSample-2+iS, t, s);
  	gFit->SetPoint(gFit->GetN(), t, s);
    }
  }

  // TF1 * fitPulse = new TF1 ("fitPulse", ExpoSum, 0, 1024, 5);
  // fitPulse->SetParameters(iMaxSample,0.,0.7,10.,40);
  // for (int iS = -3; iS <10; iS++)
  // {
  //   Double_t t,s;
  //   if (iMaxSample+iS>0 && iMaxSample+iS<1024) 
  //   {
  // 	g->GetPoint(iMaxSample+iS, t, s);
  // 	gFit->SetPoint(gFit->GetN(), t, s);
  // 	//	std::cout << "i = " << gFit->GetN()-1 << " :: t = " << t << " :: s = " << s << std::endl;
  //   }
  // }


  gFit->Fit(fitPulse, "Q");
  //  std::cout << "fitPulse->GetMaximum() <<  " << fitPulse->GetMaximum() << std::endl;
  //  gFit->Delete();
  return fitPulse;
}



int main(int argc, char** argv)
{
  //--- parse parameters
  //  CfgManager opts;
  // opts.ParseConfigFile(argv[1]);

  //  std::string inFileName = opts.GetOpt<std::string>("Input.inFileName");
  // std::vector<int> channels = opts.GetOpt<std::vector<int> >("Input.channels");

  int  runId;// = opts.GetOpt<int>("Input.trigBar");
  if (argc>1) runId = atoi(argv[1]);

  int  trigBar = 1;// = opts.GetOpt<int>("Input.trigBar");
  if (argc>2) trigBar = atoi(argv[2]);

  bool dumpWf = false;
  if (argc>3) dumpWf = atoi(argv[3]);
  
  //  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
   std::string inFileName = Form("/storage/DT5742/run%04d_wave", runId);
  //  std::string inFileName = Form("/storage/DT5742/prova_wave", runId);
  //  std::string inFileName = Form("/storage/DT5742/prova3_wave", runId);


  //  std::string outFileName = Form("/home/cmsdaq/Programs/Lab5015OscilloscopeAnalysis/reco_lautaro/reco_run%04d.root", runId);
  std::string outFileName = Form("/home/cmsdaq/Programs/Lab5015OscilloscopeAnalysis/reco_lautaro/reco_run%04d.root", runId);
  //std::string outFileName = Form("/home/cmsdaq/Programs/Lab5015OscilloscopeAnalysis/reco_lautaro/prova_reco.root", runId);
  TFile* outFile = new TFile(Form("%s",outFileName.c_str()), "RECREATE");

  std::cout << "Processing run " << runId << " with trigBar " << trigBar << std::endl;
  TTree* tree = new TTree(Form("data"),Form("data"));
  std::map<int,TTree*> chTrees;


  //--- mapping from bar to ch
  std::map<int, std::pair<int,int>> chFromBar;
  //-- bottom side
  chFromBar[0]  = std::make_pair(3,12);
  chFromBar[1]  = std::make_pair(3,12);
  chFromBar[2]  = std::make_pair(2,13);
  chFromBar[3]  = std::make_pair(2,13);
  chFromBar[4]  = std::make_pair(1,14);
  chFromBar[5]  = std::make_pair(1,14);
  chFromBar[6]  = std::make_pair(0,15);
  chFromBar[7]  = std::make_pair(0,15);

  //-- top side
  chFromBar[8]  = std::make_pair(4,11);
  chFromBar[9]  = std::make_pair(4,11);
  chFromBar[10] = std::make_pair(5,10);
  chFromBar[11] = std::make_pair(5,10);
  chFromBar[12] = std::make_pair(6,9);
  chFromBar[13] = std::make_pair(6,9);
  chFromBar[14] = std::make_pair(7,8);
  chFromBar[15] = std::make_pair(7,8);


  //--- define histograms and other outputs
  float VUnit = 1./20.;
  float tUnit = 1E09;
  float tMinSample = 240;
  float tMaxSample = 365;

  float tMin = tMinSample*tUnit;
  float tMax = tMaxSample*tUnit;

  //histos for next-left, central and next-right bar when available
  const int NBARS = 1;
  TH1F* h_ped_L[NBARS];  // left SiPM side
  TH1F* h_ped_R[NBARS];  // right SiPM side
  // TH1F* h_ped_sub_L[NBARS];
  // TH1F* h_ped_sub_R[NBARS];
  // TH1F* h_precharge_L[NBARS];
  // TH1F* h_precharge_R[NBARS];
  // TH1F* h_postcharge_L[NBARS];
  //  TH1F* h_postcharge_R[NBARS];
  TH1F* h_charge_L[NBARS];
  TH1F* h_charge_R[NBARS];
  TH1F* h_charge_Ave[NBARS];
  TH1F* h_maxamp_L[NBARS];
  TH1F* h_maxamp_R[NBARS];
  TH1F* h_maxamp_Ave[NBARS];  
  TH2F* h2_maxAmpLR[NBARS];

  TProfile* p_avgWf_L[NBARS];
  TProfile* p_avgWf_R[NBARS];
  TProfile* p_avgWf_L_norm[NBARS];
  TProfile* p_avgWf_R_norm[NBARS];
  TProfile* p_avgWf_L_norm_aligned[NBARS];
  TProfile* p_avgWf_R_norm_aligned[NBARS];

  for (int i = 0; i<NBARS; i++)
  {
    h_ped_L[i] = new TH1F(Form("h_ped_L_%d",i),Form("h_ped_L_%d",i),4000,0.*VUnit,40.*VUnit);
    h_ped_L[i] -> SetTitle(";pedestal charge [V#upointns];entries");
    h_ped_R[i] = new TH1F(Form("h_ped_R_%d",i),Form("h_ped_R_%d",i),4000,0.*VUnit,40.*VUnit);
    h_ped_R[i] -> SetTitle(";pedestal charge [V#upointns];entries");

    // h_ped_sub_L[i] = new TH1F(Form("h_ped_sub_L_%d",i),Form("h_ped_sub_L_%d",i),2000,-10.*VUnit,2.*VUnit);
    // h_ped_sub_L[i] -> SetTitle(";pedestal charge ped sub. [V#upointns];entries");
    // h_ped_sub_R[i] = new TH1F(Form("h_ped_sub_R_%d",i),Form("h_ped_sub_R_%d",i),2000,-10.*VUnit,2.*VUnit);
    // h_ped_sub_R[i] -> SetTitle(";pedestal charge ped sub. [V#upointns];entries");

    // h_precharge_L[i] = new TH1F(Form("h_precharge_L_%d",i),Form("h_precharge_L_%d",i),1024, 0., 400.*VUnit);
    // h_precharge_L[i] -> SetTitle(";charge [V#upointns];entries");
    // h_precharge_R[i] = new TH1F(Form("h_precharge_R_%d",i),Form("h_precharge_R_%d",i),1024, 0., 400.*VUnit);
    // h_precharge_R[i] -> SetTitle(";charge [V#upointns];entries");

    // h_postcharge_L[i] = new TH1F(Form("h_postcharge_L_%d",i),Form("h_postcharge_L_%d",i),1024, 0., 400.*VUnit);
    // h_postcharge_L[i] -> SetTitle(";charge [V#upointns];entries");
    // h_postcharge_R[i] = new TH1F(Form("h_postcharge_R_%d",i),Form("h_postcharge_R_%d",i),1024, 0., 400.*VUnit);
    // h_postcharge_R[i] -> SetTitle(";charge [V#upointns];entries");

    h_charge_L[i] = new TH1F(Form("h_charge_L_%d",i),Form("h_charge_L_%d",i),1024, 0., 500.*VUnit);
    h_charge_L[i] -> SetTitle(";charge [V#upointns];entries");
    h_charge_R[i] = new TH1F(Form("h_charge_R_%d",i),Form("h_charge_R_%d",i),1024, 0., 500.*VUnit);
    h_charge_R[i] -> SetTitle(";charge [V#upointns];entries");
    h_charge_Ave[i] = new TH1F(Form("h_charge_Ave_%d",i),Form("h_charge_Ave_%d",i),1024, 0., 500.*VUnit);
    h_charge_Ave[i] -> SetTitle(";charge [V#upointns];entries");

    h_maxamp_L[i] = new TH1F(Form("h_maxamp_L_%d",i),Form("h_maxamp_L_%d",i),1024, 0., 10.*VUnit);
    h_maxamp_L[i] -> SetTitle(";max amp. [V];entries");
    h_maxamp_R[i] = new TH1F(Form("h_maxamp_R_%d",i),Form("h_maxamp_R_%d",i),1024, 0., 10.*VUnit);
    h_maxamp_R[i] -> SetTitle(";max amp. [V];entries");
    h_maxamp_Ave[i] = new TH1F(Form("h_maxamp_Ave_%d",i),Form("h_maxamp_Ave_%d",i),1024, 0., 10.*VUnit);
    h_maxamp_Ave[i] -> SetTitle(";max amp. [V];entries");
    h2_maxAmpLR[i] = new TH2F(Form("h2_maxAmpLR_%d",i),Form("h2_maxAmpLR_%d",i),1024, 0., 10.*VUnit, 1024, 0., 10.*VUnit);
    h2_maxAmpLR[i] -> SetTitle(";max amp. left [V]; max. amp. right [V]");


    p_avgWf_L[i] = new TProfile(Form("p_avgWf_L_%d", i),Form("p_avgWf_L_%d", i), 1024, 0, 1024*tUnit);
    p_avgWf_L[i] -> SetTitle(";sample time [ns];sample value [V]");  
    p_avgWf_R[i] = new TProfile(Form("p_avgWf_R_%d", i),Form("p_avgWf_R_%d", i), 1024, 0, 1024*tUnit);
    p_avgWf_R[i] -> SetTitle(";sample time [ns];sample value [V]");  

    p_avgWf_L_norm[i] = new TProfile(Form("p_avgWf_L_norm_%d", i),Form("p_avgWf_L_norm_%d", i), 1024, 0, 1024*tUnit);
    p_avgWf_L_norm[i] -> SetTitle(";sample time [ns];sample value [V]");  
    p_avgWf_R_norm[i] = new TProfile(Form("p_avgWf_R_norm_%d", i),Form("p_avgWf_R_norm_%d", i), 1024, 0, 1024*tUnit);
    p_avgWf_R_norm[i] -> SetTitle(";sample time [ns];sample value [V]");  

    p_avgWf_L_norm_aligned[i] = new TProfile(Form("p_avgWf_L_norm_aligned_%d", i),Form("p_avgWf_L_norm_aligned_%d", i), 1024, 0, 1024*tUnit);
    p_avgWf_L_norm_aligned[i] -> SetTitle(";sample time [ns];sample value [V]");  
    p_avgWf_R_norm_aligned[i] = new TProfile(Form("p_avgWf_R_norm_aligned_%d", i),Form("p_avgWf_R_norm_aligned_%d", i), 1024, 0, 1024*tUnit);
    p_avgWf_R_norm_aligned[i] -> SetTitle(";sample time [ns];sample value [V]");  

  }

  
  //  for(unsigned int iCh = 0; iCh < channels.size(); ++iCh)
  for(unsigned int iCh = 0; iCh < 16; ++iCh)
  {

      
      int ch = iCh;//channels[iCh];


      //analyze only central and next bars if available
      //      int prevBar = trigBar-1;
      //      int nextBar = trigBar+1;
      //      if (prevBar < 0) prevBar = 0;
      //      if (nextBar >15) nextBar = 15;

      // if (       ch != chFromBar[prevBar].first && ch != chFromBar[prevBar].second
      // 	      && ch != chFromBar[trigBar].first && ch != chFromBar[trigBar].second     
      // 	      && ch != chFromBar[nextBar].first && ch != chFromBar[nextBar].second ) continue;

      if ( ch != chFromBar[trigBar].first && ch != chFromBar[trigBar].second ) continue;


      std::cout << ">>> processing channel " << iCh << std::endl;

      
      //      std::ifstream* input = new std::ifstream(Form("%s_%d.dat",inFileName.c_str(),ch), std::ios::binary);

      int digiGroup = 0;
      int digiCh = ch;
      if (ch>7)
      {
	digiGroup = 1;
	digiCh -= 8;
      }
      std::cout << "left sipm of trigBar("      << trigBar << ")  is ch: " << chFromBar[trigBar].first << " (digiGr = " << digiGroup << ", digiCh = "<< digiCh 
		<< ") :: right SiPM in trigBar(" << trigBar << ") is ch: " << chFromBar[trigBar].second << " (digiGr = " << digiGroup << ", digiCh = "<< digiCh  << " )" << std::endl; 

      std::ifstream* input = new std::ifstream(Form("%s_gr%d_ch%d.dat",inFileName.c_str(),digiGroup, digiCh), std::ios::binary);
      input->seekg(0, std::ios::end); // move getter to the end of file
      int fsize = input->tellg();// get input file size
      input->seekg(0, std::ios::beg); // move getter back to the beginning
      
      
      int len, boardId, pattern, channel, evtCounter, triggerTimeTag;
      int n;
      int barId = trigBar;
      float t[1024];
      float s[1024];
      float pulseInt;
      float maxAmp;
      int maxSample=0;
      float maxAmpFit;
      float timeTrig;
      chTrees[ch] = new TTree(Form("ch%d",ch),Form("ch%d",ch));
      tree -> AddFriend(chTrees[ch]);
      
      chTrees[ch] -> Branch("evt", &evtCounter, "evt/I"); // number of samples in a waveform
      chTrees[ch] -> Branch("barId", &barId, "barId/I"); // number of samples in a waveform
      if (dumpWf)      chTrees[ch] -> Branch("t", t, "t[1024]/F"); // waveform time [n]
      if (dumpWf)      chTrees[ch] -> Branch("s", s, "s[1024]/F"); // waveform samples [n]
      chTrees[ch] -> Branch("pulseInt", &pulseInt, "pulseInt/F"); 
      chTrees[ch] -> Branch("maxAmp", &maxAmp, "maxAmp/F");
      chTrees[ch] -> Branch("maxAmpFit", &maxAmpFit, "maxAmpFit/F");
      chTrees[ch] -> Branch("timeTrig", &timeTrig, "timeTrig/F");
      chTrees[ch] -> Branch("triggerTimeTag", &triggerTimeTag, "triggerTimeTag/I"); // trigger time tag from digitizer
      
      
      int nHeaders = 6; // for DT5742
      float sampleFreq = 1E09; // DT5742 1GHz
      int ssize = 4; // length of a sample value (4 bytes per sample for DT5742)
      
      while( input->good() && input->tellg() < fsize )
	{
	  input->read(reinterpret_cast<char*>(&len),4); // size of data [bytes]
	  input->read(reinterpret_cast<char*>(&boardId),4); // board id
	  input->read(reinterpret_cast<char*>(&pattern),4); // VME specific - (meaningful only for VME boards)
	  input->read(reinterpret_cast<char*>(&channel),4); // channel id
	  input->read(reinterpret_cast<char*>(&evtCounter),4); // event id
	  input->read(reinterpret_cast<char*>(&triggerTimeTag),4); // trigger time tag
	  
	  if( evtCounter%1000 == 0 )
	  {
	      std::cout << "Processing event " << evtCounter << std::endl;
	      std::cout << "Event header: " << std::endl;
	      std::cout << "Size of data [bytes] = " << len << "   boardId = " << boardId << "  Pattern = " << pattern << "   Channel = " << channel << "  EventCounter = " << evtCounter << "  Trigger Time Tag = " << triggerTimeTag << std::endl; 
	  }
	  
	  n = (len-nHeaders*4)/ssize; // number of waveform samples = total length - header lenght


	  float baseline = 0;
	  for(int i = 0; i < n; ++i)
	  {
	    char* buffer = new char[ssize];
	    input->read(buffer,ssize);
	    t[i] = i * 1./sampleFreq * 1000000000; // ns
	    s[i] = (*(float*)&buffer[0]) * 1. / 4096.;
	      
	    if( i >= 0 && i < 100 )
	    baseline += s[i];
	  }
	  
	  baseline /= 100.;
	  pulseInt = 0.;

	  TGraph* g_wf = new TGraph();
	  maxAmp = 0;

	  for(int i = 0; i < n; ++i)
	  {
	    s[i] = -1. * (s[i]-baseline);	     
	    g_wf->SetPoint(g_wf->GetN(), t[i], s[i]);
            
	    if (ch == chFromBar[trigBar].first)  p_avgWf_L[0] -> Fill(t[i]*tUnit, s[i]);
	    if (ch == chFromBar[trigBar].second) p_avgWf_R[0] -> Fill(t[i]*tUnit, s[i]);

	    if( i >= tMinSample && i < tMaxSample )
	    {
		pulseInt += s[i];
		if (s[i]>maxAmp) 
		  {
		    maxSample = i;
		    maxAmp = s[i];
		  }
		//		if (s[i]>trigTh) triggerTime
	    }	    
	  }
	  
	  
	  TF1 * fitPulse = GetFitFunc(g_wf, maxSample);
	  maxAmpFit = fitPulse->GetMaximum();
	  if (fabs(maxAmpFit - maxAmp)/maxAmp > 2) maxAmpFit = maxAmp;

	  //	  std::cout << "maxAmp = " << maxAmp << " :: maxAmpFit = " << maxAmpFit << std::endl;


	  //--- filling histos
 	  if (ch == chFromBar[trigBar].first)  
	  {
	    //  	      std::cout << " baseline = " << baseline << std::endl;
	      h_ped_L[0]->Fill(baseline);
 	      h_maxamp_L[0]->Fill(maxAmpFit);
 	      h_charge_L[0]->Fill(pulseInt);
	  }
	  if (ch == chFromBar[trigBar].second) 
	  {
	      h_ped_R[0]->Fill(baseline);
	      h_maxamp_R[0]->Fill(maxAmpFit);
 	      h_charge_R[0]->Fill(pulseInt);
	  }

	  //	  timeTrigt = fitPulse->GetParameter(1);

	  	  	  
	  chTrees[ch]->Fill();



	  //      p_avgWf_511_norm -> Fill(wf_norm.first.at(iSample),wf_norm.second.at(iSample));
	  //	  p_avgWf_511_norm_aligned -> Fill(wf_norm_aligned.first.at(iSample),wf_norm_aligned.second.at(iSample));

	  //-- save some random waveforms
	  if( evtCounter < 1000 && evtCounter%100 == 0 )
	    {
	      outFile -> cd();
	      g_wf -> Write(Form("g_wf_ch%d_event%02d", ch, evtCounter));
	    }

	  
	}
      input->close();
  }




  //fit average waveforms
  TF1* fit_avgWf_L = new TF1("fit_avgWf_L",ExpoSum,tMin+20*tUnit,tMax-40*tUnit,5);
  fit_avgWf_L -> SetNpx(10000);
  fit_avgWf_L -> SetParameters(270.*tUnit,0.,0.7,10.*tUnit,40*tUnit);
  p_avgWf_L[0] -> Fit(fit_avgWf_L, "R");
  fit_avgWf_L -> Write();

  TF1* fit_avgWf_R = new TF1("fit_avgWf_R",ExpoSum,tMin+20*tUnit,tMax-40*tUnit,5);
  fit_avgWf_R -> SetNpx(10000);
  fit_avgWf_R -> SetParameters(270.*tUnit,0.,0.7,10.*tUnit,40*tUnit);
  p_avgWf_R[0] -> Fit(fit_avgWf_R, "R");
  fit_avgWf_R -> Write();


  // run over trees to calculate average of L+R
  int NEVENTS_L = chTrees[chFromBar[trigBar].first]->GetEntries();
  int NEVENTS_R = chTrees[chFromBar[trigBar].second]->GetEntries();



  //  float pulseInt_L;
  //  float maxAmp_L; 
  float maxAmpFit_L; 

  //  float pulseInt_R;
  //  float maxAmp_R; 
  float maxAmpFit_R;


  chTrees[chFromBar[trigBar].first]  ->SetBranchAddress("maxAmpFit", &maxAmpFit_L);    
  chTrees[chFromBar[trigBar].second] ->SetBranchAddress("maxAmpFit", &maxAmpFit_R);    


  std::cout << "NEVENTS_L = " << NEVENTS_L  << " :: NEVENTS_R = " << NEVENTS_R << std::endl;
  if (NEVENTS_L == NEVENTS_R)
  {
    for (int iEvt = 0; iEvt < NEVENTS_L; iEvt++)
    {
      //
      chTrees[chFromBar[trigBar].first]  ->GetEntry(iEvt);
      chTrees[chFromBar[trigBar].second] ->GetEntry(iEvt);

      float maxAmpFit_Ave = 0.5*(maxAmpFit_L+maxAmpFit_R);

      //      std::cout << "maxAmp_Ave = " << maxAmpFit_Ave << std::endl;
      h_maxamp_Ave[0]->Fill(maxAmpFit_Ave);
      h2_maxAmpLR[0]->Fill(maxAmpFit_L, maxAmpFit_R);
    }    
  }




  //--- fit histos

  TSpectrum* spectrum;
  int nFound;
  int nPeaks;
  double* peaks;
  std::vector<double> vec_peaks;
  
  TF1* fit_511[3];
  TF1* fit_1275[3];
  std::string cap[3]={"L","R","Ave"};

  //  TH1F* hHistoToFit = (TH1F*) h_maxamp_Ave[0];

  h_maxamp_L[0]->Rebin(2);
  h_maxamp_R[0]->Rebin(2);
  h_maxamp_Ave[0]->Rebin(2);

  TH1F* hHistoToFit[3];
  hHistoToFit[0] = h_maxamp_L[0];
  hHistoToFit[1] = h_maxamp_R[0];
  hHistoToFit[2] = h_maxamp_Ave[0];
  // = h_maxamp_Ave[0];
  for (int i = 0; i<3; i++)
  {
      nPeaks = 6;
      spectrum = new TSpectrum(nPeaks);
    //    nFound = spectrum -> Search(h_charge, 20.0, "", 0.005);
      nFound = spectrum -> Search(hHistoToFit[i], 0.05, "", 0.005);
      peaks = spectrum -> GetPositionX();
      for(int ii = 0; ii < nFound; ++ii)
      {
	vec_peaks.push_back(peaks[ii]);
	std::cout << "found peak[" << ii << "] = " << peaks[ii] << std::endl;
      }
      std::sort(vec_peaks.begin(),vec_peaks.end());
    
      if( nFound > 0 )
      {
	// fit 511
	float best_peak_511 = vec_peaks[0];
	float min_dist = 99999;
	for (int ip = 0; ip < nFound; ip++)
	{
	    if (fabs(vec_peaks[nFound-1]/1.275 - vec_peaks[ip]/0.511) < min_dist ) 
	    {
	      min_dist = fabs(vec_peaks[nFound-1]/1.275 - vec_peaks[ip]/0.511);
	      best_peak_511 = vec_peaks[ip];
	    }
	}

	std::cout << " best peak 511 = " << best_peak_511 << std::endl;
	int bin_511 = hHistoToFit[i]->FindBin(best_peak_511);
 
	float sigma_511 = 0.;
	for(int bin = bin_511; bin <= hHistoToFit[i]->GetNbinsX(); ++bin)
	{
	  if( hHistoToFit[i]->GetBinContent(bin) < 0.5*hHistoToFit[i]->GetBinContent(bin_511) )
	  {
	    sigma_511 = (hHistoToFit[i]->GetBinCenter(bin) - hHistoToFit[i]->GetBinCenter(bin_511))/(2.35/2.);
	    break;
	  }
	}
      
	fit_511[i] = new TF1 (Form("fit_511_%s", cap[i].c_str()), "gaus", best_peak_511-1.*best_peak_511*0.05,best_peak_511+1.*best_peak_511*0.05);
	fit_511[i] -> SetLineColor(kGreen+1);

	//	fit_511 -> SetParameters(hHistoToFit->GetBinContent(bin_511),best_peak_511,sigma_511);
	fit_511[i] -> SetParameters(hHistoToFit[i]->GetBinContent(bin_511),best_peak_511,best_peak_511*0.05);
	std::cout << "fit range: " <<  best_peak_511 << "+/-" << best_peak_511*0.05 << std::endl;
	hHistoToFit[i] -> Fit(fit_511[i], "QNRS");
	std::cout << "fit range: " <<  fit_511[i]->GetParameter(1) << "+/-" << fit_511[i]->GetParameter(2) << std::endl;
	hHistoToFit[i] -> Fit(fit_511[i], "RS+", "", fit_511[i]->GetParameter(1)-fit_511[i]->GetParameter(2)*0.5, fit_511[i]->GetParameter(1)+fit_511[i]->GetParameter(2)*0.5);


//---------------- 
	fit_511[i] -> SetParLimits(0, 0, 2000);
	fit_511[i] -> SetParLimits(1, 0.01, 60);
	fit_511[i] -> SetParLimits(2, 0., 0.2);
//-------------------


	float peak511 = fit_511[i]->GetParameter(1);
	float sigma511 = fit_511[i]->GetParameter(2);
	
	TF1* fit_511_gausSum = new TF1("fit_511_gausSum","gaus(0)+gaus(3)",best_peak_511-1.*sigma_511,best_peak_511+1.*sigma_511);
	fit_511_gausSum -> SetParameters(hHistoToFit[i]->GetBinContent(bin_511)*2./3.,best_peak_511,sigma_511/2.,hHistoToFit[i]->GetBinContent(bin_511)/3.,best_peak_511,3.*sigma_511);
	fit_511_gausSum->SetLineColor(kYellow+2);
	hHistoToFit[i] -> Fit(fit_511_gausSum,"RS+");

//	hHistoToFit -> Fit(fit_511,"RS+");
	
	fit_1275[i] = new TF1 (Form("fit_1275_%s", cap[i].c_str()), "gaus", peak511/0.511*1.275 - sigma511/0.511*1.275, peak511/0.511*1.275 + sigma511/0.511*1.275);
	fit_1275[i] ->SetLineColor(kBlue+1);

	hHistoToFit[i] -> Fit(fit_1275[i], "QNR");
	hHistoToFit[i] -> Fit(fit_1275[i], "RS+", "", fit_1275[i]->GetParameter(1)-1*fit_1275[i]->GetParameter(2), fit_1275[i]->GetParameter(1)+1.5*fit_1275[i]->GetParameter(2));
      }
  }





  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  outFile -> Close();
}

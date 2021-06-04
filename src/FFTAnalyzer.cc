#include "interface/FFTAnalyzer.h"



FFTClass* T2F(TGraph* g_wf)
{
  int n_samples = g_wf -> GetN();
  
  double iRe[n_samples], iIm[n_samples], oRe[n_samples], oIm[n_samples];
  auto fftr2c = TVirtualFFT::FFT(1, &n_samples, "C2CF M");
  
  for(int point = 0; point < n_samples; ++point)
  {
    iRe[point] = g_wf -> GetPointY(point);
    iIm[point] = 0.;
  }
  
  fftr2c->SetPointsComplex(iRe, iIm);
  fftr2c->Transform();
  fftr2c->GetPointsComplex(oRe, oIm);
  
  for(int i = 0; i < n_samples; ++i)
  {
    oRe[i] /= n_samples;
    oIm[i] /= n_samples;
  }
  
  FFTClass* FFT = new FFTClass();
  FFT -> SetPointsComplex(n_samples, oRe, oIm);
  return FFT;
}

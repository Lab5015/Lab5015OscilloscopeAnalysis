#! /usr/bin/python
import sys
import ROOT

'''
tunes = [10., 75., 90.]

keys = []
runs = {}
runs['277'] = [75.,  30.,194.65]
keys.append('277')
runs['271'] = [75.,  50.,179.90]
keys.append('271')
runs['272'] = [75., 100.,176.00]
keys.append('272')
runs['273'] = [75., 300.,176.27]
keys.append('273')
runs['274'] = [75., 500.,171.47]
keys.append('274')
runs['275'] = [75.,1000.,164.30]
keys.append('275')
runs['276'] = [75.,3000.,134.85]
keys.append('276')

runs['279'] = [90.,  30.,88.92]
keys.append('279')
runs['280'] = [90.,  50.,84.91]
keys.append('280')
runs['281'] = [90., 100.,81.49]
keys.append('281')
runs['282'] = [90., 300.,77.86]
keys.append('282')
runs['283'] = [90., 500.,75.51]
keys.append('283')
runs['284'] = [90.,1000.,70.94]
keys.append('284')
runs['285'] = [90.,3000.,55.15]
keys.append('285')

runs['286'] = [10.,  30.,278.15]
keys.append('286')
runs['287'] = [10.,  50.,274.89]
keys.append('287')
runs['288'] = [10., 100.,269.69]
keys.append('288')
runs['289'] = [10., 300.,266.27]
keys.append('289')
runs['290'] = [10., 500.,259.31]
keys.append('290')
runs['291'] = [10.,1000.,244.49]
keys.append('291')
runs['292'] = [10.,3000.,202.46]
keys.append('292')
'''

tunes = [10., 70., 95.]

keys = []
runs = {}

runs = { '293' : [ 10,   30,  850.71], 
         '294' : [ 10,   50,  843.87], 
         '295' : [ 10,  100,  839.97], 
         '296' : [ 10,  300,  833.03], 
         '297' : [ 10,  500,  825.14], 
         '298' : [ 10, 1000,  767.16], 
         '299' : [ 10, 3000,  611.94], 
         
         '300' : [ 70,   30,  421.83], 
         '301' : [ 70,   50,  421.29], 
         '302' : [ 70,  100,  430.23], 
         '303' : [ 70,  300,  415.59], 
         '304' : [ 70,  500,  385.99], 
         '305' : [ 70, 1000,  354.03], 
         '306' : [ 70, 3000,  295.31], 

         '307' : [ 95,   30,  214.17], 
         '308' : [ 95,   50,  209.83], 
         '309' : [ 95,  100,  209.02], 
         '310' : [ 95,  300,  202.51], 
         '311' : [ 95,  500,  200.23], 
         '312' : [ 95, 1000,  191.37], 
         '313' : [ 95, 3000,  89.81]
   }

for i in range(293, 314):
    keys.append(str(i))


charge1pe = 0.0318 # from run 269


g_Npe_vs_rate_current = {}
g_Npe_vs_rate_charge = {}
g_Npe_vs_rate_ratio = {}

for key in keys:
    run = key
    inFile = ROOT.TFile.Open('/home/cmsdaq/OscilloscopeAnalysis/plots/run'+str(run)+'_Ch2.root')
    h_charge = inFile.Get("h_charge")
    fit = h_charge.GetFunction("fit")
    
    tune = runs[run][0]
    rate = runs[run][1]
    Npe_current = runs[run][2]
    
    if tune not in g_Npe_vs_rate_current:
        g_Npe_vs_rate_current[tune] = ROOT.TGraphErrors()
        g_Npe_vs_rate_charge[tune] = ROOT.TGraphErrors()
        g_Npe_vs_rate_ratio[tune] = ROOT.TGraphErrors()
    g_Npe_vs_rate_current[tune].SetPoint(g_Npe_vs_rate_current[tune].GetN(),rate,runs[run][2])
    g_Npe_vs_rate_charge[tune].SetPoint(g_Npe_vs_rate_charge[tune].GetN(),rate,fit.GetParameter(1)/charge1pe)
    g_Npe_vs_rate_ratio[tune].SetPoint(g_Npe_vs_rate_ratio[tune].GetN(),rate,runs[run][2]/(fit.GetParameter(1)/charge1pe))


c1 = ROOT.TCanvas('c1','c1',1200,600)
c1.Divide(2,1)

c1.cd(1)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogx()
hPad1 = ROOT.gPad.DrawFrame(10.,0.,10000.,1000.)
hPad1.SetTitle(";rate [kHz];N_{pe}")
hPad1.Draw()

it = 0
for tune in tunes:
    g_Npe_vs_rate_current[tune].SetMarkerStyle(20)
    g_Npe_vs_rate_current[tune].SetMarkerColor(1+it)
    g_Npe_vs_rate_current[tune].SetLineColor(1+it)
    g_Npe_vs_rate_current[tune].Draw("PL,same")
    g_Npe_vs_rate_charge[tune].SetMarkerStyle(24)
    g_Npe_vs_rate_charge[tune].SetMarkerColor(1+it)
    g_Npe_vs_rate_charge[tune].SetLineColor(1+it)
    g_Npe_vs_rate_charge[tune].SetLineStyle(7)
    g_Npe_vs_rate_charge[tune].Draw("PL,same")
    it += 1

c1.cd(2)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogx()
hPad2 = ROOT.gPad.DrawFrame(10.,0.9,10000.,1.3)
hPad2.SetTitle(";rate [kHz];N_{pe} current/charge ratio")
hPad2.Draw()

it = 0
for tune in tunes:
    g_Npe_vs_rate_ratio[tune].SetMarkerStyle(20)
    g_Npe_vs_rate_ratio[tune].SetMarkerColor(1+it)
    g_Npe_vs_rate_ratio[tune].SetLineColor(1+it)
    g_Npe_vs_rate_ratio[tune].Draw("PL,same")
    it += 1

raw_input('ok?')

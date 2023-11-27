from ROOT import TCanvas, TGraph, TH1F, TF1, TRandom3,TFile, gPad, TLegend, gROOT, gStyle, TH2F, TGraphErrors, TLatex, gPad
from array import array
from tabulate import tabulate
import pandas as pd
import numpy as np
import math

# First tru to use def in build flux 
def flux_gfm_v1( x, par):
    """"
    from ROOT import TCanvas,TH1F,TFile,TGraph, TF1
    from ROOT import gROOT, gBenchmark, gRandom, gSystem
    from array import array
    """
    Enu = x [0]
    normfact = par[0]
    isotope  = par[1]

        #-- To check the normalization factors set the following to "true" and
        #execute the macro build_flux_gfm.C
        #  > root -l 
        #  root [0] .x build_flux_gfm.C
        #The integrals of each flux component will be printed to the screen. 
        #Use these to define norm1 - 5 and normalize to the appropriate neutrino 
        #yield per fission. 

    #chk_FluxNormFactors = False

    # - Normalization factors.

    norm1 = 0.56*6.14 * 1.0/5.671342
    norm2 = 0.08*7.08 * 1.0/6.603334
    norm3 = 0.30*5.58 * 1.0/5.081123
    norm4 = 0.06*6.42 * 1.0/5.932303
    norm5 = 0.60*2.00 * 1.0/4.897445e-02

    #if (chk_FluxNormFactors == True):
    #    norm1 = 1.0
    #    norm2 = 1.0
    #    norm3 = 1.0
    #    norm4 = 1.0
    #    norm5 = 1.0
        
    # Fisile Isotopes Fluxes below 2 MeV
    flux_E     = array('d', [0.0, 7.813e-3, 1.563e-2, 3.12e-2, 6.25e-2, 0.125, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0])

    flux_U235  = array('d', [0, 0.024, 0.092, 0.35, 0.61, 1.98, 2.16, 2.66, 2.66, 2.41, 1.69, 1.26])
    flux_U238  = array('d', [0, 0.089, 0.35 , 1.32, 0.65, 2.02, 2.18, 2.91, 2.96, 2.75, 1.97, 1.50])
    flux_Pu239 = array('d', [0, 0.14 , 0.56 , 2.13, 0.64, 1.99, 2.08, 2.63, 2.58, 2.32, 1.48, 1.08])
    flux_Pu241 = array('d', [0, 0.2  , 0.79 , 3.0 , 0.59, 1.85, 2.14, 2.82, 2.90, 2.63, 1.75, 1.32])

    # Normalize to proper neutrino yield 
    for i in range(0,12):
        flux_U235[i]  = flux_U235[i] *norm1
        flux_U238[i]  = flux_U238[i] *norm2
        flux_Pu239[i] = flux_Pu239[i]*norm3
        flux_Pu241[i] = flux_Pu241[i]*norm4

    #Graph sybtaxis in pyROOt 
    #gr = TGraph (n , x , y) n is lens of array,x is the array x, y is an array y gr is the name of the TGraph
    flux_gr_Ele2Mev_U235  = TGraph(12, flux_E, flux_U235)
    flux_gr_Ele2Mev_U238  = TGraph(12,flux_E,flux_U238 )
    flux_gr_Ele2Mev_Pu239 = TGraph(12,flux_E,flux_Pu239)
    flux_gr_Ele2Mev_Pu241 = TGraph(12,flux_E,flux_Pu241)
    flux_gr_Ele2Mev_fiTotal = TGraph(12)
    for i in range(0,12):
        flux_gr_Ele2Mev_fiTotal.SetPoint(i, flux_E[i], flux_U235[i]+flux_U238[i]+flux_Pu239[i]+flux_Pu241[i])

    #Fissile Isotope Fluxes above 2 Mev
    flux_par_Ege2MeV_U235 = TF1("flux_par_Ege2MeV_U235", "[4]*[3]*exp([0]+[1]*x+[2]*x^2)",2,12)
    flux_par_Ege2MeV_U235.SetParameter(0,  0.87 )#(0, 1.260)
    flux_par_Ege2MeV_U235.SetParameter(1, -0.160)
    flux_par_Ege2MeV_U235.SetParameter(2, -0.091)
    flux_par_Ege2MeV_U235.SetParameter(3,  norm1)
    flux_par_Ege2MeV_U235.SetParameter(4,1.04410)# Smooth factor for U235 

    flux_par_Ege2MeV_U238 = TF1("flux_par_Ege2MeV_U238","[4]*[3]*exp([0]+[1]*x+[2]*x^2)",2,12)
    flux_par_Ege2MeV_U238.SetParameter(0, 0.976)#(0, 1.500)
    flux_par_Ege2MeV_U238.SetParameter(1,-0.162)
    flux_par_Ege2MeV_U238.SetParameter(2,-0.079)
    flux_par_Ege2MeV_U238.SetParameter(3,norm2)
    flux_par_Ege2MeV_U238.SetParameter(4,1.06710)# Smooth factor for U238

    flux_par_Ege2MeV_Pu239 = TF1("flux_par_Ege2MeV_Pu239","[4]*[3]*exp([0]+[1]*x+[2]*x^2)",2,12)
    flux_par_Ege2MeV_Pu239.SetParameter(0, 0.896)#(0, 1.0800)
    flux_par_Ege2MeV_Pu239.SetParameter(1,-0.2390)
    flux_par_Ege2MeV_Pu239.SetParameter(2,-0.0981)
    flux_par_Ege2MeV_Pu239.SetParameter(3,norm3)
    flux_par_Ege2MeV_Pu239.SetParameter(4,1.05006)# Smooth factor for Pu239

    flux_par_Ege2MeV_Pu241 = TF1("flux_par_Ege2MeV_Pu241","[4]*[3]*exp([0]+[1]*x+[2]*x^2)",2,12)
    flux_par_Ege2MeV_Pu241.SetParameter(0, 0.793 )#(0, 1.3200)
    flux_par_Ege2MeV_Pu241.SetParameter(1,-0.0800)
    flux_par_Ege2MeV_Pu241.SetParameter(2,-0.1085)
    flux_par_Ege2MeV_Pu241.SetParameter(3,norm4)
    flux_par_Ege2MeV_Pu241.SetParameter(4,1.07561)# Smooth factor for Pu241

    flux_par_Ege2MeV_Total = TF1("flux_par_Ege2MeV_Total","flux_par_Ege2MeV_U235+flux_par_Ege2MeV_U238+flux_par_Ege2MeV_Pu239+flux_par_Ege2MeV_Pu241",2,12)

    # Neutron capture flux, digitized from TEXONO, Phys. Rev.D 75 021001 (2007)
    n_capture = open("ncaptureSpectrum-texono.csv",'r')
    nline = 0
    flux_gr_Ele2MeV_ncU238v2 = TGraph()
    for line in n_capture:
        field = line.strip().split()
        ee = float(field[0])
        ff = float(field[1])
        # output field[0] is the first value in file and field[1] is the second value in file
        flux_gr_Ele2MeV_ncU238v2.SetPoint(nline, ee, norm5*pow(10,ff))
        nline = nline + 1

    # n_capture = open('ncaptureSpectrum-texono.csv', 'r')
    # nline = 0
    # flux_gr_Ele2MeV_ncU238v2 = TGraph()
    # for line in n_capture:
    #     ee = float(line.split(' ')[0])
    #     ff = float(line.split(' ')[1])
    #     flux_gr_Ele2MeV_ncU238v2.SetPoint(nline, ee, norm5*pow(10,ff))
    #     nline = nline + 1

    ff = 0
    
    if isotope == 0: #Total flux (fission + n-capture)
        if Enu < 1.3:
            ff = flux_gr_Ele2MeV_ncU238v2.Eval(Enu)+flux_gr_Ele2Mev_fiTotal.Eval(Enu)
        elif Enu < 2:
            ff = flux_gr_Ele2Mev_fiTotal.Eval(Enu)
        else:
            ff = flux_par_Ege2MeV_Total.Eval(Enu)
    elif isotope == 1: #Flux from U235 fission
        if Enu < 2:
            ff = flux_gr_Ele2Mev_U235.Eval(Enu)
        else:
            ff = flux_par_Ege2MeV_U235.Eval(Enu)
    elif isotope == 2: #Flux from U238 fisiion
        if Enu < 2:
            ff = flux_gr_Ele2Mev_U238.Eval(Enu)
        else:
            ff = flux_par_Ege2MeV_U238.Eval(Enu)
    elif isotope == 3: #Flux from Pu239 fission
        if Enu < 2:
            ff = flux_gr_Ele2Mev_Pu239.Eval(Enu)
        else:
            ff = flux_par_Ege2MeV_Pu239.Eval(Enu)
    elif isotope == 4: #Flux from Pu241 fission
        if Enu < 2:
            ff = flux_gr_Ele2Mev_Pu241.Eval(Enu)
        else:
            ff = flux_par_Ege2MeV_Pu241.Eval(Enu)
    elif isotope == 5: #Flux from n-capture on U238
        if Enu < 1.2:
            ff = flux_gr_Ele2MeV_ncU238v2.Eval(Enu)
        else:
            ff = 0.0
    elif isotope == 6:
        if Enu < 2:
            ff = flux_gr_Ele2Mev_fiTotal.Eval(Enu)
        else:
            ff = flux_par_Ege2MeV_Total.Eval(Enu)

    return ff 



#Begin build flux in py
#gROOT.ProcessLine(".L /home/jcesar/UAM/PT/AntineutrinoFlux/Scripts/pyScripts/flux_gfm_v1.py")

f0 = TF1("f0", flux_gfm_v1,  0, 12, 2)
f1 = TF1("f1", flux_gfm_v1,  0, 12, 2)
f2 = TF1("f2", flux_gfm_v1,  0, 12, 2)
f3 = TF1("f3", flux_gfm_v1,  0, 12, 2)
f4 = TF1("f4", flux_gfm_v1,  0, 12, 2)
f5 = TF1("f5", flux_gfm_v1,  0, 12, 2)
f6 = TF1("f6", flux_gfm_v1,  0, 12, 2)

f0.SetParameters(1,0) #Total
f1.SetParameters(1,1) #U235
f2.SetParameters(1,2) #U238
f3.SetParameters(1,3) #Pu239
f4.SetParameters(1,4) #Pu241
f5.SetParameters(1,5) #n -Cap
f6.SetParameters(1,6) #Tot - fis 

np = 5000

x = array('d', [0.0]*np)
w = array('d', [0.0]*np)
#print(x, w)

f1.CalcGaussLegendreSamplingPoints(np, x, w, 1e-15)
f2.CalcGaussLegendreSamplingPoints(np, x, w, 1e-15)
f3.CalcGaussLegendreSamplingPoints(np, x, w, 1e-15)
f4.CalcGaussLegendreSamplingPoints(np, x, w, 1e-15)
f5.CalcGaussLegendreSamplingPoints(np, x, w, 1e-15)
f6.CalcGaussLegendreSamplingPoints(np, x, w, 1e-15)

if1 = f1.IntegralFast(np, x, w, 0.0, 8.0)
#if1 = f1.IntegralFast(np-1, x, w, 0, 8,0,9.9E-13)
if2 = f2.IntegralFast(np, x, w, 0, 8)
if3 = f3.IntegralFast(np, x, w, 0, 8)
if4 = f4.IntegralFast(np, x, w, 0, 8)
if5 = f5.IntegralFast(np, x, w, 0, 8)
if6 = f6.IntegralFast(np, x, w, 0, 8)

print("if1 =", if1 )
print("if2 =", if2 )
print("if3 =", if3 )
print("if4 =", if4 )
print("if5 =", if5 )


f0.SetNpx(600)
f1.SetNpx(600)
f2.SetNpx(600)
f3.SetNpx(600)
f4.SetNpx(600)
f5.SetNpx(400)
f6.SetNpx(600)


#---------Style-----------
gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)
#---------Style-----------

#-Scaled version of f5 (n-Cap), for drawing
sf = 1/20.
f5sg = TGraph()
fin = f5.GetNpx()
for i in range(0,fin):
    xx = 0.0 + i*(1.3-0.)/(200-1)
    f5sg.SetPoint(i, xx, sf*f5.Eval(xx))
    
f5sg.SetPoint(fin, 1.3, 0)

f0.SetLineColor( 1)
f0.SetLineStyle( 1)
f0.SetLineWidth( 2)

f1.SetLineColor(50)
f1.SetLineStyle( 3)
f1.SetLineWidth( 1)

f2.SetLineColor(50)
f2.SetLineStyle( 4)
f2.SetLineWidth( 1)

f3.SetLineColor( 9)
f3.SetLineStyle( 3)
f3.SetLineWidth( 1)

f4.SetLineColor( 9)
f4.SetLineStyle( 4)
f4.SetLineWidth( 1)

f5.SetLineColor( 1)
f5.SetLineStyle( 1)
f5.SetLineWidth( 1)

f5sg.SetLineColor( 1)
f5sg.SetLineStyle( 1)
f5sg.SetLineWidth( 1)

f6.SetLineColor( 1)
f6.SetLineStyle( 2)
f6.SetLineWidth( 2)

flux_frame = TH2F("flux_frame","",1,-0.2,9.2,1,2e-4,10)
flux_frame.GetXaxis().SetTitle("E_{#bar{#nu}_{e}} (MeV)")
flux_frame.GetYaxis().SetTitle("dN_{#bar{#nu}_{e}}/dE_{#bar{#nu}_{e}} ( #bar{#nu}_{e} MeV^{-1} fis^{-1})")
flux_frame.GetYaxis().SetTitleOffset(1.3)

leg = TLegend(0.6,0.57,0.9,0.87)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(f0,"Total flux","l")
leg.AddEntry(f6,"fisile isotopes","l")
leg.AddEntry(f1,"^{235}U  fission","l")
leg.AddEntry(f2,"^{238}U  fission","l")
leg.AddEntry(f3,"^{239}Pu fission ","l")
leg.AddEntry(f4,"^{241}Pu fission","l")
leg.AddEntry(f5sg,"^{238}U n-capture (#times1/20)","l")

canv = TCanvas("canv", "", 650, 450)
canv.cd()

gPad.SetLeftMargin(0.12)
gPad.SetRightMargin(0.08)
gPad.SetLogy(1)
gPad.SetTicks(1,1)
flux_frame.Draw()
f1.Draw(" same")
f2.Draw(" same")
f3.Draw(" same")
f4.Draw(" same")
f5sg.Draw("C same")
f6.Draw(" same")
f0.Draw(" same")
leg.Draw()

canv.Print("flux-gfm_smooth.pdf")

input()
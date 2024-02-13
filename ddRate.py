import numpy
from ROOT import TCanvas,TH1F,TFile,TGraph, TF1, TMath, TF2, TF12, TLegend
from ROOT import gROOT, gBenchmark, gRandom, gSystem, gStyle, gPad
from array import array
from math import sqrt, exp, sin, cos

def flux_gfm_v1( x, par):

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

#Some physical constants

c_pi        = 3.141592653589793238
DeltaMnp    = 1.2933 #MeV
mel         = 0.511 #MeV
N_avo       = 6.022e23
MeV_per_amu = 931.494061
G_F         = 1.166378e-11 #(Mev^-2)
hbarc       = 197.326e-13   #Mev cm
sinsqqw     = 0.231 #(sin qw)^2 approx
Si_amu      = 28.085  # g / mol
M_Si        = 2.91085247e-13 #MeV Mass of the Si in MeV
knn = (G_F**2) / (2*c_pi)
Z = Si_amu -14
N = 14


#Function for Nt in Double-differential Rate 
def Nt(Mdet):
    N_t   = (Mdet * N_avo)/ Si_amu

    return N_t

#Function for Form Factor 
def formfact(x, par):
       
    Er    = x[0]
    kappa = sqrt(2*M_Si*Er) / (hbarc*1e13)
    s     = 1.0
    R     = 1.2*pow(28,1.0/3.0)
    r     = sqrt(R**2 - 5*(s**2))
    ff    = 3*exp(- 0.5 * kappa**2 * s**2) * (sin(kappa*r) - kappa*r*cos(kappa*r)) / (kappa*r)**3

    return ff
"""
#Function of dSigms in dT  
def dSigdT_CEvNS (x, *par):
    Enu = x[0]

    T       = par[0]
    knn     = par[1]
    sinsqqw = par[2]
    mN      = par[3]
    N       = par[4]
    Z       = par[5]

    #Form Factor 
    fofa = TF1 ("fofa", formfact, 0, 1e3, 2)
    fofa.SetParameter(0,mN)
    fofa.SetParameter(1, Z+N)
    Ff = fofa.Eval(T)

    val = knn*mN*(2 - 2*(T/Enu) + (T/Enu)**2 - (mN*T/Enu)) * (N*Z*(1-4*sinsqqw))**2 * Ff**2 * (1/4)

    return val
"""
#Function Double-Differential event rate in a TF2
def dSigdT_f2 (x, par):
    Enu = x [0]
    #T   = x[1]
    T  = par[0]
    """
    knn     = par[1]
    sinsqqw = par[2]
    mN      = par[3]
    N       = par[4]
    Z       = par[5]
    """    
    #Form Factor 
    fofa = TF1 ("fofa", formfact, 0.1, 1e3)
    #fofa.SetParameter(0,mN)
    #fofa.SetParameter(1, Z+N)
    Ff = fofa.Eval(T)
    #print(Ff)

    Ff = 1

    
    val = knn * M_Si * (hbarc**2) * (2 - 2*(T/Enu) + (T/Enu)**2 - (M_Si*T/Enu**2)) * ((N-Z*(1-4*sinsqqw))**2) * (Ff**2) * (1/4)
    
    return val


def ddRate(x, par):
    Enu = x[0]
    #T   = x[1]
    T   = par[0]
    """
    knn     = par[1]
    sinsqqw = par[2]
    mN      = par[3]
    N       = par[4]
    Z       = par[5]
    """

    loE = 0.5*( T + sqrt(pow(T,2)+2*T*M_Si))
    hiE = 10


    #---Cross Section--------

    dSigdT = TF1("dSigdT", dSigdT_f2,  loE, hiE, 1)
    dSigdT.SetParameter(0, T)
    """
    dSigdT.SetParameter(1, knn)
    dSigdT.SetParameter(2, sinsqqw)
    dSigdT.SetParameter(3, mN) 
    dSigdT.SetParameter(4, N)
    dSigdT.SetParameter(5, Z)
    """

    #----Flux-----------------    
    fluxE = TF1 ("fluxE", flux_gfm_v1, loE, hiE, 2)
    fluxE.SetParameter(1,0)


    ddr = 0.0   

    if Enu > loE:
        ddr = Nt(1000) * fluxE.Eval(Enu) * dSigdT.Eval(Enu)

    return ddr

#N = Nt(1000)
#print(N)


fddRate1 = TF1("fddRate1", ddRate, 0, 12, 1)
fddRate2 = TF1("fddRate2", ddRate, 0, 12, 1)
fddRate3 = TF1("fddRate3", ddRate, 0, 12, 1)
fddRate4 = TF1("fddRate4", ddRate, 0, 12, 1)
fddRate5 = TF1("fddRate5", ddRate, 0, 12, 1)
fddRate6 = TF1("fddRate6", ddRate, 0, 12, 1)
fddRate7 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate8 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate9 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate10 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate11 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate12 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate13 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate14 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate15 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate16 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate17 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate18 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate19 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate20 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate21 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate22 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate23 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate24 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate25 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate26 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate27 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate28 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate29 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate30 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate31 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate32 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate33 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate34 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate35 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate36 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate37 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate38 = TF1("fddRate7", ddRate, 0, 12, 1)
fddRate39 = TF1("fddRate7", ddRate, 0, 12, 1)


fddRate1.SetParameter(0, 0.05)
fddRate2.SetParameter(0, 0.1)
fddRate3.SetParameter(0, 0.15)
fddRate4.SetParameter(0, 0.20)
fddRate5.SetParameter(0, 0.25)
fddRate6.SetParameter(0, 0.30)
fddRate7.SetParameter(0, 0.35)
fddRate8.SetParameter(0, 0.40)
fddRate9.SetParameter(0, 0.45)
fddRate10.SetParameter(0, 0.50)
fddRate11.SetParameter(0, 0.55)
fddRate12.SetParameter(0, 0.60)
fddRate13.SetParameter(0, 0.65)
fddRate14.SetParameter(0, 0.70)
fddRate15.SetParameter(0, 0.75)
fddRate16.SetParameter(0, 0.80)
fddRate17.SetParameter(0, 0.85)
fddRate18.SetParameter(0, 0.90)
fddRate19.SetParameter(0, 0.95)
fddRate20.SetParameter(0, 1)
fddRate21.SetParameter(0, 1.5)
fddRate22.SetParameter(0, 2)
fddRate23.SetParameter(0, 2.5)
fddRate24.SetParameter(0, 3)
fddRate25.SetParameter(0, 3.5)
fddRate26.SetParameter(0, 4)
fddRate27.SetParameter(0, 4.5)
fddRate28.SetParameter(0, 5)
fddRate29.SetParameter(0, 5.5)
fddRate30.SetParameter(0, 6)
fddRate31.SetParameter(0, 6.5)
fddRate32.SetParameter(0, 7)
fddRate33.SetParameter(0, 7.5)
fddRate34.SetParameter(0, 8)
fddRate35.SetParameter(0, 8.5)
fddRate36.SetParameter(0, 9)
fddRate37.SetParameter(0, 9.5)
fddRate38.SetParameter(0, 10)
fddRate39.SetParameter(0, 10.5)



#------Style-----------
gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)
gStyle.SetTitleFontSize
#------Style------------
"""
canv = TCanvas("canv", "", 650, 450)
canv.cd()

gPad.SetLeftMargin(0.12)
gPad.SetRightMargin(0.08)
gPad.SetLogy(1)
gPad.SetTicks(1,1)

"""

#fddRate1.SetTitle("Double-diffrential Rate")
fddRate1.SetLineColor(1)
fddRate1.SetLineStyle(1)
fddRate1.SetLineWidth(1)

fddRate2.SetLineColor(1)
fddRate2.SetLineStyle(1)
fddRate2.SetLineWidth(1)

fddRate3.SetLineColor(1)
fddRate3.SetLineStyle(1)
fddRate3.SetLineWidth(1)

fddRate4.SetLineColor(1)
fddRate4.SetLineStyle(1)
fddRate4.SetLineWidth(1)

fddRate5.SetLineColor(1)
fddRate5.SetLineStyle(1)
fddRate5.SetLineWidth(1)

fddRate6.SetLineColor(1)
fddRate6.SetLineStyle(1)
fddRate6.SetLineWidth(1)

fddRate7.SetLineColor(1)
fddRate7.SetLineStyle(1)
fddRate7.SetLineWidth(1)

fddRate8.SetLineColor(1)
fddRate8.SetLineStyle(1)
fddRate8.SetLineWidth(1)

fddRate9.SetLineColor(1)
fddRate9.SetLineStyle(1)
fddRate9.SetLineWidth(1)

fddRate10.SetLineColor(1)
fddRate10.SetLineStyle(1)
fddRate10.SetLineWidth(1)

fddRate11.SetLineColor(1)
fddRate11.SetLineStyle(1)
fddRate11.SetLineWidth(1)

fddRate12.SetLineColor(1)
fddRate12.SetLineStyle(1)
fddRate12.SetLineWidth(1)

fddRate13.SetLineColor(1)
fddRate13.SetLineStyle(1)
fddRate13.SetLineWidth(1)

fddRate14.SetLineColor(1)
fddRate14.SetLineStyle(1)
fddRate14.SetLineWidth(1)

fddRate15.SetLineColor(1)
fddRate15.SetLineStyle(1)
fddRate15.SetLineWidth(1)

fddRate16.SetLineColor(1)
fddRate16.SetLineStyle(1)
fddRate16.SetLineWidth(1)

fddRate17.SetLineColor(1)
fddRate17.SetLineStyle(1)
fddRate17.SetLineWidth(1)

fddRate18.SetLineColor(1)
fddRate18.SetLineStyle(1)
fddRate18.SetLineWidth(1)

fddRate19.SetLineColor(1)
fddRate19.SetLineStyle(1)
fddRate19.SetLineWidth(1)

fddRate20.SetLineColor(1)
fddRate20.SetLineStyle(1)
fddRate20.SetLineWidth(1)

fddRate21.SetLineColor(1)
fddRate21.SetLineStyle(1)
fddRate21.SetLineWidth(1)

fddRate22.SetLineColor(1)
fddRate22.SetLineStyle(1)
fddRate22.SetLineWidth(1)

fddRate23.SetLineColor(1)
fddRate23.SetLineStyle(1)
fddRate23.SetLineWidth(1)

fddRate24.SetLineColor(1)
fddRate24.SetLineStyle(1)
fddRate24.SetLineWidth(1)

fddRate25.SetLineColor(1)
fddRate25.SetLineStyle(1)
fddRate25.SetLineWidth(1)

fddRate26.SetLineColor(1)
fddRate26.SetLineStyle(1)
fddRate26.SetLineWidth(1)

fddRate27.SetLineColor(1)
fddRate27.SetLineStyle(1)
fddRate27.SetLineWidth(1)

fddRate28.SetLineColor(1)
fddRate28.SetLineStyle(1)
fddRate28.SetLineWidth(1)

fddRate29.SetLineColor(1)
fddRate29.SetLineStyle(1)
fddRate29.SetLineWidth(1)

fddRate30.SetLineColor(1)
fddRate30.SetLineStyle(1)
fddRate30.SetLineWidth(1)

fddRate31.SetLineColor(1)
fddRate31.SetLineStyle(1)
fddRate31.SetLineWidth(1)

fddRate32.SetLineColor(1)
fddRate32.SetLineStyle(1)
fddRate32.SetLineWidth(1)

fddRate33.SetLineColor(1)
fddRate33.SetLineStyle(1)
fddRate33.SetLineWidth(1)

fddRate34.SetLineColor(1)
fddRate34.SetLineStyle(1)
fddRate34.SetLineWidth(1)

fddRate35.SetLineColor(1)
fddRate35.SetLineStyle(1)
fddRate35.SetLineWidth(1)

fddRate36.SetLineColor(1)
fddRate36.SetLineStyle(1)
fddRate36.SetLineWidth(1)

fddRate37.SetLineColor(2)
fddRate37.SetLineStyle(1)
fddRate37.SetLineWidth(1)

fddRate38.SetLineColor(4)
fddRate38.SetLineStyle(1)
fddRate38.SetLineWidth(1)

fddRate39.SetLineColor(8)
fddRate39.SetLineStyle(1)
fddRate39.SetLineWidth(1)


#ddRate_frame = TH2F("ddRate_frame","",1,-0.2,9.2,1,2e-4,10)
fddRate1.SetTitle("Double-Differential event rate")
fddRate1.GetXaxis().SetTitle("E_{#bar{#nu}_{e}} [MeV]")
fddRate1.GetYaxis().SetTitle("d^{2}R/dE_{#bar{#nu}_{e}}dT  [ #bar{#nu}_{e} cm^{2} / MeV^{2} fis]")
fddRate1.GetYaxis().SetTitleOffset(1.5)

leg = TLegend(0.71,0.71,0.87,0.85)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.03)
leg.SetHeader("T_{3} < T_{2} < T_{1}","C")
leg.AddEntry(fddRate37, "T_{3}")
leg.AddEntry(fddRate38, "T_{2}")
leg.AddEntry(fddRate39, "T_{1}")

canv = TCanvas("canv")


gPad.SetLeftMargin(0.13)
gPad.SetRightMargin(0.07)
gPad.SetLogy(1)
gPad.SetTicks(1,1)


fddRate1.Draw()
fddRate2.Draw("same")
fddRate3.Draw("same")
fddRate4.Draw("same")
fddRate5.Draw("same")
fddRate6.Draw("same")
fddRate7.Draw("same")
fddRate8.Draw("same")
fddRate9.Draw("same")
fddRate10.Draw("same")
fddRate11.Draw("same")
fddRate12.Draw("same")
fddRate13.Draw("same")
fddRate14.Draw("same")
fddRate15.Draw("same")
fddRate16.Draw("same")
fddRate17.Draw("same")
fddRate18.Draw("same")
fddRate19.Draw("same")
fddRate20.Draw("same")
fddRate21.Draw("same")
fddRate22.Draw("same")
fddRate23.Draw("same")
fddRate24.Draw("same")
fddRate25.Draw("same")
fddRate26.Draw("same")
fddRate27.Draw("same")
fddRate28.Draw("same")
fddRate29.Draw("same")
fddRate30.Draw("same")
fddRate31.Draw("same")
fddRate32.Draw("same")
fddRate33.Draw("same")
fddRate34.Draw("same")
fddRate35.Draw("same")
fddRate36.Draw("same")
fddRate37.Draw("same")
fddRate38.Draw("same")
fddRate39.Draw("same")
leg.Draw()


#canv.cd()

print(knn)

#fddRate.Draw()

canv.Print("ddRate_v1.pdf")

input()
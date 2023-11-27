def flux_gfm_v1( x, *par):
    
    from ROOT import TCanvas,TH1F,TFile,TGraph, TF1
    from ROOT import gROOT, gBenchmark, gRandom, gSystem
    from array import array
    
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

    return ff * normfact



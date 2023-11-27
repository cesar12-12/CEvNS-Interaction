import numpy
from ROOT import TCanvas,TH1F,TFile,TGraph, TF1, TMath, TF2
from ROOT import gROOT, gBenchmark, gRandom, gSystem, gStyle
from array import array
from math import sqrt, exp, sin, cos


#Some physical constants

c_pi        =3.141592653589793238
DeltaMnp    = 1.2933 #MeV
mel         = 0.511 #MeV
N_avo       = 6.022e23
MeV_per_amu = 931.494061
G_F         = 1.166378e-11 #(Mev^-2)
hbarc       = 197.326e-13   #Mev cm
sinsqqw     = 0.231 #(sin qw)^2 approx
Si_amu      = 28.085  # g / mol
M_Si        = 2.910525e-13 #MeV Mass of the Si in MeV
knn = G_F**2 / (2*c_pi)
Z = Si_amu -14
N = 14

    
Er    = 1
kappa = sqrt(2*M_Si*Er) / (hbarc*1e13)
s     = 1.0
R     = 1.2*pow(Si_amu,1.0/3.0)
#print(R)
r     = sqrt(R**2 - 5*(s**2))
ff    = 3*exp(- 0.5 * kappa**2 * s**2) * (sin(kappa*r) - kappa*r*cos(kappa*r)) / (kappa*r)**3

T = 0.5
loE = 0.5*( T + sqrt(pow(T,2)+2*T*M_Si))
print(T, loE)

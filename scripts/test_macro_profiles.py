# configure (get file name etc.)
import calo_init
calo_init.add_defaults()
calo_init.parser.add_argument("sf", help="SF", type = float)
calo_init.parse_args()
calo_init.print_config()
SF = calo_init.args.sf
#Setup ROOT
from ROOT import gSystem
gSystem.Load("libCaloAnalysis")
from ROOT import ShowerProfiles, gStyle, TCanvas, TFile, TF1, gPad, TMath
#import draw functions
from draw_functions import draw_1histogram, draw_2histograms

# use this script for only 1 file
energy = calo_init.energies[0]
filename = calo_init.filenamesIn[0]
if len(calo_init.filenamesIn) > 1:
    print "WARNING: analysis of the first input file, ignoring the rest"


ma = ShowerProfiles(energy, SF)
ma.loop(filename, calo_init.verbose)
print "Mean hit energy: ", ma.h_hitEnergy.GetMean()
print "1/SF calculated: ", energy/(ma.h_hitEnergy.GetMean())

  #Longo-Sestili formula
  #http://arxiv.org/pdf/hep-ex/0001020.pdf
fit = TF1("fit","[0]*(pow([2]*x,[1]-1)*[2]*exp(-[2]*(x))/TMath::Gamma([1]))")
fit.SetParName(0,"A") #normalization factor
fit.SetParName(1,"#alpha")
fit.SetParName(2,"#beta")

gStyle.SetOptStat("emr");

c1 = TCanvas("c1","c1",1000,1000)
c1.Divide(3,2)
c1.cd(1)
draw_1histogram(ma.h_ptGen, "p_{T}^{gen} [GeV]","")
c1.cd(2)
draw_1histogram(ma.h_pdgGen, "PDG code","")
c1.cd(4)
ma.h_cellEnergy.Rebin(2)
draw_1histogram(ma.h_cellEnergy, "Total cell energy [GeV]","")
if (ma.h_cellEnergy.GetEntries()>50):
    ma.h_cellEnergy.Fit("gaus")
c1.cd(5)
draw_2histograms(ma.h_longProfile_particle, ma.h_longProfile, "Longitudinal distance [mm]", "Energy [GeV]", "Particle dir.", "Hits in 1st layer")
fit.SetParameters(100, 8.15/1.15,1.0/1.15);
ma.h_longProfile.Fit("fit")
gPad.SetLogy(1)
c1.cd(6)
draw_2histograms(ma.h_radialProfile_particle, ma.h_radialProfile, "Radial distance [mm]", "Energy [GeV]", "Particle dir.", "Hits in 1st layer")

raw_input("Press ENTER to exit")
#c1.SaveAs("plots_"+PARTICLE+str(ENERGY)+".gif")

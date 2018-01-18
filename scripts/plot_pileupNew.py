import calo_init
## add arguments relevant only for that script
calo_init.add_defaults()
calo_init.parser.add_argument("--histoname", help="Name of the pile-up histograms", type = str)
calo_init.parser.add_argument("--nLayers", help="Number of layers", type = int)
calo_init.parser.add_argument("--nBins", help="Number of |eta|-bins", type = int)
calo_init.parser.add_argument("--mu", help="Number of pileup collisions to scale data to", type = int)
calo_init.parse_args()

from math import sqrt
from ROOT import gSystem
gSystem.Load("libCaloAnalysis")
from ROOT import TCanvas, TFile, gStyle, gPad, kGreen, kRed, kBlue, TColor, TPad, TH1, TH2
from draw_functions import *

# use this script for multiple files
# gStyle.SetPalette(56) # kInvertedDarkBodyRadiator
gStyle.SetPalette(73) # kCMYK
gStyle.SetOptFit(1)
gStyle.SetOptStat(0)

nlayers = 8
nbins = 20
mu = 1000
histoname = "energyVsAbsEta"

if calo_init.args.histoname:
    histoname = calo_init.args.histoname
if calo_init.args.nLayers:
    nlayers = calo_init.args.nLayers
if calo_init.args.nBins:
    nbins = calo_init.args.nBins
if calo_init.args.mu:
    mu = calo_init.args.mu

print "Histogram name: ", histoname

# 1.6: correction for out-of-time pile-up
#scale = 1.6*sqrt(mu)
# only in-time pile-up
scale = sqrt(mu)
print "Scaling factor: ", scale

hEnVsAbsEta = []
hPileup = []

for ifile, filename in enumerate(calo_init.filenamesIn):
    energy = calo_init.energy(ifile)
    print "Files with pileup: " + filename
    f = ROOT.TFile(filename,"r")
    # get the TH2 histograms, make a projection to extract the RMS per eta bin
    for index in xrange(nlayers):
        histoname_i = histoname+str(index)
        hEnVsAbsEta.append(f.Get(histoname_i))
        if (hEnVsAbsEta[index].GetNbinsX()<nbins):
            nbins = hEnVsAbsEta[index].GetNbinsX()
            print "More bins than available required!!! Setting nbins to ", nbins
        etaMax = nbins*hEnVsAbsEta[index].GetXaxis().GetBinWidth(1)
        hPile = TH1F("h_pileup_layer"+str(index+1),"Pileup in "+str(index+1)+". layer", nbins, 0, etaMax)
        hPile.GetXaxis().SetTitle("|#eta|")
        hPile.GetYaxis().SetTitle("Pileup [GeV]")
        if (index!=0):
            hPile.SetMaximum(0.2)
        for i in range(1, nbins):
            hProj = hEnVsAbsEta[index].ProjectionY( "hprojection_"+str(index), i, i,"e")
            hPile.SetBinContent(i,hProj.GetRMS()*scale)
        hPileup.append(hPile)

    #Draw pileup per layer
    canv = TCanvas('Minimum bias event in ECal', 'ECal', 1600, 1000 )
    canv.Divide(3,3)
    for i in xrange(nlayers):
        canv.cd(i+1)
        hPileup[i].Draw()

    if calo_init.output(ifile):
        canv.Print(calo_init.output(ifile)+".gif")
        plots = TFile(calo_init.output(ifile)+".root","RECREATE")
    else:
        canv.Print("TestPileup_mu"+str(mu)+".gif")
        plots = TFile("TestPileup_mu"+str(mu)+".root","RECREATE")
        
    plots.cd()
    for hset in [hPileup]:
        for h in hset:
            h.Write()
    plots.Close()
    
#raw_input("Press ENTER to exit")

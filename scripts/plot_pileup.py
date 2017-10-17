import calo_init
## add arguments relevant only for that script
calo_init.add_defaults()
calo_init.parser.add_argument("--nLayers", help="Number of layers", type = int)
calo_init.parser.add_argument("--nBins", help="Number of |eta|-bins", type = int)
calo_init.parser.add_argument("--mu", help="Number of pileup collisions to scale data to", type = int)
calo_init.parse_args()

from math import sqrt
from ROOT import gSystem
gSystem.Load("libCaloAnalysis")
from ROOT import TCanvas, TFile, gStyle, gPad, kGreen, kRed, kBlue, TColor, TPad, TH1, TH2, TLegend
from draw_functions import *

# use this script for multiple files
# gStyle.SetPalette(56) # kInvertedDarkBodyRadiator
gStyle.SetPalette(73) # kCMYK
gStyle.SetOptFit(1)
gStyle.SetOptStat(0)

nlayers = 8
nbins = 20
mu = 1000

if calo_init.args.nLayers:
    nlayers = calo_init.args.nLayers
if calo_init.args.nBins:
    nbins = calo_init.args.nBins
if calo_init.args.mu:
    mu = calo_init.args.mu

# 1.6: correction for out-of-time pile-up
#scale = 1.6*sqrt(mu)
scale = 1.0*sqrt(mu)

print "Scaling factor: ", scale

hEnVsAbsEta = []
hProjectionFirstEta = []
hPileup = []

for ifile, filename in enumerate(calo_init.filenamesIn):
    energy = calo_init.energy(ifile)
    print "Files with pileup: " + filename
    f = ROOT.TFile(filename,"r")
    # retrieve histograms to draw them
    hEnCell = f.Get("cellEnergy")
    hEnCellTest = f.Get("cellEnergyTest")
    # get the TH2 histograms, make a projection to extract the RMS per eta bin
    for index in xrange(nlayers):
        histoname = "EnergyVsAbsEta_"+str(index)
        hEnVsAbsEta.append(f.Get(histoname))
        if (hEnVsAbsEta[index].GetNbinsX()<nbins):
            nbins = hEnVsAbsEta[index].GetNbinsX()
            print "More bins than available required!!! Setting nbins to ", nbins
        etaMax = nbins*hEnVsAbsEta[index].GetXaxis().GetBinWidth(1)
        #print "Bin width ", hEnVsAbsEta[index].GetXaxis().GetBinWidth(1)
        hProjectionFirstEta.append(hEnVsAbsEta[index].ProjectionY( "hprojection_"+str(index), 1, 1,"e"))
        hPile = TH1F("h_pileup_layer"+str(index+1),"", nbins, 0, etaMax)
        hPile.GetXaxis().SetTitle("|#eta|")
        hPile.GetYaxis().SetTitle("Pileup [GeV]")
        hPile.SetLineColor(index+1)
        hPile.SetLineWidth(2)
        if (index!=0):
            hPile.SetMaximum(0.2)
        for i in range(1, nbins):
            hProj = hEnVsAbsEta[index].ProjectionY( "hprojection_"+str(index), i, i,"e")
            if index==0:
                print "layer ", index, "bin ", i, " rms (1 event) ", hProj.GetRMS() 
            hPile.SetBinContent(i,hProj.GetRMS()*scale)
        hPileup.append(hPile)

#Draw pileup per layer
canv = prepare_single_canvas("pileup","Pile-up")
canv.cd()
legend = TLegend(0.2,0.6,0.4,0.85)
legend.SetBorderSize(0)

# Prepare graphs: set colours, axis range, ...
colour = [c for c in range(1,100)]
print(colour)
for i,g in enumerate(hPileup):
    prepare_graph(g, g.GetName(), g.GetTitle(), colour[i])
    if (i==0):
        g.Draw()
    else:
        g.Draw("same")
    g.GetYaxis().SetRangeUser(0,0.31)

    legend.AddEntry(g,"layer " + str(i+1),"l")

#legend.Draw()
draw_text(["FCC-hh simulation"], [0.67,0.88, 0.95,0.98], kGray+3, 0).SetTextSize(0.05)
draw_text(["Pile-up noise per cell"], [0.25,0.88, 0.45,0.98], 1, 0).SetTextSize(0.05)
draw_text(["in-time pile-up, #LT#mu#GT = "+str(mu)], [0.33,0.78, 0.45,0.88], 1, 0).SetTextSize(0.05)
canv.Update()

print "cellEnergyTest", hEnCellTest.GetEntries(), hEnCellTest.GetMean(), hEnCellTest.GetRMS(), hEnCellTest.Integral(0,161)
#print "projection", hProjectionFirstEta[0].GetEntries(), hProjectionFirstEta[0].GetMean(), hProjectionFirstEta[0].GetRMS()

if calo_init.output(ifile):
    canv.Print(calo_init.output(ifile)+".gif")
    plots = TFile(calo_init.output(ifile)+".root","RECREATE")
else:
    canv.Print("FinalPileup_mu"+str(mu)+".gif")
    plots = TFile("FinalPileup_mu"+str(mu)+".root","RECREATE")
    
plots.cd()
for hset in [hPileup]:
    for h in hset:
        h.Write()
plots.Close()
    
#raw_input("Press ENTER to exit")

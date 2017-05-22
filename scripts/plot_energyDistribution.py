import calo_init
calo_init.add_defaults()
calo_init.parser.add_argument("-n","--histogramName", default=["total"], help="Name(s) of the energy distribution histogram in input file(s)", type = str, nargs='+')
calo_init.parser.add_argument("--legend","-l",default=[],type=str,nargs='+')
calo_init.parser.add_argument("-m","--axisMax", help="Maximum of the axis", type = float)
calo_init.parser.add_argument("--roundBrackets", help="Use round brackets for unit", action = 'store_true')
calo_init.parser.add_argument("--preview", help="Plot preview of fits", action = 'store_true')
calo_init.parser.add_argument("--specialLabel", help="Additional label to be plotted", type=str, default = "FCC-hh simulation")
calo_init.parse_args()
calo_init.print_config()

if calo_init.args.roundBrackets:
    energyUnit = '(GeV)'
else:
    energyUnit = '[GeV]'
if len(calo_init.args.histogramName) == 1:
    histogramNames = [calo_init.args.histogramName[0] for i in range (0, len(calo_init.filenamesIn))]
else:
    histogramNames = calo_init.args.histogramName
if len(histogramNames) > 1 and len(histogramNames) != len(calo_init.filenamesIn):
    print("If names of the histograms differ in each input file, the same amount of filenames and histogram names needs to be given.")
    exit()

if len(calo_init.args.legend) == 0:
    legend = ["" for i in range (0, len(calo_init.filenamesIn))]
elif len(calo_init.args.legend) == 1:
    legend, check = calo_init.substitute(calo_init.args.legend[0])
    legend = [leg + ": " for leg in legend]
else:
    legend = calo_init.args.legend
    legend = [leg + ": " for leg in legend]
print (legend)

from ROOT import gSystem, gROOT, TCanvas, TGraphErrors, TF1, gStyle, kRed, kBlue, kGray, TFile, TTree, TPad, TGaxis, gPad, TLine, kOrange, kGreen, SetOwnership
from draw_functions import prepare_histogram, prepare_divided_canvas,  prepare_single_canvas, draw_text, draw_1histogram
import numpy
gStyle.SetOptFit(0)
gStyle.SetOptStat(0)

colours = [ kRed + 2, 8, kGray + 3 , kOrange + 1 , 8, 9]
coloursFit = [kRed, kGreen+2, kGray+3, kOrange + 2 , kGreen + 2, kBlue + 2]

if len(calo_init.energies) > 1:
    print("Plotting of the enrgy distribution can be performed only for one energy value")
    exit()

energy = calo_init.energy(0)
canv = prepare_single_canvas( 'energyDistribution_'+str(energy)+'GeV', 'Energy distribution for '+str(energy)+'GeV' )
canv.cd()
energy = calo_init.energy(0)
resolution_list = []
mean_list = []
# files represent different types for the same particle energy:
for ifile, filename in enumerate(calo_init.filenamesIn):
    f = TFile(filename, "READ")
    hEn = f.Get(histogramNames[ifile])
    hEn.Sumw2()
    if hEn.GetEntries() > 0:
        hEn.Scale(1./hEn.GetEntries())
    hEn.Rebin(4)
    hEn.SetTitle('')
    if calo_init.args.axisMax:
        hEn.GetYaxis().SetRangeUser(0, calo_init.args.axisMax)
    fitPre = TF1("fitPre","gaus", hEn.GetMean() - 1. * hEn.GetRMS(), hEn.GetMean() + 1. * hEn.GetRMS())
    resultPre = hEn.Fit(fitPre, "SQRN")
    fit = TF1("fit","gaus", resultPre.Get().Parameter(1) - 2. * resultPre.Get().Parameter(2), resultPre.Get().Parameter(1) + 2. * resultPre.Get().Parameter(2))
    fit.SetLineColor(coloursFit[ifile])
    result = hEn.Fit(fit, "SQRN")
    mean = result.Get().Parameter(1)
    sigma = result.Get().Parameter(2)
    dMean = result.Get().Error(1)
    dSigma = result.Get().Error(2)
    resolution_list.append("#color["+str(colours[ifile])+"]{"+legend[ifile]+str(round(sigma / mean * 100 ,2))+" %}")
    mean_list.append("#color["+str(colours[ifile])+"]{"+legend[ifile]+str(round((mean-energy)/energy*100 ,1))+" %}")
    prepare_histogram(hEn, ';E_{beam} '+energyUnit+';fraction of events', colours[ifile])
    hEn.GetXaxis().SetRangeUser(0.9 * energy, 1.1* energy)
    if ifile == 0:
        hEn.DrawCopy()
    else:
        hEn.DrawCopy("sameep")
    fit.Draw("same")

canv.Update()

draw_text(["energy resolution", "#color[1]{#sigma_{E_{rec}}/#LTE_{rec}#GT}"]
          + resolution_list,
          [0.17,0.8-0.05*len(resolution_list),0.42,0.88], 1, 1)
draw_text(["linearity", "#color[1]{(#LTE_{rec}#GT-E_{beam})/E_{beam}}"]
          + mean_list,
          [0.7,0.8-0.05*len(resolution_list),0.95,0.88], 1, 1)
if calo_init.args.specialLabel:
    draw_text([calo_init.args.specialLabel], [0.57,0.88, 0.85,0.98], kGray+3, 0).SetTextSize(0.05)
canv.Update()

if calo_init.output(ifile):
    canv.Print("energyDistribution_"+calo_init.output(0)+".png")
else:
    canv.Print("energyDistribution_"+str(energy)+"GeV.png")


input("Press ENTER to exit")

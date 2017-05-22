import calo_init
calo_init.add_defaults()
calo_init.parser.add_argument("--merge", help="merge layers", default = [1] * 32, type = int, nargs='+') # bin 0 is empty! (before calo)
calo_init.parser.add_argument("-t","--title", default="Sampling fraction", help="Graph title", type=str)
calo_init.parser.add_argument("-n","--histogramName", default="ecal_sf_layer", help="Name of the histogram with sampling fraction (postfixed with number of layer)", type = str)
calo_init.parser.add_argument("--histogramNameMean", default="ecal_sf", help="Name of the histogram with sampling fraction (sufixed with number of layer)", type = str)
calo_init.parser.add_argument("-max","--axisMax", help="Maximum of the axis", type = float)
calo_init.parser.add_argument("-min","--axisMin", help="Minimum of the axis", type = float)
calo_init.parser.add_argument("--totalNumLayers", default = 32, help="Total number of the layers used in simulation", type = int)
calo_init.parser.add_argument("--numFirstLayer", default = 1, help="ID of first layer used in histograms name", type = int)
calo_init.parser.add_argument("--layerWidth", default = 2, help="Width of the layer (cm)", type = float)
calo_init.parser.add_argument("--X0density", default = 0.422, help="Xo density of a current detector (X0/cm)", type = float)
calo_init.parser.add_argument("--roundBrackets", help="Use round brackets for unit", action = 'store_true')
calo_init.parser.add_argument("--preview", help="Plot preview of fits", action = 'store_true')
calo_init.parser.add_argument("--specialLabel", help="Additional label to be plotted", type=str, default = "FCC-hh simulation")
calo_init.parse_args()
calo_init.print_config()

histName = calo_init.args.histogramName
histNameMean = calo_init.args.histogramNameMean

from ROOT import gSystem, gROOT, TCanvas, TGraphErrors, TF1, gStyle, kRed, kBlue, kGray, TFile, TTree, TPad, TGaxis, gPad, TLine
from draw_functions import prepare_graph, prepare_divided_canvas,  prepare_single_canvas, draw_text, draw_1histogram
import numpy
from math import sqrt, ceil, floor

merge = [sum(calo_init.args.merge[:i]) for i in range(0,len(calo_init.args.merge))]
sliceWidth = calo_init.args.layerWidth  # cm
startIndex = calo_init.args.numFirstLayer
Nslices = calo_init.args.totalNumLayers
if sum(calo_init.args.merge) != Nslices:
    print('Number of total layers (',Nslices,') is not the same as a sum of "--merge" arguments (sum = ',sum(calo_init.args.merge),')')
    exit(0)
Nslicesmerged = len(merge)
all_graphs = []
avgSF = []
avgSFerr = []

# first get all the resolutions and prepare graphs
for ifile, filename in enumerate(calo_init.filenamesIn):
    energy = calo_init.energy(ifile)
    f = TFile(filename, "READ")
    # mean value of the sampling fraction
    hMean = f.Get(histNameMean)
    if hMean:
        fitPre = TF1("fitPre","gaus", hMean.GetMean() - 1. * hMean.GetRMS(), hMean.GetMean() + 1. * hMean.GetRMS())
        resultPre = hMean.Fit(fitPre, "SQRN")
        fit = TF1("fit","gaus",resultPre.Get().Parameter(1) - 2. * resultPre.Get().Parameter(2), resultPre.Get().Parameter(1) + 2. * resultPre.Get().Parameter(2) )
        result = hMean.Fit(fit, "SQRN")
        avgSF.append(result.Get().Parameter(1))
        avgSFerr.append(result.Get().Parameter(2))
    hmerged = []
    # first merge adjacent layers and get histograms of SF
    for islice in range(startIndex, Nslices + startIndex):
        h = f.Get(histName+str(islice))
        # if first hist to be merged
        lastIm = -1
        if islice - startIndex in merge:
            lastIm += 1
            hmerged.append(h)
        else:
            hmerged[lastIm].Add(h)
    gSF = TGraphErrors()
    # now fit SF with Gaussians
    if calo_init.args.preview:
        cPreview = prepare_divided_canvas('preview_e'+str(energy)+'GeV', 'Preview for '+str(energy)+'GeV', Nslicesmerged)
        fitoptions = "SQR"
    else:
        fitoptions = "SQRN"
    for islice, h in enumerate(hmerged):
        fitPre = TF1("fitPre","gaus", h.GetMean() - 1. * h.GetRMS(), h.GetMean() + 1. * h.GetRMS())
        resultPre = h.Fit(fitPre, fitoptions)
        fit = TF1("fit","gaus",resultPre.Get().Parameter(1) - 2. * resultPre.Get().Parameter(2), resultPre.Get().Parameter(1) + 2. * resultPre.Get().Parameter(2) )
        result = h.Fit(fit, fitoptions)
        if result and result.Ndf() > 0:
            # if it fits terribly, try to fit in narrower range
            if result.Chi2() / result.Ndf() > 10:
                refit = TF1("refit","gaus",resultPre.Get().Parameter(1) - resultPre.Get().Parameter(2), resultPre.Get().Parameter(1) + resultPre.Get().Parameter(2) )
                result = h.Fit(refit, fitoptions)
        if calo_init.args.preview:
            cPreview.cd(islice + 1)
            draw_1histogram(h,"","")
        # make graph
        if result:
            if islice < len(merge) - 1:
                gSF.SetPoint(islice, (merge[islice] + 0.5 * (merge[islice + 1] - merge[islice])) * sliceWidth, result.Get().Parameter(1))
                gSF.SetPointError(islice, 0.5 * (merge[islice + 1] - merge[islice]) * sliceWidth , result.Get().Parameter(2))
            else:
                gSF.SetPoint(islice, (merge[islice] + 0.5 * (merge[islice] - merge[islice - 1])) * sliceWidth, result.Get().Parameter(1))
                gSF.SetPointError(islice, 0.5 * (merge[islice] - merge[islice - 1]) * sliceWidth , result.Get().Parameter(2))
    all_graphs.append(gSF)

canv = prepare_single_canvas('sf_e'+str(energy)+'GeV', 'Sampling fraction for '+str(energy)+'GeV')

# Draw graph and all labels
prepare_graph(gSF, 'sf_'+str(len(merge))+'layers', ';radial depth;sampling fraction', ifile+9)
all_graphs[0].Draw("ape")
for g in all_graphs[1:]:
    g.Draw("pe")
if calo_init.args.axisMax:
    all_graphs[0].SetMaximum(calo_init.args.axisMax)
if calo_init.args.axisMin:
    all_graphs[0].SetMinimum(calo_init.args.axisMin)
canv.Update()

lines = []
for iLine, line in enumerate(avgSF):
    lines.append(TLine(0, avgSF[iLine], Nslices * sliceWidth, avgSF[iLine]))
    lines[iLine].SetLineColor(iLine+9)
    lines[iLine].Draw('same')

# add second axis
canv.SetRightMargin(0.1)
all_graphs[0].GetXaxis().SetRangeUser(0,68)
all_graphs[0].GetXaxis().SetLabelOffset(0.02)
gPad.RangeAxis(0,gPad.GetUymin(),68,gPad.GetUxmax())
axis = TGaxis(gPad.GetUxmin(),
              gPad.GetUymin(),
              gPad.GetUxmax(),
              gPad.GetUymin(),
              gPad.GetUxmin() * calo_init.args.X0density,
              gPad.GetUxmax() * calo_init.args.X0density,506,"-")
axis.SetLabelSize(0.05)
axis.SetTitleSize(0.07)
axis.SetTitleOffset(0.9)
axis.SetLabelOffset(0.02)
axis.Draw()
axis.SetLabelFont(42)
unit1 = draw_text(['(cm)'],[0.91,0.09,1,0.14] , 1, 0)
unit1.SetTextSize(0.05)
unit1.SetTextFont(42)
unit2 = draw_text(['(X_{0})'],[0.91,0.16,1,0.21] , 1, 0)
unit2.SetTextSize(0.05)
unit2.SetTextFont(42)
canv.Update()


# Draw all labels
if calo_init.args.specialLabel:
    draw_text([calo_init.args.specialLabel], [0.57,0.88, 0.85,0.98], kGray+3, 0).SetTextSize(0.05)
canv.Update()

# Save canvas and root file with graph, const term and sampling term
if calo_init.output(0):
    canv.SaveAs(calo_init.output(0)+".pdf")
    canv.SaveAs(calo_init.output(0)+".png")
    plots = TFile(calo_init.output(0)+".root","RECREATE")
    if calo_init.args.preview:
        cPreview.SaveAs("preview_"+calo_init.output(0)+".png")
else:
    canv.SaveAs("sampling_fraction_plots.pdf")
    canv.SaveAs("sampling_fraction_plots.png")
    plots = TFile("sampling_fraction.root","RECREATE")
    if calo_init.args.preview:
        cPreview.SaveAs("preview_sampling_fraction.png")
for g in all_graphs:
    g.Write()

mean = numpy.zeros(1, dtype=float)
std = numpy.zeros(1, dtype=float)
t = TTree("samplingFraction", "Sampling fraction for detector layers")
t.Branch("mean", mean, "mean/D");
t.Branch("std", std, "std/D");
for islice in range(0, Nslicesmerged):
    for ilay in range(0, calo_init.args.merge[islice]):
        mean[0] = gSF.GetY()[islice]
        std[0] = gSF.GetErrorY(islice)
        t.Fill()
plots.Write()
plots.Close()

print("============================================================")
print("== to be used in FCCSW, with CalibrateInLayers algorithm: ==")
print("============================================================")
print("samplingFraction = ",end="")
for islice in range(0, Nslicesmerged):
    if islice > 0:
        print(" + ",end="")
    print("["+str(gSF.GetY()[islice])+"] * "+str(calo_init.args.merge[islice]),end="")
print()
print("============================================================")

input("Press ENTER to exit")

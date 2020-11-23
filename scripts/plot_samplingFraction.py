from __future__ import print_function
import calo_init
calo_init.add_defaults()
calo_init.parser.add_argument("--merge", help="merge layers", default = [1] * 8, type = int, nargs='+') # bin 0 is empty! (before calo)
calo_init.parser.add_argument("-t","--title", default="Sampling fraction", help="Graph title", type=str)
calo_init.parser.add_argument("-n","--histogramName", default="ecal_sf_layer", help="Name of the histogram with sampling fraction (postfixed with number of layer)", type = str)
calo_init.parser.add_argument("--histogramNameMean", default="ecal_sf", help="Name of the histogram with sampling fraction (sufixed with number of layer)", type = str)
calo_init.parser.add_argument("-max","--axisMax", help="Maximum of the axis", type = float)
calo_init.parser.add_argument("-min","--axisMin", help="Minimum of the axis", type = float)
calo_init.parser.add_argument("--totalNumLayers", default = 8, help="Total number of the layers used in simulation", type = int)
calo_init.parser.add_argument("--numFirstLayer", default = 0, help="ID of first layer used in histograms name", type = int)
calo_init.parser.add_argument("--layerWidth", default = [2] , nargs='+', help="Width of the layers (cm). One value for identical widths of each layer.", type = float)
calo_init.parser.add_argument("--X0density", default = 0.422, help="Xo density of a current detector (X0/cm)", type = float)
calo_init.parser.add_argument("--roundBrackets", help="Use round brackets for unit", action = 'store_true')
calo_init.parser.add_argument("--preview", help="Plot preview of fits", action = 'store_true')
calo_init.parser.add_argument("--specialLabel", help="Additional label to be plotted", type=str, default = "")
calo_init.parse_args()
print('START')
calo_init.print_config()

histName = calo_init.args.histogramName
histNameMean = calo_init.args.histogramNameMean

from ROOT import gSystem, gROOT, TCanvas, TH1F, TGraphErrors, TF1, gStyle, kRed, kBlue, kGray, TFile, TTree, TPad, TGaxis, gPad, TLine, TColor
from draw_functions import prepare_graph, prepare_divided_canvas,  prepare_single_canvas, draw_text, draw_1histogram
import numpy
from math import sqrt, ceil, floor

merge = [sum(calo_init.args.merge[:i]) for i in range(0,len(calo_init.args.merge))]
print('merge',merge)
sliceWidth = calo_init.args.layerWidth  # cm
if len(sliceWidth) == 1:
    sliceWidth = sliceWidth * len(merge)
sliceSum = []
sumWidths = 0
for width in sliceWidth:
    sumWidths += width
    sliceSum.append(sumWidths)
print('sliceWidths',sliceWidth)
startIndex = calo_init.args.numFirstLayer
Nslices = calo_init.args.totalNumLayers
if sum(calo_init.args.merge) != Nslices:
    print('Number of total layers (',Nslices,') is not the same as a sum of "--merge" arguments (sum = ',sum(calo_init.args.merge),')')
    exit(0)
Nslicesmerged = len(merge)
all_graphs = []
graphTitles = []
avgSF = []
avgSFerr = []

colour = ['#4169E1','#D2691E','#228B22','#DC143C','#696969','#9932CC','#D2B48C', 1, 2, 3, 4, 5, 6, 7, 8, 9,10]
colour = [TColor.GetColor(c) for c in colour]

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
        h = TH1F() 
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
        h.Rebin(1)
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
            gSF.SetPoint(islice, sliceSum[islice]-sliceWidth[islice]*0.5, result.Get().Parameter(1))
            gSF.SetPointError(islice, sliceWidth[islice]*0.5 , result.Get().Parameter(2))
    prepare_graph(gSF, 'sf_'+str(len(merge))+'layers', ';radial depth [cm];sampling fraction', ifile+9)
    all_graphs.append(gSF)
    graphTitles.append('#color['+str(colour[ifile])+']{'+str(energy)+' GeV e^{-}}')
    print("samplFractMap ["+str(energy)+"]= ", end='')
    for islice in range(0, Nslicesmerged):
        if islice == 0:
            print("{ ", end='')
        if islice > 0:
            print(", ", end='')
        print(str(gSF.GetY()[islice]), end='')
    print("};")


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
    lines.append(TLine(0, avgSF[iLine], 65, avgSF[iLine]))
    lines[iLine].SetLineColor(colour[iLine])
    all_graphs[iLine].SetMarkerColor(colour[iLine])
    all_graphs[iLine].SetLineColor(colour[iLine])
    lines[iLine].Draw('same')

if len(graphTitles) > 1:
    draw_text(graphTitles, [0.18,0.9 - 0.07 * len(graphTitles),0.4,0.95], 0.4, 0).SetTextSize(0.06)

# Draw all labels
if calo_init.args.specialLabel:
    draw_text([calo_init.args.specialLabel], [0.57,0.88, 0.85,0.98], kGray+3, 0).SetTextSize(0.05)
canv.Update()

# Save canvas and root file with graph, const term and sampling term
if calo_init.output(0):
    canv.SaveAs(calo_init.output(0)+".eps")
    canv.SaveAs(calo_init.output(0)+".png")
    plots = TFile(calo_init.output(0)+".root","RECREATE")
    if calo_init.args.preview:
        cPreview.SaveAs("preview_"+calo_init.output(0)+".png")
else:
    canv.SaveAs("sampling_fraction_plots.eps")
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

raw_input("Press ENTER to exit")

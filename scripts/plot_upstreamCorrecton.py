import calo_init
calo_init.add_defaults()
calo_init.parser.add_argument("--numLayers", help="Number of layers that create the first readout layer", default = 1, type = int)
calo_init.parser.add_argument("--simulation", help="Path to the simulation file, to extract the total deposited energy in calorimeter. If undefined, the beam energy is used", type = str)
calo_init.parser.add_argument("--etaValues", help="Values of eta for which the upstream material correction was calculated. They are used in the file names.", default = [0], type = float, nargs='+')
calo_init.parser.add_argument("-n","--histogramName", default="upstreamEnergy_presamplerEnergy2D", help="Name of the histogram with the histogram of energy in first layer (X-axis) and the upstrem energy (Y-axis)", type = str)
calo_init.parser.add_argument("--layerWidth", default = 2, help="Width of the layer (cm)", type = float)
calo_init.parser.add_argument("--roundBrackets", help="Use round brackets for unit", action = 'store_true')
calo_init.parser.add_argument("--preview", help="Plot preview of fits", action = 'store_true')
calo_init.parser.add_argument("--specialLabel", help="Additional label to be plotted", type=str, default = "FCC-hh simulation")
calo_init.parse_args()
calo_init.print_config()

layer = calo_init.args.numLayers
if calo_init.args.simulation:
    ifActive = True
    activeNames, check = calo_init.substitute(calo_init.args.activeFile)
    print(activeNames)
else:
    ifActive = False
width = calo_init.args.layerWidth  # cm
if calo_init.args.roundBrackets:
    energyUnit = '(GeV)'
    angleUnit = '(rad)'
else:
    energyUnit = '[GeV]'
    angleUnit = '[rad]'

from ROOT import gSystem, gROOT, TCanvas, TGraphErrors, TF1, gStyle, kRed, kBlue, kGray, TFile, TTree, TPad, TGaxis, gPad, TLine
from draw_functions import prepare_graph, prepare_divided_canvas,  prepare_single_canvas, draw_text, draw_1histogram
import numpy
from math import sqrt, ceil, floor
gStyle.SetOptFit(0000)
gStyle.SetOptStat(0000)

# graphs for the eta-dependent correction parameters
par0par0 = TGraphErrors()
par0par1 = TGraphErrors()
par1par0 = TGraphErrors()
par1par1 = TGraphErrors()

for ieta, eta in enumerate(calo_init.args.etaValues):
    # graphs for the energy-dependent correction parameters
    param0 = TGraphErrors()
    param1 = TGraphErrors()
    # calculate width of current eta bin
    if ieta < len(calo_init.args.etaValues) - 1:
        etaWidth = (calo_init.args.etaValues[ieta+1] - eta) / 2.
    else:
        etaWidth = (eta - calo_init.args.etaValues[ieta-1]) / 2.
    # loop over all energies within eta bin
    for ifile, filename in enumerate(calo_init.filenamesIn):
        # get energy values: beam energy ...
        energy = calo_init.energy(ifile)
        # ... and if possible, total energy deposited in calorimeter (which is the actual energy of reconstructed cluster before this correction is applied)
        if ifActive:
            fActive = TFile(activeNames[ifile], "READ")
            hActive = fActive.Get('total')
            resultActive = hActive.Fit("gaus","SN")
            energyGraph = resultActive.Get().Parameter(1)
            energyErrorGraph = resultActive.Get().Parameter(2)
            axisName = "Energy deposited in the calorimeter, EM scale (GeV)"
        else:
            energyGraph = energy
            energyErrorGraph = 0
            axisName = "Energy of the particle (GeV)"
        if eta>0:
            filename = filename.replace('eta0','eta'+str(eta))
        # get phi distribution and profile of upstream energy vs deposited in first layer
        f = TFile(filename, "READ")
        hCellPhi = f.Get('presamplerEnergy_phi'+str(layer))
        hCellPhi.Scale(1. / hCellPhi.GetEntries())
        hUpstremCellProfile = f.Get(calo_init.args.histogramName+str(layer)).ProfileX("profile")
        hUpstremCellProfile.SetTitle('')
        hCellPhi.SetTitle('')
        hCellPhi.GetXaxis().SetRangeUser(-0.02,0.02)

        # hUpstremCellProfile.Rebin(layer)
        fitProfile = TF1("fitProfile","pol1", 0, energy)
        result = hUpstremCellProfile.Fit(fitProfile, "SN")

        # fill energy-dependent graphs
        param0.SetPoint(ifile,energyGraph,result.Get().Parameter(0))
        param0.SetPointError(ifile,energyErrorGraph,result.Get().Error(0))
        param1.SetPoint(ifile,energyGraph,result.Get().Parameter(1))
        param1.SetPointError(ifile,energyErrorGraph,result.Get().Error(1))

        if calo_init.args.preview:
            canvProfile = prepare_single_canvas( 'upstreamEnergy_e'+str(energy)+'GeV_eta'+str(eta), 'Upstream energy vs energy in first layer for #eta = '+str(eta))
            draw_1histogram(hUpstremCellProfile, 'E_{firstLayer} '+energyUnit,'upstream energy '+energyUnit)
            hUpstremCellProfile.Fit(fitProfile, "S")
            for ibin in reversed(range(0,hUpstremCellProfile.GetXaxis().GetNbins())):
                if hUpstremCellProfile.GetBinContent(ibin) > 0:
                    lastNonEmpty = ibin
                    break
            hUpstremCellProfile.GetXaxis().SetRange(0, lastNonEmpty + 3)
            canvProfile.Update()
            canvPhi = prepare_single_canvas( 'phi_e'+str(energy)+'GeV__eta'+str(eta), 'Phi distribution of first layer deposits for #eta = '+str(eta))
            draw_1histogram(hCellPhi, '#varphi ' + angleUnit,'E_{firstLayer} ' + energyUnit)
            # Draw all labels
            if calo_init.args.specialLabel:
                canvProfile.cd()
                draw_text([calo_init.args.specialLabel], [0.57,0.88, 0.85,0.98], kGray+3, 0).SetTextSize(0.05)
                canvPhi.cd()
                draw_text([calo_init.args.specialLabel], [0.57,0.88, 0.85,0.98], kGray+3, 0).SetTextSize(0.05)
            canvProfile.cd()
            legendProfile = draw_text([str(energy)+' GeV e^{-}, B = 4T'], [0.67,0.18, 0.85,0.28], 9, 0)
            legendProfile.SetTextSize(0.05)
            legendProfile.SetTextFont(42)
            canvPhi.cd()
            legendPhi = draw_text([str(energy)+' GeV e^{-}, B = 4T'], [0.67,0.18, 0.85,0.28], 9, 0)
            legendPhi.SetTextSize(0.05)
            legendPhi.SetTextFont(42)
            canvProfile.Update()
            canvPhi.Update()
            # save canvases filled for each energy and eta
            if calo_init.output(ifile):
                canvPhi.SaveAs(calo_init.output(ifile) + "_previewPhi_eta" + str(eta) + ".png")
                canvProfile.SaveAs(calo_init.output(ifile) + "_previewProfile_eta" + str(eta) + ".png")
            else:
                canvPhi.SaveAs("upstremCorrection_previewPhi_eta"+str(eta)+"_"+str(layer*width)+"cm.png")
                canvProfile.SaveAs("upstremCorrection_previewProfile_eta"+str(eta)+"_"+str(layer*width)+"cm.png")

    # fit energy-dependent parameters
    fitP0 = TF1("fitP0","pol1", 0, energy)
    par0result = param0.Fit(fitP0, "S")
    fitP1 = TF1("fitP1","[0]+[1]/sqrt(x)", 0, energy)
    par1result = param1.Fit(fitP1, "S")
    cEnergy = prepare_divided_canvas('upstreamParams_eta'+str(eta), 'Energy upstream E=p0+p1E for eta='+str(eta), 2)
    cEnergy.cd(1)
    prepare_graph(param0, "param0", 'P0 (E);'+axisName+'; parameter P0')
    param0.Draw("aep")
    param0.GetYaxis().SetRangeUser(param0.GetYaxis().GetXmin(), param0.GetYaxis().GetXmax() * 1.2)
    cEnergy.cd(2)
    prepare_graph(param1, "param1", 'P1 (E);'+axisName+'; parameter P1')
    param1.Draw("aep")
    param1.GetYaxis().SetRangeUser(param1.GetYaxis().GetXmin(), param1.GetYaxis().GetXmax() * 1.2)
    cEnergy.Update()

    # fill energy-dependent graphs
    par0par0.SetPoint(ieta,eta,par0result.Get().Parameter(0))
    par0par0.SetPointError(ieta,etaWidth,par0result.Get().Error(0))
    par0par1.SetPoint(ieta,eta,par0result.Get().Parameter(1))
    par0par1.SetPointError(ieta,etaWidth,par0result.Get().Error(1))
    par1par0.SetPoint(ieta,eta,par1result.Get().Parameter(0))
    par1par0.SetPointError(ieta,etaWidth,par1result.Get().Error(0))
    par1par1.SetPoint(ieta,eta,par1result.Get().Parameter(1))
    par1par1.SetPointError(ieta,etaWidth,par1result.Get().Error(1))

    # save canvases filled for each eta
    if calo_init.output(ifile):
        cEnergy.SaveAs(calo_init.output(ifile) + "_energyDependent_eta"+str(eta) + ".png")
    else:
        cEnergy.SaveAs("upstreamCorrection_energyDependent_eta"+str(eta)+"_"+str(layer*width)+"cm.png")

# plot eta-dependent parameters
cParams = prepare_divided_canvas( 'etaParameters', 'Upstream energy correction parameters E_{upstream}=(p00(#eta)+p01(#eta)#cdot E) + (p10(#eta)+ p11(#eta)/#sqrt{E})#cdot E', 4 )
paramGraphs = [par0par0, par0par1, par1par0, par1par1]
parTitles = ['P00', 'P01', 'P10', 'P11']
paramValues = []
for igraph, g in enumerate(paramGraphs):
    pad = cParams.cd(igraph + 1)
    prepare_graph(g, parTitles[igraph], parTitles[igraph]+'(#eta);#eta; parameter '+parTitles[igraph])
    g.Draw("aep")
    paramValues.append(list(g.GetY()))
    pad.Update()

# save
if calo_init.output(0):
    cParams.SaveAs(calo_init.output(0)+".png")
    fileParCorrection = TFile(calo_init.output(0)+".root","RECREATE")
else:
    cParams.SaveAs("upstreamCorrection_etaDependent_"+str(layer*width)+"cm.png")
    fileParCorrection = TFile("upstreamCorrection_etaDependent_"+str(layer*width)+"cm.root","RECREATE")
fileParCorrection.cd()
for par in paramGraphs:
    par.Write()
fileParCorrection.Close()

print("===================================")
print("correction parameters: " +str(layer*width)+"cm    ")
if calo_init.output(0):
    print(calo_init.output(0))
print("===================================")
print("etaValues = ", calo_init.args.etaValues)
print("P00 = ", paramValues[0])
print("P01 = ", paramValues[1])
print("P10 = ", paramValues[2])
print("P11 = ", paramValues[3])
print("===================================")

input("Press ENTER to exit")

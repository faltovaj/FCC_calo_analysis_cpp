import calo_init
## add arguments relevant only for that script
calo_init.add_defaults()
calo_init.parser.add_argument("--particleColl", help="Name of the MC particle collection (fcc::MCParticleCollection)", type = str)
calo_init.parser.add_argument("--clusterColl", help="Name of the clusters collection (fcc::CaloClusterCollection)", type = str)
calo_init.parser.add_argument("--dEta", help="Size of the tower in eta", type = float, nargs=2)
calo_init.parser.add_argument("--maxEta", help="Maximum eta", type = float)
group = calo_init.parser.add_mutually_exclusive_group()
group.add_argument("--dPhi", help="Size of the tower in phi", type = float)
group.add_argument("--numPhi", help="Number of the towers in phi", type = int)
calo_init.parser.add_argument("inputSim", help="Additional input file name with the simulated events (to get MC particle phi,eta,E)", type = str)
calo_init.parser.add_argument("--cellColl", help="Name of the cells collection (fcc::CaloHitCollection) for the upstream energy correction", type = str)
calo_init.parser.add_argument("--correctionParams", help="Parameters for the upstream energy correction", type = float, nargs=4)
calo_init.parser.add_argument("--bitfield", help="Bitfield used to encode the IDs (from DD4hep xml, e.g. \"system:4,x:4,y:4\"", type = str)
calo_init.parse_args()

from math import pi, floor
# set of default parameters
maxEta = 1.716
maxPhi = pi-(pi/512.)
nPhi = 629 # artificially increase by 1 (odd number) - to make plots look OK
dEta = 0.01
nEta = int(2*maxEta/dEta + 1)
dPhi = 2*pi/nPhi
nameClusterCollection = "caloClusters"
nameParticlesCollection = "GenParticles"
filenamesSim, checkRegexInSimInput = calo_init.substitute(calo_init.args.inputSim)

# get parameters if passed from command line
if calo_init.args.maxEta:
    maxEta = calo_init.args.maxEta
if calo_init.args.numPhi:
    nPhi = calo_init.args.numPhi
    dPhi = 2*pi/nPhi
if calo_init.args.dPhi:
    dPhi = calo_init.args.dPhi
    nPhi = int(dPhi*2*pi)
if calo_init.args.clusterColl:
    nameClusterCollection = calo_init.args.clusterColl
if calo_init.args.particleColl:
    nameParticlesCollection = calo_init.args.particleColl
if calo_init.args.correctionParams and calo_init.args.cellColl and calo_init.args.bitfield:
    nameCellCollection = calo_init.args.cellColl
    bitfield = calo_init.args.bitfield
    doMaterialInFrontCorrection = True
    par00 = calo_init.args.correctionParams[0]
    par01 = calo_init.args.correctionParams[1]
    par10 = calo_init.args.correctionParams[2]
    par11 = calo_init.args.correctionParams[3]
else:
    doMaterialInFrontCorrection = False

from ROOT import gSystem
gSystem.Load("libCaloAnalysis")
from ROOT import SingleParticleRecoMonitors, TCanvas, TFile, gStyle, gPad, kGreen, kRed, kBlue, TColor, TF1
from draw_functions import *

# use this script for multiple files
# gStyle.SetPalette(56) # kInvertedDarkBodyRadiator
gStyle.SetPalette(73) # kCMYK
gStyle.SetOptFit(1)

for ifile, filename in enumerate(calo_init.filenamesIn):
    energy = calo_init.energy(ifile)
    filenameSim = filenamesSim[ifile] if checkRegexInSimInput else filenamesSim[0]
    print("Initial particle energy: " + str(energy) + "GeV")
    print("File with simulation results: " + filenameSim)
    print("File with reconstruction results: " + filename)
    if doMaterialInFrontCorrection:
        analysis = SingleParticleRecoMonitors(nameClusterCollection,
                                              nameParticlesCollection,
                                              energy,
                                              maxEta, # max eta
                                              maxPhi,
                                              nEta, # number of bins in eta
                                              nPhi, # number of bins in phi
                                              dEta, # tower size in eta
                                              dPhi, # tower size in phi
                                              nameCellCollection,
                                              bitfield,
                                              "cell", # layer field name in the bitfield
                                              1, # Id of first layer
                                              4, # Id of last layer that counts as first (= 4*2cm = 8cm layer)
                                              0.168, # sampling fraction of the first layer, if calibrated cells were given
                                              par00,
                                              par01,
                                              par10,
                                              par11)
    else:
        analysis = SingleParticleRecoMonitors(nameClusterCollection,
                                              nameParticlesCollection,
                                              energy,
                                              maxEta, # max eta
                                              maxPhi,
                                              nEta, # number of bins in eta
                                              nPhi, # number of bins in phi
                                              dEta, # tower size in eta
                                              dPhi)# tower size in phi
    analysis.loop(filenameSim, filename, calo_init.verbose)
    # retrieve histograms to draw them
    hEn = analysis.hEn
    hEnTotal = analysis.hEnTotal
    hEnCorrected = analysis.hEnCorr
    hEnFirstLayer = analysis.hEnFirstLayer
    hEnUpstream = analysis.hEnUpstream
    hEnFncPhi = analysis.hEnFncPhi
    hEta = analysis.hEta
    hPhi = analysis.hPhi
    hPhiFncPhi = analysis.hPhiFncPhi
    hEtaFncEta = analysis.hEtaFncEta
    hNo = analysis.hNo
    hNoFncPhi = analysis.hNoFncPhi
    hNoFncEta = analysis.hNoFncEta
    hEnMoreClu = analysis.hEnMoreClu
    hEnDiffMoreClu = analysis.hEnDiffMoreClu
    hEtaMoreClu = analysis.hEtaMoreClu
    hEtaDiffMoreClu = analysis.hEtaDiffMoreClu
    hPhiMoreClu = analysis.hPhiMoreClu
    hPhiDiffMoreClu = analysis.hPhiDiffMoreClu
    hRDiffMoreClu = analysis.hRDiffMoreClu
    h1dset1 = [hEn, hEnCorrected,  hEnFirstLayer, hEnUpstream, hEnTotal, hEta, hPhi, hNo, hEnMoreClu, hEtaMoreClu, hPhiMoreClu]
    for h in h1dset1:
        h.SetMarkerColor(kBlue+3)
        h.SetFillColor(39)
        h.SetLineColor(39)
        h.SetMarkerSize(2.2)
    h1dset2 = [hEnDiffMoreClu, hEtaDiffMoreClu, hPhiDiffMoreClu, hRDiffMoreClu]
    for h in h1dset2:
        h.SetMarkerColor(39)
        h.SetMarkerStyle(21)
        h.SetLineColor(39)
        h.SetMarkerSize(1.5)
    h2d = [hEnFncPhi, hPhiFncPhi, hEtaFncEta, hNoFncPhi, hNoFncEta]
    hNo.GetXaxis().SetNdivisions(7)
    hNoFncPhi.GetYaxis().SetNdivisions(7)
    hNoFncEta.GetYaxis().SetNdivisions(7)
    hPhi.GetXaxis().SetRangeUser(-50*dPhi,50*dPhi)
    hPhiFncPhi.GetYaxis().SetRangeUser(-50*dPhi,50*dPhi)

    # fit functions
    fitEnergy = TF1('fitEnergy','gaus',0.8*energy,1.2*energy)
    resultFitEn = hEn.Fit('fitEnergy','S')
    resultFitEnCorr = hEnCorrected.Fit('fitEnergy','SRN')
    fitEnergy2 = TF1('fitEnergy2','gaus',resultFitEn.Get().Parameter(1)-2.*resultFitEn.Get().Parameter(2),
                     resultFitEn.Get().Parameter(1)+2.*resultFitEn.Get().Parameter(2))
    resultFitEn2 = hEn.Fit('fitEnergy2','SR')
    fitEnergy2Corr = TF1('fitEnergy2Corr','gaus',resultFitEnCorr.Get().Parameter(1)-2.*resultFitEnCorr.Get().Parameter(2),
                     resultFitEnCorr.Get().Parameter(1)+2.*resultFitEnCorr.Get().Parameter(2))
    resultFitEn2Corr = hEnCorrected.Fit('fitEnergy2Corr','SR')
    fitEta = TF1('fitEta','gaus',-10*dEta,10*dEta)
    resEta = hEta.Fit('fitEta','SR')
    fitPhi = TF1('fitPhi','gaus',-10*dPhi,10*dPhi)
    resPhi = hPhi.Fit('fitPhi','SR')


    canv = TCanvas('ECal_monitor_plots_e'+str(energy)+'GeV', 'ECal', 2000, 1600 )
    canv.Divide(4,4)
    canv.cd(1)
    if doMaterialInFrontCorrection:
        hEnCorrected.Draw('bar')
        draw_text(["energy: "+str(round(resultFitEn2Corr.Get().Parameter(1),1))+" GeV",
                   "           "+str(round(resultFitEn2Corr.Get().Parameter(1),1)/energy*100.)+" %",
                   "resolution: "+str(round(resultFitEn2Corr.Get().Parameter(2)/resultFitEn2Corr.Get().Parameter(1)*100,1))+" %"],
                  [0.1,0.6,0.4,0.9])
    else:
        hEn.Draw('bar')
        draw_text(["energy: "+str(round(resultFitEn2.Get().Parameter(1),1))+" GeV",
                   "           "+str(round(resultFitEn2.Get().Parameter(1),1)/energy*100.)+" %",
                   "resolution: "+str(round(resultFitEn2.Get().Parameter(2)/resultFitEn2.Get().Parameter(1)*100,1))+" %"],
                  [0.1,0.6,0.4,0.9])
    canv.cd(2)
    draw_hist2d(hEnFncPhi)
    canv.cd(5)
    hEta.Draw('bar')
    canv.cd(6)
    draw_hist2d(hEtaFncEta)
    canv.cd(9)
    hPhi.Draw('bar')
    canv.cd(10)
    draw_hist2d(hPhiFncPhi)
    canv.cd(13)
    hNo.Draw('hist')
    hNo.Draw('text0 same')
    canv.cd(14)
    draw_hist2d(hNoFncPhi)
    canv.cd(15)
    draw_hist2d(hNoFncEta)
    canv.cd(3)
    hEnMoreClu.Draw('bar')
    canv.cd(7)
    hEtaMoreClu.Draw('bar')
    canv.cd(11)
    hPhiMoreClu.Draw('bar')
    canv.cd(4)
    hEnDiffMoreClu.Draw()
    canv.cd(8)
    hEtaDiffMoreClu.Draw()
    canv.cd(12)
    hPhiDiffMoreClu.Draw()
    canv.cd(16)
    hRDiffMoreClu.Draw()

    if calo_init.output(ifile):
        canv.Print(calo_init.output(ifile)+".gif")
        plots = TFile(calo_init.output(ifile)+".root","RECREATE")
    else:
        canv.Print("Reco_monitor.gif")
        plots = TFile("Reco_monitor_plots.root","RECREATE")

    plots.cd()
    for hset in [h1dset1, h1dset2, h2d]:
        for h in hset:
            h.Write()
    plots.Close()

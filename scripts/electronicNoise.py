from ROOT import TH1F, TCanvas, TLegend, TFile, gStyle, gPad
import itertools

#filename = "50Ohm_default"
filename = "50Ohm_default_differentTraces"
fIn = TFile("capacitances_"+filename+".root","r")

nLayers = 8
etaMax = 1.5
nbins = 150
noisePerCapacitance = 0.0397372130805
print "noise per capacitance ", noisePerCapacitance, " MeV/pF"
GeV = 1000.

gStyle.SetOptStat(0)

#Graphs with capacitances
hCapShield = []
hCapTrace = []
hCapDetector = []
hCapTotal = []
# electronic noise histograms 
h_elecNoise_fcc = [] # total noise from shield + trace + detector capacitance
h_elecNoise_withoutTraceCap = [] # total noise without trace capacitance (as in ATLAS - trace cap. can be neglected): from shield + detector capacitance
h_elecNoise_shield = []
h_elecNoise_trace = []
h_elecNoise_detector = []

for i in range (0, nLayers):
    #Prepare electronic noise histograms
    h_elecNoise_fcc.append( TH1F() )
    h_elecNoise_fcc[i].SetLineWidth(3)
    h_elecNoise_fcc[i].SetLineColor(i+1)
    h_elecNoise_fcc[i].SetBins(nbins, 0., etaMax)
    h_elecNoise_fcc[i].SetTitle("Total electronic noise; |#eta|; Electronic noise [GeV]")
    h_elecNoise_fcc[i].SetName("h_elecNoise_fcc_"+str(i+1))

    h_elecNoise_withoutTraceCap.append( TH1F() )
    h_elecNoise_withoutTraceCap[i].SetLineWidth(3)
    h_elecNoise_withoutTraceCap[i].SetLineColor(i+1)
    h_elecNoise_withoutTraceCap[i].SetBins(nbins, 0., etaMax)
    h_elecNoise_withoutTraceCap[i].SetTitle("Electronic noise without trace capacitance; |#eta|; Electronic noise [GeV]")
    h_elecNoise_withoutTraceCap[i].SetName("h_elecNoise_withoutTraceCap_"+str(i+1))

    h_elecNoise_shield.append( TH1F() )
    h_elecNoise_shield[i].SetLineWidth(3)
    h_elecNoise_shield[i].SetLineColor(i+1)
    h_elecNoise_shield[i].SetBins(nbins, 0., etaMax)
    h_elecNoise_shield[i].SetTitle("Electronic noise - shields; |#eta|; Electronic noise [GeV]")
    h_elecNoise_shield[i].SetName("h_elecNoise_shield_"+str(i+1))

    h_elecNoise_trace.append( TH1F() )
    h_elecNoise_trace[i].SetLineWidth(3)
    h_elecNoise_trace[i].SetLineColor(i+1)
    h_elecNoise_trace[i].SetBins(nbins, 0., etaMax)
    h_elecNoise_trace[i].SetTitle("Electronic noise -traces; |#eta|; Electronic noise [GeV]")
    h_elecNoise_trace[i].SetName("h_elecNoise_trace_"+str(i+1))

    h_elecNoise_detector.append( TH1F() )
    h_elecNoise_detector[i].SetLineWidth(3)
    h_elecNoise_detector[i].SetLineColor(i+1)
    h_elecNoise_detector[i].SetBins(nbins, 0., etaMax)
    h_elecNoise_detector[i].SetTitle("Electronic noise - detector; |#eta|; Electronic noise [GeV]")
    h_elecNoise_detector[i].SetName("h_elecNoise_detector_"+str(i+1))
  
    #Read graphs from files
    nameShield = "hCapacitance_shields"+str(i)
    hCapShield.append(fIn.Get(nameShield))
    nameTrace = "hCapacitance_traces"+str(i)
    hCapTrace.append(fIn.Get(nameTrace))
    nameDetector = "hCapacitance_area"+str(i)
    hCapDetector.append(fIn.Get(nameDetector))
    #Sum the capacitances
    hCapTotal.append( TH1F() )
    hCapTotal[i].SetBins(nbins, 0, etaMax)
    hCapTotal[i].SetLineColor(i+1)
    hCapTotal[i].SetLineWidth(3)
    hCapTotal[i].SetTitle("Total capacitance; |#eta|; Capacitance [pF]")
    hCapTotal[i].SetName("hCapacitance"+str(i))
    #hCapTotal[i] = hCapShield[i] + hCapTrace[i] + hCapDetector[i]

    for ibin in range(0, hCapShield[i].GetNbinsX()):
        capShield = hCapShield[i].GetBinContent(ibin)
        capTrace = hCapTrace[i].GetBinContent(ibin)
        capDetector = hCapDetector[i].GetBinContent(ibin)
        #total capacitance
        hCapTotal[i].SetBinContent( ibin, capShield + capTrace + capDetector )
        #noise
        noise = noisePerCapacitance * ( capShield + capTrace + capDetector ) / GeV
        noiseWithoutTrace = noisePerCapacitance * ( capShield + capDetector ) / GeV
        noiseShield = noisePerCapacitance * capShield / GeV
        noiseTrace = noisePerCapacitance * capTrace / GeV
        noiseDetector = noisePerCapacitance * capDetector / GeV
        #fill histogram
        h_elecNoise_fcc[i].SetBinContent(ibin, noise)
        h_elecNoise_withoutTraceCap[i].SetBinContent(ibin, noiseWithoutTrace)
        h_elecNoise_shield[i].SetBinContent(ibin, noiseShield)
        h_elecNoise_trace[i].SetBinContent(ibin, noiseTrace)
        h_elecNoise_detector[i].SetBinContent(ibin, noiseDetector)
        if ibin==1:
            print "layer %d" %i, "eta==0: capacitance %.0f pF," %( capShield + capTrace + capDetector ), "total elec. noise %.4f GeV" %noise, "elec. noise without trace cap. %.4f GeV" %noiseWithoutTrace

maximumCap = 1500.
maximumNoise = 0.08

fSave = TFile("elecNoise_"+filename+".root","RECREATE")

cCapacitance = TCanvas("cCapacitance","Capacitance per cell",800,600)
cCapacitance.cd()
legend = TLegend(0.1,0.5,0.3,0.9)
legend.SetHeader("Longitudinal layers")
for i, h in enumerate(hCapTotal):
    h.SetMinimum(0)
    h.SetMaximum(maximumCap)
    h.GetYaxis().SetTitleOffset(1.4)
    if i == 0:
        h.Draw("")
    else:
        h.Draw("same")
    legend.AddEntry(h,str(i)+" layer","l")

print maximumCap

for h in hCapTotal:
    h.SetMinimum(0.)
    h.SetMaximum(maximumCap*1.2)
    h.GetYaxis().SetTitleOffset(1.4)
    h.Write()

legend.Draw()
cCapacitance.Update()
cCapacitance.Write()

for h in itertools.chain(h_elecNoise_fcc, h_elecNoise_withoutTraceCap, h_elecNoise_shield, h_elecNoise_trace, h_elecNoise_detector):
    h.SetMinimum(0.)
    h.SetMaximum(maximumNoise)
    h.GetYaxis().SetTitleOffset(1.4)
    h.Write()

cNoise = TCanvas("cNoise","Electronic noise per cell",800,600)
cNoise.cd()
for i, h in enumerate(h_elecNoise_fcc):
    if i == 0:
        h.Draw("")
    else:
        h.Draw("same")

legend.Draw()
cNoise.Update()
cNoise.Write()

cNoiseWithoutTrace = TCanvas("cNoiseWithoutTrace","Electronic noise without trace cap. per cell",800,600)
cNoiseWithoutTrace.cd()
for i, h in enumerate(h_elecNoise_withoutTraceCap):
    if i == 0:
        h.Draw("")
    else:
        h.Draw("same")

legend.Draw()
cNoiseWithoutTrace.Update()
cNoiseWithoutTrace.Write()

legendP = TLegend(0.1,0.6,0.43,0.9)
legendP.SetHeader("Capacitance")

cCapParts = TCanvas("cCapParts","",1200,1000)
cCapParts.Divide(3,3)    
for i in range (0, nLayers):
    cCapParts.cd(i+1)
    if i<2:
        hCapTotal[i].SetMaximum(400)
    hCapTotal[i].SetTitle("Layer "+str(i))
    hCapTotal[i].GetXaxis().SetTitleSize(0.045)
    hCapTotal[i].GetYaxis().SetTitleSize(0.045)
    hCapTotal[i].GetYaxis().SetTitleOffset(1.15)
    hCapTotal[i].SetLineColor(1)
    hCapShield[i].SetLineColor(2)
    hCapTrace[i].SetLineColor(3)
    hCapDetector[i].SetLineColor(4)
    hCapTotal[i].Draw()
    hCapShield[i].Draw("same")
    hCapTrace[i].Draw("same")
    hCapDetector[i].Draw("same")
    gPad.Update()
    if i==0:
        legendP.AddEntry(hCapTotal[i],"total cap.","l")
        legendP.AddEntry(hCapShield[i],"shield cap.","l")
        legendP.AddEntry(hCapTrace[i],"trace cap.","l")
        legendP.AddEntry(hCapDetector[i],"detector cap.","l")
    legendP.Draw()

cCapParts.Write()

closeInput = raw_input("Press ENTER to exit")

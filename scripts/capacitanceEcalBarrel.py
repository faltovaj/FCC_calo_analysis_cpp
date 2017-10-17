from ROOT import TH1F, TF1, TF2, TCanvas, TLegend, TFile, Double, gStyle
from math import ceil, sin, cos, atan, exp, log, tan, pi

deltaEta = 0.01
maxEta = 1.5
numEta = int(ceil(maxEta/deltaEta))
layerThickness = [26.9] + [120.9] * 7 #mm
tracesPerLayer = [1, 0, 0, 1, 2, 3, 4, 5]
tracesLength0 = [26.9 * 0, 26.9 * 1,
          120.9 * 5, 120.9 * 4, 120.9 * 3,
          120.9 * 2, 120.9 * 1, 120.9 * 0]

pcbThickness = 1.2 #mm
pcbThicknessNoHV = 1.0 #mm (0.1 mm between PCB inside & HV layer on both sides)
Nplanes = 1408
angle = 50./180. * pi #inclination angle
passiveThickness = 2.0 #mm

filename = "capacitances_50Ohm_default.root"

#constants:
# distance from signal trace to shield (HS) - from impedance vs. trace width vs. distance to ground layer 2D plot (Z = 50 Ohm)
# trace width (W) - min value
# trace thickness (T) - min value
# distance from shield to the edge of PCB
#http://www.analog.com/media/en/training-seminars/design-handbooks/Basic-Linear-Design/Chapter12.pdf, page 40
#signal trace
hs=0.17 # for 50 Ohm
#hs = 0.09 # for 33 Ohm
w=0.127
t=0.035
#shield
nmult_ws = 2 # how many times are shield larger than gap hs (consider 2 and 4 - extreme case)
ws = nmult_ws*hs
hm = (pcbThicknessNoHV - 5*t - 2*hs)/2.

epsilonR = 4.0 # PCB

#conversion factor: 1 inch = 25/4 mm
inch2mm = 25.4

#multiplicative factor
# factor 2 for traces, 2 for 2 widths between planes (dphi_segmentation = 0.01, dphi_planes = 0.004)
nmult = 2*2

#Detector area calculation
radius = [1920,
          1920 + 26.9,
          1920 + 26.9 + 120.9,
          1920 + 26.9 + 120.9 * 2,
          1920 + 26.9 + 120.9 * 3,
          1920 + 26.9 + 120.9 * 4,
          1920 + 26.9 + 120.9 * 5,
          1920 + 26.9 + 120.9 * 6,
          1920 + 26.9 + 120.9 * 7]
epsilonRLAr = 1.5 # LAr at 88 K
epsilon0 = 8.854/1000. #pF/mm

gStyle.SetOptStat(0)

cImpedance = TCanvas("cImpedance","",600,800)
cImpedance.Divide(1,2)
cImpedance.cd(1)
fImpedance = TF2("fImpedance","60/sqrt([0])*log(1.9*(2*x+[1])/(0.8*y+[1]))",0.04,0.2,0.04,0.2)
fImpedance.SetTitle("Impedance vs trace width and distance to ground")
fImpedance.SetParameters(epsilonR, t)
fImpedance.Draw("colz")
fImpedance.GetXaxis().SetTitle("Distance to ground [mm]")
fImpedance.GetYaxis().SetTitle("Trace width [mm]")
cImpedance.cd(2)
fImpedance1D = TF1("fImpedance1D","60/sqrt([0])*log(1.9*(2*x+[1])/(0.8*[2]+[1]))",0.04,0.2)
fImpedance1D.SetTitle("Impedance vs distance to ground")
fImpedance1D.SetParameters(epsilonR, t, w)
fImpedance1D.Draw()
fImpedance1D.GetXaxis().SetTitle("Distance to ground [mm]")
fImpedance1D.GetYaxis().SetTitle("Impedance [#Omega]")

hCapTrace = []
hCapShield = []
hCapArea = []
for i in range (0, len(layerThickness)):
    #traces
    hCapTrace.append(TH1F())
    hCapTrace[i].SetBins(numEta, 0, maxEta)
    hCapTrace[i].SetLineColor(i+1)
    hCapTrace[i].SetLineWidth(2)
    hCapTrace[i].SetTitle("Capacitance of traces; |#eta|; Capacitance [pF]")
    hCapTrace[i].SetName("hCapacitance_traces"+str(i))
    #shields
    hCapShield.append(TH1F())
    hCapShield[i].SetBins(numEta, 0, maxEta)
    hCapShield[i].SetLineColor(i+1)
    hCapShield[i].SetLineWidth(2)
    hCapShield[i].SetTitle("Capacitance of shields; |#eta|; Capacitance [pF]")
    hCapShield[i].SetName("hCapacitance_shields"+str(i))
    #area
    hCapArea.append(TH1F())
    hCapArea[i].SetBins(numEta, 0, maxEta)
    hCapArea[i].SetLineColor(i+1)
    hCapArea[i].SetLineWidth(2)
    hCapArea[i].SetTitle("Capacitance of detector area; |#eta|; Capacitance [pF]")
    hCapArea[i].SetName("hCapacitance_area"+str(i))

cTrace = TCanvas("cTrace","",600,400)
cShield = TCanvas("cShield","",600,400)
cArea = TCanvas("cArea","",600,400)

legend = TLegend(0.1,0.5,0.3,0.9)
legend.SetHeader("Longitudinal layers")
for i in range (0, len(layerThickness)):
    for index in range(0, numEta):
        eta = index * deltaEta
        #Trace capacitance (stripline)
        traceLength = tracesLength0[i] / (sin(2. * atan(exp(-eta))))
        logStripline = log(3.1 * hs / (0.8 * w + t))
        capacitanceTrace = nmult * 1 / inch2mm * 1.41 * epsilonR / logStripline * traceLength
        hCapTrace[i].SetBinContent(index+1, capacitanceTrace)
    
        #Shield capacitance (microstrip)
        shieldLength = tracesPerLayer[i] * layerThickness[i] / (sin(2. * atan(exp(-eta))))
        logMicrostrip = log(5.98 * hm / (0.8 * ws + t))
        capacitanceShield = 2 * nmult * 1 / inch2mm * 0.67 * (epsilonR + 1.41) / logMicrostrip * shieldLength
        hCapShield[i].SetBinContent(index+1, capacitanceShield)

        #Detector area (C = epsilon*A/d)
        area = ( radius[i] * ( 1 / (tan(2. * atan(exp(- (index + 1) * deltaEta)))) -  1 / (tan(2. * atan(exp(- index * deltaEta))) ) )
                 + radius[i + 1] * ( 1 / (tan(2. * atan(exp(- (index + 1) * deltaEta)))) -  1 / (tan(2. * atan(exp(- index * deltaEta))) ) )
                 ) / 2. * (radius[i+1] - radius[i])
        distance = (radius[i+1] + radius[i]) / 2. * pi / Nplanes * cos (angle) - pcbThickness / 2. - passiveThickness / 2.
        capacitanceArea = nmult * epsilon0 * epsilonRLAr * area / distance
        hCapArea[i].SetBinContent(index+1, capacitanceArea)
        if index==0:
            print "layer %d" %i, "eta==0: capacitanceTrace %.0f mm," %capacitanceTrace, "capacitanceShield %.0f pF" %capacitanceShield, "capacitanceArea %.0f pF" %capacitanceArea
            #, "distance %.1f mm" %distance

    #Draw
    cTrace.cd()
    if i == 0:
        hCapTrace[i].Draw()
    else:
        hCapTrace[i].Draw("same")
    legend.AddEntry(hCapTrace[i],str(i)+" layer","l")
    cShield.cd()
    if i == 0:
        hCapShield[i].Draw()
    else:
        hCapShield[i].Draw("same")
    cArea.cd()
    if i == 0:
        hCapArea[i].Draw()
    else:
        hCapArea[i].Draw("same")

maximum = 1600

plots = TFile(filename,"RECREATE")

for i in range (0, len(layerThickness)):
    hCapTrace[i].SetMinimum(0)
    hCapTrace[i].SetMaximum(maximum*1.2)
    hCapTrace[i].Write()
    hCapShield[i].SetMinimum(0)
    hCapShield[i].SetMaximum(maximum*1.2)
    hCapShield[i].Write()
    hCapArea[i].SetMinimum(0)
    hCapArea[i].SetMaximum(maximum*1.2)
    hCapArea[i].Write()

cTrace.cd()
legend.Draw()
cTrace.Update()
cTrace.Write()
cShield.cd()
legend.Draw()
cShield.Update()
cShield.Write()
cArea.cd()
legend.Draw()
cArea.Update()
cArea.Write()

fImpedance.Write()
fImpedance1D.Write()

closeInput = raw_input("Press ENTER to exit")

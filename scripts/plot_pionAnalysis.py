import calo_init
## add arguments relevant only for that script
calo_init.add_defaults()
calo_init.parser.add_argument("--particleColl", help="Name of the MC particle collection (fcc::MCParticleCollection)", type = str)
calo_init.parser.add_argument("--clusterColl", help="Name of the clusters collection (fcc::CaloClusterCollection)", type = str)
calo_init.parser.add_argument("--cellColl", help="Name of the cells collection (fcc::CaloHitCollection) for the upstream energy correction", type = str)
calo_init.parser.add_argument("--bitfield", help="Bitfield used to encode the IDs (from DD4hep xml, e.g. \"system:4,x:4,y:4\"", type = str)
calo_init.parse_args()


nameClusterCollection = "EcalClusters"
nameParticlesCollection = "GenParticles"
nameCellCollection = "ECalPositions"
bitfield = "system:4,cryo:1,type:3,subtype:3,cell:6,eta:9,phi:10"

from math import pi, floor
# set of default parameters
maxEta = 1.716
maxPhi = pi-(pi/512.)
nPhi = 512 # artificially increase by 1 (odd number) - to make plots look OK
dEta = 0.01
nEta = int(2*maxEta/dEta + 1)
dPhi = 2*pi/nPhi

# get parameters if passed from command line
if calo_init.args.clusterColl:
    nameClusterCollection = calo_init.args.clusterColl
if calo_init.args.particleColl:
    nameParticlesCollection = calo_init.args.particleColl
if calo_init.args.cellColl:
    nameCellCollection = calo_init.args.cellColl
if calo_init.args.bitfield:
    bitfield = calo_init.args.bitfield
doMaterialInFrontCorrection = True
par00 = 0.1314
par01 = 0.000498
par10 = 4.312
par11 = 4.068

from ROOT import gSystem
gSystem.Load("libCaloAnalysis")
from ROOT import PionAnalysis, TCanvas, TFile, gStyle, gPad, kGreen, kRed, kBlue, TColor, TF1
from draw_functions import *

# use this script for multiple files
# gStyle.SetPalette(56) # kInvertedDarkBodyRadiator
gStyle.SetPalette(73) # kCMYK
gStyle.SetOptFit(1)

for ifile, filename in enumerate(calo_init.filenamesIn):
    energy = calo_init.energy(ifile)
    print "Initial particle energy: " + str(energy) + "GeV"
    print "File with reconstruction results: " + filename
    if doMaterialInFrontCorrection:
        analysis = PionAnalysis(nameClusterCollection,
                                              nameParticlesCollection,
                                              energy,
                                              maxEta, # max eta
                                              maxPhi,
#                                               nEta, # number of bins in eta                                                                                                                         
#                                               nPhi, # number of bins in phi                                                                                                                                                     
#                                               dEta, # tower size in eta                                                                                                                                                      
#                                               dPhi, # tower size in phi                                                                                                                                                   
                                              nameCellCollection,
                                              bitfield,
                                              "cell", # layer field name in the bitfield
                                              1, # Id of first layer
                                              1, # Id of last layer that counts as first (= 4*2cm = 8cm layer)
                                              0.168, # sampling fraction of the first layer, if calibrated cells were given
                                              par00,
                                              par01,
                                              par10,
                                              par11)
    else:
        analysis = PionAnalysis(nameClusterCollection,
                                              nameParticlesCollection,
                                              energy,
                                              maxEta, # max eta
                                              maxPhi,
                                              nEta, # number of bins in eta
                                              nPhi, # number of bins in phi
                                              dEta, # tower size in eta
                                              dPhi)# tower size in phi
    analysis.loop(filename, calo_init.verbose)
    # retrieve histograms to draw them
    hEnPi = analysis.hPiEnergy
#    hEnPi.Draw()
    print hEnPi.GetEntries(), "  ", hEnPi.GetMean()

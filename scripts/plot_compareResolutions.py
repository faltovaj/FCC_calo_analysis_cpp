import calo_init
calo_init.add_defaults()
calo_init.parser.add_argument("--legend","-l",default=[],type=str,nargs='+')
calo_init.parser.add_argument("--title","-t",default="Energy resolution",type=str)
calo_init.parser.add_argument("-m","--axisMax", help="Maximum of the axis", type = float)
calo_init.parser.add_argument("--sequentialColours", "--colours", help="If Gradient of colours should be used insted of ROOT standard", default=False, action='store_true')
calo_init.parser.add_argument("--noLinearity", help="If linearity plot should not be drawn", action = 'store_true')
calo_init.parser.add_argument("--specialLabel", help="Additional label to be plotted", type=str, default = "FCC-hh simulation")
calo_init.parse_args()
calo_init.print_config()

print("Draw linearity: ", calo_init.args.noLinearity)

from ROOT import gSystem, gROOT, TCanvas, TGraphErrors, TF1, TFile, TColor, TPad, TGaxis, TAxis
from draw_functions import prepare_graph, prepare_second_graph, prepare_single_canvas, prepare_double_canvas, draw_text

graphsRes=[]
graphsLin=[]
constTerm=[]
samplTerm=[]

# first loop over all files to get the resolutions
for ifile, filename in enumerate(calo_init.filenamesIn):
    energy = calo_init.energy(ifile)
    f = TFile(filename, "READ")
    graphsRes.append(f.Get("resolution"))
    if not calo_init.args.noLinearity:
        graphsLin.append(f.Get("linearity"))
    for param in f.Get("params"):
        constTerm.append(param.const)
        samplTerm.append(param.sampl)

# Prepare canvas
if not calo_init.args.noLinearity:
    factor = 2 # meaning the upper plot is twice smaller than the bottom plot
    cRes, padRes, padLin = prepare_double_canvas("resolution","Energy resolution", factor)
else:
    cRes = prepare_single_canvas("resolution","Energy resolution")

# Prepare graphs: set colours, axis range, ...
if calo_init.args.sequentialColours and len(graphsRes) < 8:
    if calo_init.filenamesIn[0].find("Bfield1"):
        colour = ['#b84341','#bc5e42','#be7a42','#c09644','#c2b246','#b9c448','#a0c54a']
        colour = [colour[len(colour)-i-1] for i,c in enumerate(colour)]
        colour = [TColor.GetColor(c) for c in colour]
    else:
        colour = ['#4aa7bf','#4c90c0','#4d79c2','#4f62c4','#5851c6','#7253c7','#8d55c9']
        colour = [TColor.GetColor(c) for c in colour]
else:
    colour = [c for c in range(1,100)]
print(colour)
minima=[]
maxima=[]
for i,g in enumerate(graphsRes):
    prepare_graph(g, g.GetName(), g.GetTitle(), colour[i])
    g.GetFunction("res").SetLineColor(colour[i])
    minima.append(g.GetYaxis().GetXmin())
    maxima.append(g.GetYaxis().GetXmax())
if calo_init.args.axisMax:
    graphsRes[0].GetYaxis().SetRangeUser(0, calo_init.args.axisMax)
else:
    graphsRes[0].GetYaxis().SetRangeUser(0.8*min(minima),0.55*max(maxima))
if not calo_init.args.noLinearity:
    for gLin, gRes in zip(graphsLin, graphsRes):
        prepare_second_graph(gLin, gRes, gLin.GetName(), gLin.GetTitle(), factor)

# Draw graphs
graphsRes[0].Draw("aep")
for g in graphsRes[1:]:
    g.Draw("sameep")
if not calo_init.args.noLinearity:
    padLin.cd()
    graphsLin[0].Draw("aep")
    for g in graphsLin[1:]:
        g.Draw("sameep")

# Set proper names for the plot legend
if len(calo_init.args.legend) > 1:
    graphTitles = calo_init.args.legend
elif calo_init.args.legend:
    graphTitles, check = calo_init.substitute(calo_init.args.legend[0])
else:
    graphTitles=[]
for ileg, legend in enumerate(graphTitles):
    graphTitles[ileg] = legend.replace("formula","#frac{#sigma_{E}}{E} = #frac{"+str(round(samplTerm[ileg]*100.,2))+"%}{#sqrt{E}} + "+str(round(constTerm[ileg]*100.,2))+"%")
graphTitles = ['#color['+str(colour[i])+']{'+t+'}' for i,t in enumerate(graphTitles)]

# Draw all labels
if not calo_init.args.noLinearity:
    padRes.cd()
    draw_text(graphTitles, [0.65,0.95 - 0.08 * len(graphTitles),0.95,0.95], 1, 0).SetTextSize(0.05)
else:
    draw_text(graphTitles, [0.65,0.85 - 0.08 * len(graphTitles),0.95,0.85], 1, 0).SetTextSize(0.05)
draw_text(["energy resolution"], [0.2,0.88, 0.4,0.98], 1, 0).SetTextSize(0.05)
if calo_init.args.noLinearity and calo_init.args.specialLabel:
    draw_text([calo_init.args.specialLabel], [0.67,0.88, 0.95,0.98], 1, 0).SetTextSize(0.05)
if not calo_init.args.noLinearity:
    padLin.cd()
    draw_text(["linearity"], [0.2,0.78, 0.4,0.88], 1, 0).SetTextSize(0.05*factor)
    if calo_init.args.specialLabel:
        draw_text([calo_init.args.specialLabel], [0.67,0.78, 0.95,0.88], 1, 0).SetTextSize(0.05*factor)

cRes.Update()

# Save canvas
if calo_init.output(0):
    cRes.SaveAs(calo_init.output(0)+".gif")
else:
    cRes.SaveAs("energy_resolution_plots.gif")

raw_input("Press ENTER to exit")

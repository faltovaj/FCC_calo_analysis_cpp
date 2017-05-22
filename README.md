FCC_calo_analysis
=====================

Analysis tools for calorimetry:
 - hits/cells/cluster read from EDM in c++
 - setup and results of analysis (drawing etc.) handled from python

# Initialization

Dependencies:
 - ROOT
 - PODIO
 - fcc-edm
taken from the FCC software stack:

~~~{.sh}
source init.sh
~~~

# Build

~~~{.sh}
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install ..
make -j 8 install
cd ..
~~~

# Run analysis

## Common: configuration

To see the help message:

~~~{.sh}
python scripts/calo_init.py --help
~~~

The required argument is the input file name and the energy (vector):

~~~{.sh}
python scripts/calo_init.py <INPUT_FILE_NAME>.root [energy]
~~~

User may also specify any of the optional arguments.

In particular, if user wants to analyse multiple files:

~~~{.sh}
python scripts/calo_init.py <INPUT_FILE_NAME_CONTAINING_?_CHARACTER>.root ENERGY -r [LIST_OF_SUBSTITUTIONS]
~~~

where `-r` indicates that a regex should be used, substituting the `?` character in the input (or output) file name with the elements of the list of substitutions.

For instance,
~~~{.sh}
python scripts/calo_init.py 'output_10GeV_?.root' 10 -r firstPart secondPart
~~~

results in:
~~~{.sh}
Input file: output_10GeV_firstPart.root
Input file: output_10GeV_secondPart.root
Energy of initial particle: [10] GeV
~~~

If instead of of a list, "energy" is used after "-r" option, the '?' character gets substituted with the energy of the particles.

For instance,
~~~{.sh}
python scripts/calo_init.py 'output_?GeV.root' 10 50 100 -r energy
~~~

results in:
~~~{.sh}
Input file: output_10GeV.root
Input file: output_50GeV.root
Input file: output_100GeV.root
Energy of initial particle: [10, 50, 100] GeV
~~~


## Example: `SimpleAnalysis`

Required arguments:
 - input file name
 - energy
 - sampling factor

~~~{.sh}
python scripts/test_macro_simple.py <INPUT_FILE_NAME>.root ENERGY SF
~~~

## Shower profiles: `ShowerProfiles`

Required arguments:
 - input file name
 - energy
 - sampling factor

Analysis is done for one file.

~~~{.sh}
python scripts/test_macro_profiles.py <INPUT_FILE_NAME>.root ENERGY SF
~~~

## Cells: `CellAnalysis`

Required arguments:
 - input file name
 - energy
 - sf

Analysis is done for one file.

~~~{.sh}
python scripts/test_macro_cells.py <INPUT_FILE_NAME>.root ENERGY SF
~~~

## Reconstructed event preview: `ReconstructionExample`

Preview of the reconstructed cluster.
The first plot contains all towers, the second the centre of the cluster and the third one all the cells associated to the cluster.
Additionally the windows used in the reconstrucion are drawn.

Required arguments:
 - input file name
 - energy

Analysis is done for all the files (all energies), but for one event from each (default: first event).

~~~{.sh}
python scripts/plot_recoExample.py <INPUT_FILE_NAME>.root ENERGY
~~~

Input file name and the energy are required.

List of additional options:
~~~{.sh}
   -r REGEX [REGEX ...], --regex REGEX [REGEX ...]
                        String to insert in place of '?' character in file
                        names ("energy" inserts the values of energies)
  -o OUTPUT, --output OUTPUT
                        Output file name
  -v, --verbose         Verbose
    --clusterColl CLUSTERCOLL
                        Name of the clusters collection
                        (fcc::CaloClusterCollection)
  --positionColl POSITIONCOLL
                        Name of the positioned cells collection
                        (fcc::PositionedCaloHitCollection)
  --windowSeed WINDOWSEED WINDOWSEED
                        Size of the window used for seeding [eta,phi]
  --windowPos WINDOWPOS WINDOWPOS
                        Size of the window used for berycentre coalculation
                        [eta,phi]
  --windowDupl WINDOWDUPL WINDOWDUPL
                        Size of the window used for duplicate removal
                        [eta,phi]
  --dEta DETA DETA      Size of the tower in eta
  --maxEta MAXETA       Maximum eta
  --dPhi DPHI           Size of the tower in phi
  --numPhi NUMPHI       Number of the towers in phi
  --zoom ZOOM ZOOM      How many bins around centre should be visible
  --event EVENT         Number of an event to draw
~~~

## Reconstruction monitor: `SingleParticleRecoMonitors`

Monitor plots of the reconstruction procedure.
The plots are created based on the clusters and generated particles (for single particle events only).

~~~{.sh}
python scripts/plot_recoMonitor.py <INPUT_FILE_NAME>.root ENERGY <INPUT_FILE_WITH_MC_INFORMATION>.root
~~~

The plots are saved to root file and as a png. The first row includes the energy distribution, energy distribution as a function of phi, number of reconstructed clusters and the number of reconstructed clusters as a function of phi.
The second row contains information about any reconstructed duplicates (difference of energy, eta, phi and R=sqrt(eta^2+phi^2) with respect to the most energetic cluster).
The last row contains comparison to the MC particle: difference in eta distribution, difference in eta as a function of eta, difference in phi distribution and difference in phi as a function of phi.


## Reconstruction monitor with upstream energy correction

python scripts/plot_recoMonitor.py root://eospublic.cern.ch//eos/fcc/users/n/novaj/combCaloForBerlin/output_combCalo_reconstructionSW_e50GeV_part1.root 50 root://eospublic.cern.ch//eos/fcc/users/n/novaj/combCaloForBerlin/output_combCalo_e50GeV_part1_v3.root --clusterColl EcalClusters --correctionParams 0.1369 0.004587 0.1692 0.6769 --cellColl ECalCellsForSW --bitfield system:4,cryo:1,type:3,subtype:3,cell:6,eta:9,phi:10

--- Everything in reconstruction file (GenParticles, GenVertices, Cells, Clusters), simulations for Berlin
**** Upstream material 1.2 X0, B = 4 T
python scripts/plot_clusterAnalysis.py root://eospublic.cern.ch//eos/fcc/users/n/novaj/combCaloForBerlin/output_combCalo_reconstructionSW_e?GeV_bfield1_eta0_v2.root 20 50 100 500 1000 --clusterColl CombinedClusters --correctionParams 0.1037 0.0007507 0.1382 1.002 --cellColl ECalCellsForSW --bitfield system:4,cryo:1,type:3,subtype:3,cell:6,eta:9,phi:10 -o combCalo_?GeV_bfield1_eta0 -r energy
**** Upstream material 1.2 X0, B = 0 T
python scripts/plot_clusterAnalysis.py root://eospublic.cern.ch//eos/fcc/users/n/novaj/combCaloForBerlin/output_combCalo_reconstructionSW_e?GeV_bfield0_eta0_v2.root 20 50 100 500 1000 --clusterColl CombinedClusters --correctionParams 0.05774 0.000315 0.0621 0.8285 --cellColl ECalCellsForSW --bitfield system:4,cryo:1,type:3,subtype:3,cell:6,eta:9,phi:10 -o combCalo_?GeV_bfield0_eta0 -r energy


List of additional options:
~~~{.sh}
  -r REGEX [REGEX ...], --regex REGEX [REGEX ...]
                        String to insert in place of '?' character in file
                        names ("energy" inserts the values of energies)
  -o OUTPUT, --output OUTPUT
                        Output file name
  -v, --verbose         Verbose
  --particleColl PARTICLECOLL
                        Name of the MC particle collection
                        (fcc::MCParticleCollection)
  --clusterColl CLUSTERCOLL
                        Name of the clusters collection
                        (fcc::CaloClusterCollection)
  --dEta DETA DETA      Size of the tower in eta
  --maxEta MAXETA       Maximum eta
  --dPhi DPHI           Size of the tower in phi
  --numPhi NUMPHI       Number of the towers in phi
~~~


## Energy resolution

Input files: ROOT files with energy distribution saved under name "energy" (can be changed using option `-n`). Each ROOT file should contain distribution for one energy. Such input files may be obtained e.g. by [`plot_recoMonitor.py`](scripts/plot_recoMonitor.py) macro.
Energy distributions are fitted twice with Gaussian, and the energy resolution plot is fitted with p0+p1/sqrt(E) function.

~~~{.sh}
python scripts/plot_enResolution.py energy?GeV.root 20 50 100 200 500 1000 -r energy
~~~
python scripts/plot_enResolution.py combCalo_?GeV_bfield1_eta0.root 20 50 100 500 -n 'energyCorrected' -o resolution_bfield1 -r energy
-n: name of the histogram

## Energy resolution - comparison

Input files: ROOT files produced by `plot_enResolution.py`.

~~~{.sh}
python scripts/plot_compareResolution.py '?/energy_resolution_plots.root' 0 -r noBfield withBfield
~~~

Will produce a plot comparing two resolutions, from files

~~~{.sh}
./noBfield/energy_resolution_plots.root
./withBfield/energy_resolution_plots.root
~~~

Example of use of additional parameters:

~~~{.sh}
python scripts/plot_compareResolution.py '?/energy_resolution_plots.root' 0 -r noBfield withBfield -t 'Histogram title' -m 0.1 -l '?' --sequentialColours
python scripts/plot_compareResolution.py '?/energy_resolution_plots.root' 0 -r noBfield withBfield -t 'Histogram title' -m 0.1 -l 'no B field: formula' 'B = 4 T: formula'
~~~

`formula` inside the legend will be replaced with energy resolution formula with fit parameters.



## Sampling fraction

Input files: ROOT files with sampling fraction histograms, as generated in [`FCCSW`](https://github.com/HEP-FCC/FCCSW) using [`SamplingFractionInLayers` algorithm](https://github.com/HEP-FCC/FCCSW/blob/master/Detector/DetStudies/src/components/SamplingFractionInLayers.h).

If simulation was performed initially for more layers than the number of layers to study, one can use `--merge` option, specifying the number of consecutive layers to be merged into one layer. E.g. if simulation was run for 10 layers, using `--merge 4 6` means that he sampling fraction is ploted for 2 layers (layers 1-4 as the first layer, and layers 5-10 as the second one).

### Example
Using the output of an [`example from FCCSW`](https://github.com/HEP-FCC/FCCSW/blob/master/Detector/DetStudies/tests/options/samplingFraction_inclinedEcal.py), one may generate the sampling fraction for the calibration to EM scale.

~~~{.sh}
python scripts/plot_samplingFraction.py histSF_inclined_e50GeV_eta0_1events.root 0
~~~
-> to extract sampling fraction for every existing layer in the geometry;


~~~{.sh}
python scripts/plot_samplingFraction.py histSF_inclined_e50GeV_eta0_1events.root 0 --merge 4 4 4 4 4 4 4 4
~~~
-> to extract sampling fraction for 8 layers, created by merging 4 adjacent layers into one;

The output ROOT file, as well as the printout of the macro contain the values of the sampling fractions for detector layers.
The printout `samplingFraction = [ ... ]` may be directly copied to the option of the [`CalibrateInLayers` algorithm](https://github.com/HEP-FCC/FCCSW/blob/master/Reconstruction/RecCalorimeter/src/components/CalibrateInLayersTool.h), as in the [example](https://github.com/HEP-FCC/FCCSW/blob/master/Reconstruction/RecCalorimeter/tests/options/runEcalInclinedDigitisation.py#L28).


## Correction for the upstream material

Input: ROOT file with histograms of the energy in the first layer of the calorimeter (X-axis) and the correspoding energy deposited in the upstream material.
For sufficiently thin layer it is a linear dependence.
The parameters for the upstream correction depend on both energy and direction of the particle (pseudorapidity). It is assumed there is no azimuthal angle dependence (which can be checked on the control plot, if `--preview` option is used). `--preview` option also plots the dependence of the upstream energy on the energy deposited in the first layer.

~~~{.sh}
python scripts/plot_upstreamCorrecton.py histUpstreamCorrection_e?GeV_Bfield1.root 100 -r energy -o Bfield1  --preview
~~~

Plots a preview plot.

~~~{.sh}
python scripts/plot_upstreamCorrecton.py histUpstreamCorrection_e?GeV_Bfield1.root 20 50 100 200 500 -r energy -o Bfield1  --preview
~~~

Plots the energy dependence of parameters. It is assumed that the constant parameter increases linearily with the energy. The slope parameter is fitted with [0]+[1]/sqrt(E) function.

~~~{.sh}
python scripts/plot_upstreamCorrecton.py histUpstreamCorrection_e?GeV_Bfield1.root 20 50 100 200 500 -r energy --etaValues 0 0.25 0.5 0.75 1.0 1.25 1.5 -o Bfield1  --preview
~~~

Plots the pseudorapidity dependence of the parameters describing the energy-dependance.


# How to create own analysis

1. Create Class deriving from `BaseAnalysis` or `BaseTwoFileAnalysis`. Include:
  - histograms (push back to m_histograms)
  - processEvent method
  - finishLoop method
2. Add class to `include/LinkDef.h` and `CMakeLists.txt` so it can be accessed from ROOT.
3. Create analysis script:
  - add command line arguments
  - retrieve histograms after analysis is done

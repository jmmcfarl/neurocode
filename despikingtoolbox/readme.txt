
This set of files implement the Bayesian LFP despiking algorithm derived 
in:

T Zanos, PJ Mineault and CC Pack (2010) Removal of spurious correlations
between spikes and local field potentials, J Neurophysiol (in press). 

To get started, add the toolbox folder to your Matlab path and then open 
ExampleDespiking.m in Cell Mode in Matlab and read through it. A data set
is furnished (examplelfpdespiking.mat). The files in the toolbox folder 
are copiously documented; type help _filename_ at the command line for 
usage information. A copy of the appendix of the paper containing the 
derivation of the algorithm and its implementation is included as 
paperappendix.pdf.

Other files outside of the toolbox directory are support files for the 
example; they are not required by the core methods.

The files in the toolbox folder are:

despikeLFP.m : Main despiking routine. 
despikeLFPbyChunks.m : Does despiking by chunks for signals that are too 
large for normal despiking. Calls despikeLFP.
fitLFPpowerSpectrum.m: Can be used to find a good value for g, the 
parameter that determines the prior on the LFP. 
ksr.m: Kernel smoothing regression, by Yi Cao, File ID: #19195 on Matlab 
Central, required by fitLFPpowerSpectrum

These files are distributed under the GPL (see license.txt), save for 
ksr.m which comes with a BSD license. 

History:

Version 1.02: 04/11/2010 - Paper accepted, first public version
Version 1.01: 25/08/2010 - Fixed bug in despikeLFP when zero spikes are present
Version 1.0 : 26/06/2010 - Original version

Author:

Patrick Mineault
patrick DOT mineault AT gmail DOT com

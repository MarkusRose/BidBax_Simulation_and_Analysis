# Bcl2Analysis

This repository contains the scripts used to simulate and analyze
two-color-channel confocal images of diffraction-limited spots. It is optimized
for the analysis of Bcl-2 family proteins on a planar supported lipid bilayer. The repository contains two packages:
 - ConfocalSimulation
 - IJPlugins

## Simulation package: ConfocalSimulation

The simulation is written in python.
All required simulation funcions are in the file "confocalsim.py".
"benchmarking.py" demonstates the use of the confocal simulation functions.

The principle of the simulation is demonstrated in the jupyter notebook "ConfocalSimulation_ThePrinciple.ipynb".

## Analysis package: IJPlugins

The image analysis is implemented using plugnis for ImageJ, written in Java.
 - The particle fitting is carried out in both "sixFittingMethod1.java" and
 - "sixFittingMethod2.java". The methods are called in "SixFitPlugin_v11.java" and applied to a confocal two-channel image pair loaded in ImageJ. 
 - With the wrapper plugin "SixFitFolderAnalysis_.java" all images contained in a selected folder are analyzed.



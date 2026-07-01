## MR-SHIFT

MR-SHIFT is a MATLAB App Designer application for simulating, fitting, and analysing NMR spectra using spin-system information stored in a user-defined metabolite database.

The app is designed to work with Bruker TopSpin data and supports simulation and fitting workflows for 1D spectra, T1/T2 relaxation analyses, and selected 2D simulations.

## Main files
MR_SHIFT_1p3.mlapp
Main MATLAB App Designer application. This is the file to open and run in MATLAB.
getTopSpinAddLB.m
Reads Bruker TopSpin data and extracts the experimental and processing parameters required by MR-SHIFT. It also handles line broadening settings.
database_phosphoromics_P31only_CREATE.m
Example database creation script. Users can edit this file, or use it as a template, to define their own metabolite database.
NSpin_Simulation_1p3.m
Simulation engine for quantum-mechanical spin-system calculations. This function calculates spectra by matrix diagonalization and generates the spectral skeleton used for subsequent fast simulations and fitting.
FastFIDgeneration_1p3.m
Fast spectrum-generation engine based on the skeleton produced by NSpin_Simulation_1p3.m. It is used for replotting and for amplitude/T2 fitting.
Important note for users

Most users should only modify the database creation script, for example:

database_phosphoromics_P31only_CREATE.m

The following files are core calculation engines and should not normally be changed:

NSpin_Simulation_1p3.m
FastFIDgeneration_1p3.m

Modifying these files may change the numerical behaviour of the simulations and fits.

## Requirements

This app was tested using MATLAB R2024b.

The base app requires:

MATLAB

Additional requirements for specific features:

T1/T2 fitting
Curve Fitting Toolbox
Statistics and Machine Learning Toolbox
2D simulations
SPINACH; please refer to the SPINACH documentation for its requirements.

The listed MathWorks products were identified using:

matlab.codetools.requiredFilesAndProducts
Installation
Download or clone this repository.
Open MATLAB.
Set the repository folder as the current MATLAB working directory.
Make sure all .m files are on the MATLAB path.
Open the app file:
MR_SHIFT_1p3.mlapp

Alternatively, the app can be opened from the MATLAB file browser by double-clicking MR_SHIFT_1p3.mlapp.

## Input data

MR-SHIFT is intended to read processed Bruker TopSpin data.

The TopSpin data folder should contain the standard Bruker experiment structure.

The function getTopSpinAddLB.m is responsible for reading the experimental and processing parameters from the TopSpin files.

## Database format

Metabolites and spin systems are defined in a MATLAB structure called simDB.

The example script:

database_phosphoromics_P31only_CREATE.m

creates a database for 31P metabolites. Each database entry contains information such as chemical shifts, uncertainty estimates, T2 values, and notes.

Example structure:

simDB.MetaboliteName = struct( ...
    'w_CS_J', chemical_shift_or_spin_matrix, ...
    'CSerr', chemical_shift_error, ...
    'T2', T2_value, ...
    'T2err', T2_error, ...
    'notes', 'description' ...
);

Chemical shifts are given in ppm and J-couplings are given in Hz.

Users can create their own database by editing the database creation script and saving the resulting .mat database file.

## Basic workflow
Prepare or load a metabolite database.
Open MR_SHIFT_1p3.mlapp.
Select the TopSpin data folder and experiment/processing numbers.
Plot the experimental spectrum.
Add metabolites from the database.
Simulate the selected spin systems.
Fit amplitudes, chemical shifts, T2 values, or other selected parameters as required.
Save the simulation or fit results.

## Example
This repository includes an example TopSpin-style data folder named: 1. with a 1D 31P experiment.

To open the example dataset in MR-SHIFT, use the following settings in the app:
Path: .
Experiment number: 1 

The corresponding fitresult is also included: fitresult_1.m

## Line broadening

MR-SHIFT supports additional line broadening through the TopSpin data-reading function.

The line broadening mode can be interpreted in two ways:

extra: add additional exponential line broadening on top of the TopSpin processing value.
total: use the specified value as the total exponential line broadening.

## Outputs

Depending on the workflow, MR-SHIFT can generate and save simulated spectra, fitted parameters, fit results, and database-derived simulations.

The exact output depends on the selected fitting or simulation mode.

## Citation

If you use MR-SHIFT in published work, please cite the associated publication:

[Add citation once available]

## License

[Add license information here.]

## Contact

For questions, bug reports, or suggestions, please contact:

[Add contact information here.]

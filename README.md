# MR-SHIFT

MR-SHIFT is a MATLAB App Designer application for simulating, fitting, and analysing NMR spectra using spin-system information stored in user-defined metabolite databases.

The app is designed to work with processed Bruker TopSpin data and supports simulation and fitting workflows for 1D spectra, T1/T2 relaxation analyses, lineshape distortions, and selected 2D simulations.

## Stable publication version

The validated version of MR-SHIFT associated with the publication is **v1.3.0**.

For reproducibility of the analyses reported in the manuscript, please use the [v1.3.0 release](../../releases/tag/v1.3.0).

The `main` branch may contain newer features, improvements, and ongoing developments that have not yet undergone the same level of validation.

## Repository structure

```text
MR_SHIFT/
│
├── MR_SHIFT_1p3.mlapp
├── createAutoFitRules.m
│
├── core/
│   ├── NSpin_Simulation_1p3.m
│   ├── FastFIDgeneration_1p3.m
│   └── getTopSpinAddLB.m
│
├── databases/
│   ├── database_phosphoromics_P31only.mat
│   ├── database_phosphoromics_P31only_CREATE.m
│   ├── database_phosphoromics_1H31P.mat
│   ├── database_phosphoromics_1H31P_CREATE.m
│   └── database_1H31P.xlsx
│
├── lineshape/
│   └── functions for lineshape distortion modelling
│
├── simulations_2d/
│   └── SPINACH-based 2D simulation functions
│
├── 1/
│   └── bundled example Bruker dataset
│
├── fitresult_1.mat
├── MR_SHIFT.prj
├── LICENSE
└── README.md
```

## Main user files

### `MR_SHIFT_1p3.mlapp`

Main MATLAB App Designer application. This is the file most users should open to run MR-SHIFT.

### `createAutoFitRules.m`

User-editable configuration file for automated fitting.

This file defines which experiments should be fitted and can be modified to control which metabolites are included in a full fitting workflow.

### `databases/`

Contains the metabolite databases used by MR-SHIFT together with the scripts used to create or modify them.

The database creation scripts can also be used as templates for building custom metabolite databases.

## Core calculation functions

The main simulation engine is located in `core/`.

### `core/NSpin_Simulation_1p3.m`

Quantum-mechanical spin-system simulation engine.

This function calculates spectra by matrix diagonalization and generates the spectral skeleton used for subsequent fast simulations and fitting.

### `core/FastFIDgeneration_1p3.m`

Fast spectrum-generation engine based on the spectral skeleton produced by `NSpin_Simulation_1p3.m`.

It is used extensively during spectrum generation and fitting.

### `core/getTopSpinAddLB.m`

Reads processed Bruker TopSpin data and extracts the experimental and processing parameters required by MR-SHIFT.

It also handles line-broadening settings.

These core functions should not normally need to be modified by users. Changes to them may alter the numerical behaviour of simulations and fits.

## Requirements

MR-SHIFT was tested using **MATLAB R2024b**.

The base application requires:

- MATLAB

Additional requirements depend on the functionality being used.

### T1/T2 fitting

- Curve Fitting Toolbox
- Statistics and Machine Learning Toolbox

### 2D simulations

Selected 2D simulations use **SPINACH**.

Please refer to the SPINACH documentation for its installation and requirements.

The required MathWorks products were identified using:

```matlab
matlab.codetools.requiredFilesAndProducts
```

## Installation

1. Download or clone this repository.
2. Open MATLAB.
3. Open `MR_SHIFT_1p3.mlapp`.
4. Run the app.

The application automatically adds its internal module folders to the MATLAB path when it starts.

It is therefore not necessary to manually add `core/`, `databases/`, `lineshape/`, or `simulations_2d/` to the MATLAB path before running the app.

Developers may alternatively open `MR_SHIFT.prj` to work with MR-SHIFT as a MATLAB Project.

## Input data

MR-SHIFT is intended to read processed Bruker TopSpin data.

The TopSpin data directory should contain the standard Bruker experiment and processing structure.

The function:

```text
core/getTopSpinAddLB.m
```

is responsible for reading the required experimental and processing parameters.

## Metabolite databases

Metabolites and spin systems are defined using MATLAB structures.

The supplied database creation scripts are located in:

```text
databases/
```

For example:

```text
databases/database_phosphoromics_P31only_CREATE.m
```

creates the supplied ^31P metabolite database.

A typical database entry contains information such as:

```matlab
simDB.MetaboliteName = struct( ...
    'w_CS_J', chemical_shift_or_spin_matrix, ...
    'CSerr', chemical_shift_error, ...
    'T2', T2_value, ...
    'T2err', T2_error, ...
    'notes', 'description' ...
);
```

Chemical shifts are expressed in **ppm** and scalar J-couplings in **Hz**.

The `w_CS_J` field can contain either a single chemical shift or a spin-system topology matrix containing chemical shifts and J-couplings.

Users can create their own databases by modifying one of the supplied database creation scripts.

The generated `.mat` database files are stored together with their corresponding creation scripts in the `databases/` directory.

## Automated fitting

Automated fitting is configured through:

```text
createAutoFitRules.m
```

This script is intentionally kept in the repository root because it is designed to be modified by users.

It can be used to define:

- which experiments should be fitted;
- which metabolites should be included;
- experiment-specific fitting behaviour.

## Basic workflow

A typical MR-SHIFT workflow is:

1. Open `MR_SHIFT_1p3.mlapp`.
2. Select a metabolite database.
3. Select the Bruker TopSpin data directory and experiment/processing numbers.
4. Plot the experimental spectrum.
5. Add metabolites from the database.
6. Simulate the selected spin systems.
7. Fit amplitudes, chemical shifts, T2 values, or other selected parameters as required.
8. Save the resulting simulation or fit.

## Example dataset

The repository includes a bundled TopSpin-style example dataset in:

```text
1/
```

This contains a 1D ^31P experiment and allows MR-SHIFT to display an example spectrum immediately after installation.

The default settings are:

```text
Path: .
Experiment number: 1
Processing number: 1
```

The corresponding example fit result is also included:

```text
fitresult_1.mat
```

## Lineshape distortion modelling

Additional lineshape functionality is contained in:

```text
lineshape/
```

This module supports modelling of non-ideal spectral lineshapes, including shim-related distortions and asymmetric lineshape effects.

These features are optional and do not modify the underlying core spin-system simulation engine.

## 2D simulations

Functions used for selected 2D NMR simulations are located in:

```text
simulations_2d/
```

These simulations currently use SPINACH and are kept separate from the core 1D simulation engine so that additional 2D experiments can be incorporated independently.

## Line broadening

MR-SHIFT supports additional exponential line broadening through the TopSpin data-reading functionality.

The line-broadening mode can be interpreted in two ways:

- `extra` — adds additional exponential line broadening on top of the TopSpin processing value.
- `total` — uses the specified value as the total exponential line broadening.

## Outputs

Depending on the selected workflow, MR-SHIFT can generate and save:

- simulated spectra;
- fitted parameters;
- fit results;
- relaxation-analysis results;
- database-derived simulations.

The exact output depends on the selected fitting or simulation mode.

## Citation

If you use MR-SHIFT in published work, please cite the associated publication:

**[Add final publication citation once available]**

For exact reproducibility of the published implementation, use **MR-SHIFT v1.3.0**.

## License

This software is available for academic and non-commercial research use.

Commercial use requires prior written authorization from CIC bioGUNE.

See the [LICENSE](LICENSE) file for details.

## Contact

For questions, bug reports, or suggestions, please contact:

djardon@cicbiogune.es

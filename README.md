# RETSE_Metal
This repository contains the code and datasets for "Distortion correction in TSE near titanium implants at 0.55T using reversed frequency-encoding and model-based reconstruction", by Bochao Li, Nam G. Lee, Daehyun Yoon, Kübra Keskin, Alexander R. Toews, Jay Acharya, Jordan S. Gross, Brian A. Hargreaves, Krishna S. Nayak.

Bochao Li, University of Southern California, July 2025.
## Distortion Correction Pipeline for Cartesian MRI Data

This repository contains MATLAB code to correct gradient non-linearity (GNL) and off-resonance distortion in Cartesian MRI acquisitions using dual frequency-encoding techniques. It includes tools for image reconstruction, off-resonance displacement estimation, and phase-constrained SENSE based joint reconstruction.

## Setup

### 1. Required Installation

- **[BART (Berkeley Advanced Reconstruction Toolbox)](https://mrirecon.github.io/bart/)**  

- **[ISMRMRD (MRI raw data format)](https://ismrmrd.readthedocs.io/en/latest/)**  

### 2. Required third-party dependencies
Before running the pipeline, make sure the following third-party dependencies are downloaded and placed under the directory: **/Helper_Functions/thirdparty/** (Please create it for first use)
1. **BART (MATLAB interface)**  
   https://github.com/mrirecon/bart/tree/master/matlab

2. **FINUFFT**  
   https://github.com/flatironinstitute/finufft

3. **mapVBVD**  
   https://github.com/pehses/mapVBVD

4. **export_fig**  
   https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig

## Folder Structure

- `main_recon.m` — Main entry point for the reconstruction pipeline.
- `Helper_Functions/` — Contains custom MATLAB functions and third-party dependencies.
- `results/` — Stores output images and figures (not included by default).
- `data/` — Directory for raw TWIX data (not included).

## Data

To run this pipeline, you need reverse-encoding Turbo Spin Echo (TSE) raw MRI data acquired in two frequency-encoding directions data.

A sample dataset (THA, rBW = 150 Hz/Px) is publicly available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15833570.svg)](https://doi.org/10.5281/zenodo.15833570)

### Folder Setup

After downloading the dataset, organize it into the following directory structure:
```
\raw\
  |- PE_RL\
    |- <PE_RL_data>.dat
  |- PE_LR\
    |- <PE_LR_data>.dat
```
- Place the `.dat` file from each acquisition direction in the respective folder.

> **Note:** These `.dat` files are Siemens raw TWIX files. You will convert them to ISMRMRD format as part of the pipeline using the provided `auto_convert_twix_to_ismrmrd()` function.


## How to Run

Before running `main_recon.m`, please update the following variables at the top of the script:

```matlab
ismrmrd_path = 'path_of_your_ismrmrd_installation';
bart_path    = 'path_of_your_bart_installation';
```

## Output
After ```show_results```, ```img_recon``` includes Dir.1 TSE, corrected Dir.1 TSE, Dir.2 TSE, corrected Dir.2 TSE, corrected image from PES-SENSE
# RETSE_Metal
This repository contains the code and datasets for "ADistortion correction in TSE near titanium implants at 0.55T using reversed frequency-encoding and model-based reconstruction", by Bochao Li1, Nam G. Lee, Daehyun Yoon, Kübra Keskin, Alexander R. Toews, Jay Acharya, Jordan S. Gross, Brian A. Hargreaves, Krishna S. Nayak.

Bahadir Alp Barlas, University of Southern California, May 2025.
## Distortion Correction Pipeline for Cartesian MRI Data

This repository contains MATLAB code to correct gradient non-linearity (GNL) and off-resonance distortion in Cartesian MRI acquisitions using dual frequency-encoding techniques. It includes tools for image reconstruction, B0 field estimation using the Mullen method, and PEC-SENSE based joint reconstruction.

## Prerequisites

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

## How to Run

Before running `main_recon.m`, please update the following variables at the top of the script:

```matlab
ismrmrd_path = 'path_of_your_ismrmrd_installation';
bart_path    = 'path_of_your_bart_installation';
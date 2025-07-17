# RETSE_Metal
This repository contains the code and datasets for "Distortion correction in TSE near titanium implants at 0.55T using reversed frequency-encoding and model-based reconstruction", by Bochao Li, Nam G. Lee, Daehyun Yoon, Kubra Keskin, Alexander R. Toews, Jay Acharya, Jordan S. Gross, Brian A. Hargreaves, Krishna S. Nayak.

[Bochao Li](mailto:bochaoli@usc.edu), University of Southern California, July 2025.
## Distortion Correction Pipeline for Cartesian MRI Data

This repository contains MATLAB code to correct gradient non-linearity (GNL) and off-resonance distortion in Cartesian MRI acquisitions using dual frequency-encoding techniques. It includes tools for image reconstruction, [off-resonance displacement estimation](https://data.mendeley.com/datasets/2nbpddxd8f/1), and phase-constrained SENSE based joint reconstruction.

## Setup

### 1. git clone https://github.com/usc-mrel/RETSE_Metal.git

### 2. Required Installation

- **[BART (Berkeley Advanced Reconstruction Toolbox)](https://mrirecon.github.io/bart/)**  

- **[ISMRMRD (MRI raw data format)](https://ismrmrd.readthedocs.io/en/latest/)**  

### 3. Required third-party dependencies
Before running the pipeline, make sure the following third-party dependencies are downloaded and placed under the directory: **/Helper_Functions/thirdparty/** (Please create it for first use)
1. **BART (MATLAB interface)**  
   https://github.com/mrirecon/bart/tree/master/matlab

2. **FINUFFT**  
   https://github.com/flatironinstitute/finufft
- Please install it correctly according to the [online document](https://finufft.readthedocs.io/en/latest/install.html#install) on your system.

3. **mapVBVD**  
   https://github.com/pehses/mapVBVD

4. **export_fig**  
   https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig

## Folder Structure

- `main_recon.m`: Main entry point for the reconstruction pipeline with automatic procedures.
- `Helper_Functions/`: Contains custom MATLAB functions and third-party dependencies.
- `results/`: Stores output images and figures (not included by default).
- `raw/`: Directory for raw data (available on Zenodo).

## Data

To run this pipeline, you need reverse-encoding Turbo Spin Echo (TSE) raw MRI data acquired in two frequency-encoding directions data.

A sample dataset of ISMRMRD format (THA, rBW = 150 Hz/Px) is publicly available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14752389.svg)](https://doi.org/10.5281/zenodo.14752389)

### Folder Setup

After downloading the dataset, please unzip it and copy to the code folder, then you can get the following directory structure manually:
```
\raw\
  |- PE_RL\
    |- <PE_RL_data>.dat
  |- PE_LR\
    |- <PE_LR_data>.dat
```

## How to Run

Before running `main_recon.m`, please update the following variables at the top of the script:

```matlab
ismrmrd_path = 'path_of_your_ismrmrd_installation'; (not BART MATLAB codes under \Helpfer_Functions)
bart_path    = 'path_of_your_bart_installation';
```

## Outputs
- `.json` are the configuration files that are generated in the pipeline, saved in `\raw` folder.

- Estimated off-resonance displacement maps in meter are saved under `offres` .

- Corrected image per slice using proposed phase-constrained SENSE method are saved (`img_typ1_xxx.cfl`) in `\result\dualrecon_pcs` folder that is generated automatically; `\result\meas_xxx` folder contains non-corrected, corrected image for each direction's data. 

- After ```show_results```, ```img_recon``` has 5 multi-slice TSE images, including Dir.1 TSE, corrected Dir.1 TSE, Dir.2 TSE, corrected Dir.2 TSE, corrected image from PES-SENSE.

## References and Acknowledgements
1. Mullen M, Garwood M. Dual polarity encoded MRI using high bandwidth radiofrequency pulses for robust imaging with large field inhomogeneity. Magn Reson Med. 2021;86(3):1271-1283. doi:10.1002/mrm.28771
2. Tao S, Trzasko JD, Shu Y, Huston J 3rd, Bernstein MA. Integrated image reconstruction and gradient nonlinearity correction. Magn Reson Med. 2015;74(4):1019-1031. doi:10.1002/mrm.25487
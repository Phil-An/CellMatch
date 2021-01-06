
# Cellmatch – Matlab

A computational toolbox for identifying neurons, recorded with *in vivo* calcium imaging using miniaturized epifluorescence microscopes (miniscopes), *post hoc* in histology. Identifying *in vivo* recorded neurons *post hoc* allows a cell-type specific discrimination of neuronal subtypes. *in vivo* recorded calcium imaging data was pre-processed using Inscopix Mosaic software. The best way to start is by looking at the example data provided.

![Graphical abstract - cellmatch](https://ars.els-cdn.com/content/image/1-s2.0-S0165027020301886-ga1_lrg.jpg )

# Features
* Semi-automatic alignment and tracking of the same neurons recorded across multiple *in vivo* Ca2+ imaging experiments 
* 3D histology reconstruction from laser confocal microscopy scans
* Manual cell segmentation tool in histology images
* Automatic matching of *in vivo* recorded cells, *post hoc* in histology images

# Citation and description of the method

If you use this code please cite the paper:

Anner, P., Passecker, J., Klausberger, T., & Dorffner, G. (2020). Ca2+ imaging of neurons in freely moving rats with automatic *post hoc* histological identification. Journal of Neuroscience Methods, 341(January), 108765. [https://doi.org/10.1016/j.jneumeth.2020.108765](https://doi.org/10.1016/j.jneumeth.2020.108765)


# Dependencies
The following Matlab toolboxes are required:
1.	Statistics and Machine Learning Toolbox
2.	Image processing toolbox
3.	Optimization toolbox

The following Matlab implementations are required from the Mathworks File Exchange repository:
1. [Hessian based Frangi Vesselness filter](https://de.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter) by Dirk-Jan Kroon
2. [ B-spline Grid, Image and Point based Registration](https://de.mathworks.com/matlabcentral/fileexchange/20057-b-spline-grid-image-and-point-based-registration?s_tid=prof_contriblnk) by Dirk-Jan Kroon
3. [Region Growing (2D/3D grayscale)](https://de.mathworks.com/matlabcentral/fileexchange/32532-region-growing-2d-3d-grayscale?focused=5195969&s_tid=gn_loc_drop&tab=function) by Daniel Kellner
4. [sort_nat: Natural Order Sort](https://de.mathworks.com/matlabcentral/fileexchange/10959-sort_nat-natural-order-sort) by Douglas Schwarz

# Documentation
The main script is [Main.m](https://github.com/Phil-An/CellMatch/blob/main/Main.m) that will guide you through the processing operations. [Sampledata.mat](https://github.com/Phil-An/CellMatch/blob/main/Sampledata.mat) stores manually derived parameters for processing the [sample data](https://github.com/Phil-An/CellMatch/tree/main/Data) provided.


## Processing *in vivo* recorded calcium imaging data
The script [InVivo_Align_CaImagingSessions.m](https://github.com/Phil-An/CellMatch/blob/main/InVivo_Align_CaImagingSessions.m) loads *in vivo* recorded Ca2+ imaging data. Ca2+ imaging data was pre-processed with [Inscopix](https://www.inscopix.com/) Mosaic software. For further details on required structure of *in vivo* recorded Ca2+ imaging data, please refer to section [*in vivo* Ca2+ imaging data ](#ca2+data). First, the final experiment for the detection of blood vessels (RECDSA) is loaded and processed. All following experiments are loaded and co-registered with RECDSA. [Sampledata.mat](https://github.com/Phil-An/CellMatch/blob/main/Sampledata.mat) data includes all manually defined control points for registration. 

## Processing of histology data
[Volumetric laser confocal microscope scans](https://github.com/Phil-An/CellMatch/tree/main/Data/Post%20hoc%20histology) are provided within the sample data. Histology was not labeled with antibodies and endogenous Green fluorescence protein (GFP) signals from the Calciun indicator GCaMP6f were scanned. Histology imaging data is processed within [Main.m](https://github.com/Phil-An/CellMatch/blob/main/Main.m). 

If you performed co-labeling of neurons in histology, you must export z-stacks of all channels separately and load them within [Main.m](https://github.com/Phil-An/CellMatch/blob/main/Main.m). All transformations must then also be applied to the data structures of histology channel 2. When you assess the expression of fluorescent markers of matched cells using [plot_interactive_identifiedcells.m](https://github.com/Phil-An/CellMatch/blob/main/functions/plot_interactive_identifiedcells.m), you can switch the displayed channels of the histology scans. 

## <a id="ca2+data"></a> Structure of *in vivo* Ca2+ imaging data 
In the sample data, three individual *in vivo* Ca2+ imaging recordings are provided:
1. [2411](https://github.com/Phil-An/CellMatch/tree/main/Data/In%20vivo%20Ca2%2B/2411)
2. [2811](https://github.com/Phil-An/CellMatch/tree/main/Data/In%20vivo%20Ca2%2B/2811)
3. [1412](https://github.com/Phil-An/CellMatch/tree/main/Data/In%20vivo%20Ca2%2B/1412)

Each *in vivo* Ca2+ recording consists of:
1. a Minimum intensity projection image of the original movie ([*.minrec.tif](https://github.com/Phil-An/CellMatch/blob/main/Data/In%20vivo%20Ca2%2B/1412/JP68-1412-01-minrec.tif))
2. a [.csv](https://github.com/Phil-An/CellMatch/blob/main/Data/In%20vivo%20Ca2%2B/1412/JP68-1412-01-traces.csv) file containing calcium transients and calcium events
3. and a folder [IC](https://github.com/Phil-An/CellMatch/tree/main/Data/In%20vivo%20Ca2%2B/1412/IC) containing all independent component images (IC) of all identified neurons.

Each [.csv](https://github.com/Phil-An/CellMatch/blob/main/Data/In%20vivo%20Ca2%2B/1412/JP68-1412-01-traces.csv) file must follow a specific structure:
1. The first row contains variable names (export them from Mosaic)
2. The first column includes a time code
3. The following columns contain dF/F traces, optionally calcium event traces as well. 

Further, a final [experiment](https://github.com/Phil-An/CellMatch/tree/main/Data/In%20vivo%20Ca2%2B/1912-DSA) a fluorescent dye was injected into the tail vein of the anaesthetized animal. Prior to the injection of the fluorescent dye, a three-minute-long recording was obtained. Following the injection of the fluorescent dye, another three-minute-long recording was obtained. In this experiment no neuronal activity was recorded. Instead, projection images of the original movies are provided:
1. [min-natice.tif](https://github.com/Phil-An/CellMatch/blob/main/Data/In%20vivo%20Ca2%2B/1912-DSA/JP68-1912-min_native.tif) is a minimum intensity projection image of the recording acquired prior to the injection of the fluorescent dye.
2. [mean-native.tif](https://github.com/Phil-An/CellMatch/blob/main/Data/In%20vivo%20Ca2%2B/1912-DSA/JP68-1912-mean-native.tif) is a mean intensity projection image of the recording acquired prior to the injection of the fluorescent dye.
3. [maxIP-CE.tif](https://github.com/Phil-An/CellMatch/blob/main/Data/In%20vivo%20Ca2%2B/1912-DSA/JP68-1912-MaxIP_CE.tif) is a standard deviation projection image of the recording acquired following the injection of the fluorescent dye.

# License 
This program is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) license. Non-commercial use of the software is free.

This program is distributed in the hope that it will be useful for research, but without any warranties of title, merchantability, fitness for a particular purpose, non-infringement, absence of latent or other defects, accuracy, or the presence or absence of errors, whether or not known or discoverable. See the CC BY-NC-ND License for more details.
[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

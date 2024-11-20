# FOV_alignment

Semi-automated process to determine whether spines from two time-points were retained, lost, or added by comparing distances between dendritic spine coordinates projected onto the dendritic branch and marked fiducial points (i.e branch points, stable spines) 

## Installation Instructions
Analysis run on MATLAB2021b, no specific software installation required besides downloading MATLAB (any version should suffice)

## Usage Instructions

### Overview
All code can be run through "scripts/run_FOV_alignment_by_dendriticDistance.m"

#### Setting directories/paths
Before running the code, be sure to change the following paths to your home directory and location of code: 
- homepath = '/Users/ktsimring/Documents/';
- path = '/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/2P_Analysis';


#### Input Files: data is in "demo_data/2P_data" folder
there should be at least two timepoints in order to compare the spine identities, with the following files:
1. mean_intensity_projection.tif, which this is the average intensity projection of a dendritic branch 
2. RoiSet1.zip, which is a N*3 array of coordinates for the spine ROIs where N is the number of spines (for each spine, there is also a dendritic ROI and background ROI drawn)
3. RoiSet1_dend.zip, which is contains an ROI of the dendritic skeleton and fiducial points
4. dendritic_distance folder (all files come from the output for https://github.com/surlab/dendritic-distance)

#### Running the code
While running the script, a display will show up asking if you want to check the spine distances between two timepoints: 
/var/folders/v8/z6zg9tq52c3ddwf7x3z3531r0000gn/T/TemporaryItems/NSIRD_screencaptureui_6wbpAu/Screen Shot 2024-11-20 at 2.57.29 AM.png
1. If you hit [Y]: 

-   A GUI will pop up with the two FOVs from the two timepoints and will run through each spine that was classified as missing (red spines in top image), retained (blue spines in both images), or added (green spines in bottom image)
https://www.dropbox.com/s/d0hvd0yly8y37y7/Screen%20Shot%202024-11-20%20at%203.06.56%20AM.png?dl=0
-   A M*2 (M = total number of retained, missing, added spines) editable table will also pop up with the classification for the spines identities for each timepoint (a spine that is missing or added will be noted as 0). You also can make any changes to the classifications, which will be saved under the updated "modify_mat" variable
https://www.dropbox.com/s/efy1zp4w1z8umfj/Screen%20Shot%202024-11-20%20at%203.08.05%20AM.png?dl=0

2. If you hit any other character, the GUI will exit

#### Output file
1. alignment_based_on_fiducial.mat which contains the original table with the classified ROI identities and "squeeze_mat" and the modified table "modify_mat" variable
2. Comparison_fig.fig is the saved figure of the GUI with the classified spines

#### Runtime
Depends on how many files, but typically takes < 1 minute to run



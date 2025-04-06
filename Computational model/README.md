The code was written in Matlab 2023b. It should work with other versions as well.

The figure 7 can be reproduced by simply running the script main.m on Matlab - results may slightly differ due to randomness.
Change the parameters ntrials (default = 20) to have different number of trials.
Each trial should take about 25-30 minutes to run.
Change sigma_het (default = 1.5) and sigma_heb (default = 15) to vary the heterosynaptic and the Hebbian factor, respectively. 
All the other parameters are as defined in the main paper. The code automatically saves the output for each trial in a .mat file, that can be then used to plot the results. 
To see an example without running new simulations please download the five .zip files "demoXX-YY.zip" containing results for 20 simulations (in total). The simulations were splitted to keep a size limit of 25mb. Please save all the .mat files from the .zip files in a single folder and run the other scripts as follows.

The script fig7B.m can be used to get an example as in Figure 7B.
The script fig7C_D_E.m can be used to reproduce panel C, D, and E in figure 7. 
Figure 7F can be reproduced together with Figure S8 through the script figS8.m in the following way: first, choose the parameters that determine the heterosynaptic and the Hebbian factors (sigma_het, sigma_heb, in the script main.m ). Run the main.m for each choice of (sigma_het, sigma_heb). Load the results in figS8.m and find average spine pair correlation and average mismatch in the somatic response for the given values of (sigma_het, sigma_heb) and multiple trials. 

The files preferred_orientation.m and test.m are functions to determine the soma’s preferred orientation and the somatic response for eight possible directions when plasticity is turned off, and they are used in the script figS8.m 


The .zip files contain .mat files with the simulations results used for fig S8.

The filenames are structured as follows:

data(row)(column)_(trial_range).zip --> if unzip, data(row)(column)_(trial_index).mat

Example:
data12_01-04.zip contains .mat files for trials 1, 2, 3, 4, for the row = 1 and the column = 2. data12_1.mat is the trial #1 for the row = 1, column = 2.

There are 20 trials for each row and column.

Rows and columns are indexed = 1, 2, 3 and they respectively represent the indexes for Sigma_Heb and Sigma_Het in Fig. S8.
Specifically,

Row = 1 corresponds to sigma_Heb = 0
Row = 2 corresponds to sigma_Heb = 15
Row = 3 corresponds to sigma_Heb = 30

Col = 1 corresponds to sigma_Heb = 0
Col = 2 corresponds to sigma_Heb = 1.5
Col = 3 corresponds to sigma_Heb = 5

Please unzip the folders, and use the script figS8.m in the folder 'Computational Model' to reproduce the results of FigS8. The results could differ slightly due to the 'test.m' function. 

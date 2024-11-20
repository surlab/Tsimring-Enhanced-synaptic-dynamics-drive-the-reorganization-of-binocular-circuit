function [ROIdata_fit,shaft_corr] = sub_chen(ROIdata_norm, framelength, numROIs,shaft)
%sub_chen, method from Chen et al https://www.nature.com/articles/nature12354

% (1) perform a robust fit (MATLAB) of spine DF/F against dendritic DF/F for stimulus-evoked data 
%       (a) indices of raw F that are below 2 sd of baseline become NaN
%       (getDFF_OFFperiod_NaN)
% (2) subtract a scaled version of the dendritic DF/F,where the scaling factor equals the slope from the robust fit.
warning('off','all'); %turn off warnings for duration of function
ROIdata_fit = zeros(framelength,numROIs);
for jj = 1:numROIs
    if shaft(jj,1) == 0
b = robustfit(ROIdata_norm(:,shaft(jj,2)), ROIdata_norm(:,jj));
ROIdata_fit(:,jj) = ROIdata_norm(:,jj) - b(2).*ROIdata_norm(:,shaft(jj,2));
    else
        ROIdata_fit(:,jj) = ROIdata_norm(:,jj);
     end
clear b        
end
clear jj

% Calculate remaining shaft correlation with only positive values
% (https://doi.org/10.1016/j.neuron.2017.10.017). 
ROIdata_nan = ROIdata_fit;
ROIdata_nan(ROIdata_fit<0) = NaN;
shaft_corr = nan(numROIs,1);
for jj = 1:numROIs 
    if shaft(jj,1) == 0
        shaft_corr(jj,1) = corr(ROIdata_nan(:,shaft(jj,2)), ROIdata_nan(:,jj),'Rows','complete');
    end
end
shaft_corr(isnan(shaft_corr)) = [];

%remove shaft ROIs
ROIdata_fit(:,shaft(:,1)==1) = [];
warning('on','all'); %turn warnings back on





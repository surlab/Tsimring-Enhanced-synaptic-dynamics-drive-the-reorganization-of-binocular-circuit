function [contra_ipsi_binoc] = concate_contra_ipsi_binoc(temp, features, join_key, for_fig2)
     if nargin < 4 || isempty(for_fig2)
        for_fig2 = 0; % default value
    end
    if ~isempty(find(strcmp([temp.session(:)], 'contra')))
            contra = temp(strcmp([temp.session(:)], 'contra'),features);
            ipsi = temp(strcmp([temp.session(:)], 'ipsi'),features);
            binoc = temp(strcmp([temp.session(:)], 'binoc'),features);
    else
        contra = temp(strcmp([temp.session{:}], "contra"),features);
        ipsi = temp(strcmp([temp.session{:}], "ipsi"),features);
        binoc = temp(strcmp([temp.session{:}], "binoc"),features);
    end
    contra_ipsi = outerjoin(contra, ipsi, "Keys",join_key, "MergeKeys",true);
    contra_ipsi_binoc = outerjoin(contra_ipsi, binoc, "Keys",join_key, "MergeKeys",true);
    
    % remove any missing values
    if for_fig2
        contra_ipsi_binoc=contra_ipsi_binoc(~any(ismissing(contra_ipsi_binoc),2),:);
    else
        contra_ipsi_binoc=contra_ipsi_binoc(~any(ismissing(contra_ipsi_binoc.resp_ipsi),2),:);
        contra_ipsi_binoc=contra_ipsi_binoc(~any(ismissing(contra_ipsi_binoc.resp_contra),2),:);
        contra_ipsi_binoc=contra_ipsi_binoc(~any(ismissing(contra_ipsi_binoc.resp),2),:);
    end

end
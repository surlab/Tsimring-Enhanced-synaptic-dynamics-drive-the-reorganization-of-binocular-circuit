function [contra_ipsi_binoc] = concate_contra_ipsi_binoc(temp, features, join_key)
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
    contra_ipsi_binoc=contra_ipsi_binoc(~any(ismissing(contra_ipsi_binoc),2),:);
end
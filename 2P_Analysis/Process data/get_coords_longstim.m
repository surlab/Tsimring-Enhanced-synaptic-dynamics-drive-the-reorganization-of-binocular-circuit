function get_coords_longstim(inputpath,coordname,FOV_name,savepath)

    Q=importdata(fullfile(inputpath,coordname));
    R = Q.data;
    R= R(:,2:end);
    
    x_coords=R(:,1:8:end);
    y_coords=R(:,2:8:end);
    if ~isdir(savepath)
        mkdir(savepath);
    end
   
    save(fullfile(savepath, [FOV_name,'.mat']), 'x_coords', 'y_coords');
end
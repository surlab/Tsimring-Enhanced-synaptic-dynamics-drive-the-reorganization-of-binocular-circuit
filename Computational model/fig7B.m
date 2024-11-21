n_ori = 8; lin_angles = linspace(0,360,n_ori+1);lin_angles(end) = []; rad_angles = deg2rad(lin_angles);
mytimes = [1 500 1000]; % before (1) , during (500), after (1000) plasticity
trial = 1;
load(['main' num2str(trial) '.mat']);

for indtime = 1:length(mytimes)
    timepoint= mytimes(indtime);
    % get values for the given time point
    myangles = [saving_angles{1}(timepoint,:); saving_angles{2}(timepoint,:)];
    myW = [saving_W{1}(timepoint,:); saving_W{2}(timepoint,:)];
    myview = [saving_eye{1}(timepoint,:); saving_eye{2}(timepoint,:)];
    myamp = [saving_ampVM{1}(timepoint,:); saving_ampVM{2}(timepoint,:)];
    myzVM = [saving_zVM{1}(timepoint,:); saving_zVM{2}(timepoint,:)];
    % test somatic response for 1) ipsi, 2) contra, 3) bino viewing
    [xtest_complete_L,xtest_complete_R,fr_I] = test(1,myangles,myW,myview,myamp,myzVM,dMat,dSoma,1,randVM);
    [xtest_complete_L,xtest_complete_R,fr_C] = test(2,myangles,myW,myview,myamp,myzVM,dMat,dSoma,1,randVM);
    [xtest_complete_L,xtest_complete_R,fr_all] = test(3,myangles,myW,myview,myamp,myzVM,dMat,dSoma,1,randVM);
    
    fr_abs = [fr_I; fr_C; fr_all];
    fr = fr_abs./(max(fr_abs,[],2)); 

    oris_rad = (lin_angles'.*pi)./180; % column vector
    num _oris = length(oris_rad)/2;
    oris_rad_half = oris_rad(1:num_oris);
    mean_data_half = (fr(:,1:num_oris) + fr(:,num_oris+1:end))/2;
    R_ori = mean_data_half*exp(oris_rad_half*2*1i);
    Ori_pref = mod(atan2(imag(R_ori),real(R_ori)),2*pi)*(90/pi);
    
    % Plot somatic response 
    figure() 
    for i = 1:3 % three types of view
        polarplot([0 2*deg2rad(Ori_pref(i))],[0 100],'LineWidth',3)
        hold on
    end
    thetaticks([0 90 180 270])
    thetaticklabels([0 45 90 135])
    rticks([0 25 50 75 100])
    rticklabels([])
    set(gca,'FontSize',24)
    %% Plot dendrite
    wth = 0.05;
    rng(0)
    N = size(myW,2);
    rand_pos = rand(N,2); 
    pos_test = [pos{1} pos{2}];
    
    cm = getWilsonMap(9);
    cm = [cm; 0.6 0.6 0.6]; colormap(cm);clim([0 10]);
    matrice = myangles; threshold = 1.5*wth;
    
    z = zeros(2,N);
    for ii=1:n_ori
       z(myangles==lin_angles(ii)) = ii;
    end
    z(myW<=threshold) = 10;
    
    matrice_W = 2*myW;
    figure('Renderer', 'painters', 'Position', [10 10 400 1200])
    p = bar(max(max(pos_test))+20,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7],'BarWidth',0.08);
    hold on
    p = bar(2,max(max(pos_test))+20,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7],'BarWidth',0.08);
    
    for i = 1:N
        ind = z(1,i);    col = cm(ind,:);
        plot(1+(-0.05)+0.1*rand_pos(i,1),pos_test(i,1),'o','MarkerFaceColor',col,'MarkerSize',10+matrice_W(1,i)*10,'MarkerEdgeColor','k');
    
        ind = z(2,i);
        col = cm(ind,:);
        plot(2+(-0.05)+0.1*rand_pos(i,2),pos_test(i,2),'o','MarkerFaceColor',col,'MarkerSize',10+matrice_W(2,i)*10,'MarkerEdgeColor','k');
    end
    xlim([0.2 2.8])
    set(gcf, 'Color', 'None')
    set(gca, 'Color', 'None')
    set(gca,'xtick',[])
    set(gca,'ytick',[],'FontSize',22)
end
    
function cmap = getWilsonMap(N) % colormap for plot
    rawMap = [213, 71, 146; 106, 84, 158; 64, 93, 166; 115, 199, 217; 116, 183, 72; 147, 193, 68; 248, 189, 22; 229, 42, 43; 213, 71, 146]/256;
    M = size(rawMap , 1);
    cmap = zeros(N , 3);
    cmap(: , 1) = interp1((0:M-1)/(M-1) , rawMap(:,1) , (0:N-1)/(N-1));
    cmap(: , 2) = interp1((0:M-1)/(M-1) , rawMap(:,2) , (0:N-1)/(N-1));
    cmap(: , 3) = interp1((0:M-1)/(M-1) , rawMap(:,3) , (0:N-1)/(N-1));
end
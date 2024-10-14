%%%%%%%
%
% Code for the paper 
% "Enhanced synaptic dynamics drive the reorganization 
% of binocular circuits in mouse visual cortex"
% Fig 7
%
%%%%%%%
clear all
close all
clc
%% Parameters
sigma_het = 1.5; % Spatial postsynaptic accumulator spread constant
sigma_heb = 15; % Spatial somatic accumulator spread constant
% Time parameters
dt = 1; % Time step (1 ms)
T = 1.8e+7; % Total time in ms
keep_angle = 250; 
T1 = 0:keep_angle:T; 
T1(1) = 1;
T2 = 0:dt:keep_angle;
% Time constants
tauw = 1000; taupre = 60; taupost = 30;
% Weight and distribution parameters
winit = 0.2; wmin = 0.05; wmax = 0.75;
% Lognormal distribution for weights
mu_w = 0.1; var_w = 0.05;
mu_log = log(mu_w^2 / sqrt(var_w + mu_w^2)); 
sigma_log = sqrt(log(var_w / mu_w^2 + 1));
% Lognormal distribution for spine distance
mu_dist = 5; var_dist = 5; 
mu_dist_log = log(mu_dist^2 / sqrt(var_dist + mu_dist^2)); 
sigma_dist_log = sqrt(log(var_dist / mu_dist^2 + 1));
min_dist = 0.05;
% Von Mises parameters for spine inputs
amp_min = 0.02; amp_max = 0.1; sig_amp = 0.01; 
z_min = 50; z_max = 100; sig_z = 0.75;
% Others
n_branch = 2; % Number of branches
N = 64; % Number of spines on each branch % option to create a vector with different spine numbers per branch
n_ori = 8; % Number of preferred orientations
lin_angles = 0:360/n_ori:315; % Angles in degrees
phi = 2; % scaling presynaptic inputs
rho = -0.13; % potentiation vs. depression threshold
A_thexp = 0.9; % decay constant for somatic accumulator threshold
pturnover = 0.5; % Probability that a new spine is formed

%% Initialization 
A = zeros(T,1); % somatic accumulator over time
A_branch = zeros(1,n_branch); % somatic accumulator on each branch
A_th = 0.5*ones(T,1); % initial bAP threshold 
bAP = zeros(T,1); % bAP vector
angles = ones(length(T1),1);
angles0 = cell(1,n_branch); 
eye0 = cell(1,n_branch);
angle_pref = cell(1,n_branch); 
eye_pref =  cell(1,n_branch); % over time 
spine_pref = cell(1,n_branch); % eye, ori
xt = cell(1,n_branch); r = cell(1,n_branch); % spike trains and rates
W = cell(1,n_branch);
PRE = cell(1,n_branch); POST = cell(1,n_branch); % pre- and postsyn accumulators
pos = cell(1,n_branch); dMat = cell(1,n_branch); dSoma = cell(1,n_branch); % spatial vars
% Von Mises inputs
ampVM = cell(1,n_branch); % over time
zVM = cell(1,n_branch); % over time
ampVM_now = cell(1,n_branch); 
zVM_now = cell(1,n_branch); 
randVM = cell(1,n_branch); % to determine direction selectivity
r_now = cell(1,n_branch);
angle_shown = 360*rand; 
% To keep track over time
added_times = cell(1,n_branch);
lost_times = cell(1,n_branch);
surv_times = cell(1,n_branch);
stored_ampVM = cell(1,n_branch);
stored_zVM = cell(1,n_branch);
stored_fr = cell(1,n_branch);
stored_prefangle = cell(1,n_branch);
stored_prefeye = cell(1,n_branch);
stored_PRE = cell(1,n_branch);
stored_POST = cell(1,n_branch);
for branch=1:n_branch
    % weights
    W{branch} = ones(T,N);
    W0 = winit + lognrnd(mu_log,sigma_log,1,N);
    W0 = W0.*(W0>=wmin).*(W0<=wmax) + (W0>wmax)*wmax + (W0<wmin)*1.5*wmin;
    W{branch}(1,:) = W0;
    % view preference
    eye_pref{branch} = zeros(T,N); 
    eye_temp = ones(1, N); 
    eye_temp(:, 1:floor(6.5 * N / 10)) = 2; % set 65% eye preference to 2 (representing 'contra')
    eye_temp(:, floor(6.5 * N / 10) + 1:floor(8.5 * N / 10) + 1) = 3; % Set 20% to 3 (representing 'binocular')
    eye_temp = eye_temp(randperm(numel(eye_temp))); % shuffle values
    eye_pref{branch}(1,:) = eye_temp; 
    % angle preference
    angle_pref{branch} = zeros(T,N);
    angles_temp = repmat(lin_angles, 1,[floor(N/n_ori)]);
    angles_temp = angles_temp(randperm(numel(angles_temp)));
    % tune strongest contra spines to same angle and ipsi to orthogonal
    strongthresh = quantile(W0, 0.2);
    angle_contra = lin_angles(randi(n_ori));
    angles_temp(W0 > strongthresh) = angle_contra;
    ipsispines = find((W0 > strongthresh) .* (eye_temp == 1) | (W0 > strongthresh) .* (eye_temp == 3));
    angles_temp(ipsispines) = mod(angle_contra - 90, 360);
    angle_pref{branch}(1,:) = angles_temp;
    spine_pref{branch} = [eye_temp; angles_temp];
    % inputs
    xt{branch} = zeros(T,N); % inputs
    r{branch} = ones(length(T1),N); % firing rates
    %accumulators
    PRE{branch} = zeros(T,N); 
    POST{branch} = zeros(T,N);
    % Define spine location on the branch
    pos{branch} = 30*ones(N,1); 
    dist = min_dist + lognrnd(mu_dist_log,sigma_dist_log,N-1,1);
    dist = 0.55*dist; % scaling
    for k=1:N-1
        pos{branch}(k+1) = pos{branch}(k) + dist(k);
    end
    dMat{branch} = normpdf(pdist2(pos{branch},pos{branch}), 0 , 2*sigma_het)*(sqrt(2*pi)*2*sigma_het);
    dSoma{branch} = normpdf(pdist2(0,pos{branch}), 0 , sigma_heb)*(sqrt(2*pi)*sigma_heb);
    ampVM{branch} = ones(length(T1),N);
    zVM{branch} = ones(length(T1),N);
    ampVM_now{branch} = amp_min + (amp_max-amp_min)*rand(1,N); 
    zVM_now{branch} = z_min + (z_max-z_min)*rand(1,N);
    ampVM{branch}(1,:) = ampVM_now{branch};
    zVM{branch}(1,:) = zVM_now{branch};
    randVM{branch} = rand(1,N);
    r_now{branch} = (ampVM{branch}(1,:).*exp(cos(((angles(1,:) - spine_pref{branch}(2,:)))/180*pi).*zVM{branch}(1,:))./(2*pi*besseli(0,zVM{branch}(1,:))) +...
                randVM{branch}.*ampVM{branch}(1,:).*exp(cos(((angles(1,:) - 180 - spine_pref{branch}(2,:)))/180*pi).*zVM{branch}(1,:))./(2*pi*besseli(0,zVM{branch}(1,:))) ) ./(1+randVM{branch});
    r{branch}(1,:) = r_now{branch};
    added_times{branch} = cell(1,N);
    added_times{branch} = cellfun(@(x) 1, added_times{branch}, 'UniformOutput', false);
    lost_times{branch} = cell(1,N);    surv_times{branch} = cell(1,N);
    stored_ampVM{branch} = cell(1,N);     stored_zVM{branch} = cell(1,N);     stored_fr{branch} = cell(1,N);
    stored_prefangle{branch} = cell(1,N);     stored_prefeye{branch} = cell(1,N);
    stored_PRE{branch} = cell(1,N);     stored_POST{branch} = cell(1,N);
end
% The following will track the input activity
temp_ampVM = zeros(n_branch,N);  
temp_zVM = zeros(n_branch,N); 
temp_fr = zeros(n_branch,N);

%% Start simulation
for tt = 2:T
    for branch=1:n_branch
        if mod(tt,keep_angle) == 0 % update Von Mises inputs every 250 ms
            time_ind = tt/keep_angle;
            % Update as OU process
            ampVM_now{branch} = ampVM{branch}(time_ind,:) + sig_amp*sqrt(1e-4)*randn(1,N);
            zVM_now{branch} = zVM{branch}(time_ind,:) + sig_z*sqrt(1e-2)*randn(1,N);
            angles(time_ind,:) = 360*rand; % angle shown (from continuous distribution on [0, 360] degrees
            r_now{branch} = (ampVM{branch}(time_ind,:).*exp(cos(((angles(time_ind,:) - spine_pref{branch}(2,:)))/180*pi).*zVM{branch}(time_ind,:))./(2*pi*besseli(0,zVM{branch}(time_ind,:))) +...
                randVM{branch}.*ampVM{branch}(time_ind,:).*exp(cos(((angles(time_ind,:) - 180 - spine_pref{branch}(2,:)))/180*pi).*zVM{branch}(time_ind,:))./(2*pi*besseli(0,zVM{branch}(time_ind,:))) ) ./(1+randVM{branch});
            ampVM{branch}(time_ind+1,:) = ampVM_now{branch};
            zVM{branch}(time_ind+1,:) = zVM_now{branch};
            r{branch}(time_ind+1,:) = r_now{branch};
        end
        W_now = W{branch}(tt-1,:) + (1./tauw)*(POST{branch}(tt-1,:).*(PRE{branch}(tt-1,:)+rho));
        W_now = W_now.*(W_now>=wmin);
        W_now = min(W_now, wmax);
        W{branch}(tt,:) = W_now;             
        Waug = repmat(W_now,[N,1]).*dMat{branch};
        xt{branch}(tt,rand(1,N) < r_now{branch} * dt)=1;
        Xold = xt{branch}(tt-1,:); 
        Xnew = xt{branch}(tt,:);
        PRE{branch}(tt,:) = PRE{branch}(tt-1,:) + (1./taupre)*(-PRE{branch}(tt-1,:) + phi*Xold);
        POST{branch}(tt,:) = POST{branch}(tt-1,:) +(1./taupost)*(-POST{branch}(tt-1,:) +(Waug*Xnew')'+ dSoma{branch}*bAP(tt-1));
        A_branch(branch) = W_now*POST{branch}(tt,:)';

        % Accumulate values for amp, z, and firing rates
        temp_ampVM(branch,:) = temp_ampVM(branch,:) + max(ampVM_now{branch},0);
        temp_zVM(branch,:) = temp_zVM(branch,:) + abs(zVM_now{branch});
        temp_fr(branch,:) = temp_fr(branch,:) + max(r_now{branch},0);
        % Identify lost spines (where W drops below threshold)
        idx_lost = (W{branch}(tt-1,:) >= wmin) & (W{branch}(tt,:) < wmin);
        % Process lost spines
        if any(idx_lost)
            lost_idx = find(idx_lost);  % Get indices of lost spines
            for idx = lost_idx
                lost_times{branch}{idx} = [lost_times{branch}{idx}, tt];
                lastadded = added_times{branch}{idx}(end);
                lastlost = lost_times{branch}{idx}(end);
                surv_time = lastlost - lastadded;
                surv_times{branch}{idx} = [surv_times{branch}{idx} surv_time];
                % Store the relevant data for lost spines
                stored_ampVM{branch}{idx} = [stored_ampVM{branch}{idx}, temp_ampVM(branch,idx) / surv_time];
                stored_zVM{branch}{idx} = [stored_zVM{branch}{idx}, temp_zVM(branch,idx) / surv_time];
                stored_fr{branch}{idx} = [stored_fr{branch}{idx}, temp_fr(branch,idx) / surv_time];
                stored_prefangle{branch}{idx} = [stored_prefangle{branch}{idx}, spine_pref{branch}(2,idx)];
                stored_prefeye{branch}{idx} = [stored_prefeye{branch}{idx}, spine_pref{branch}(1,idx)];
                stored_PRE{branch}{idx} = [stored_PRE{branch}{idx}, trapz(PRE{branch}(lastadded:tt,idx))];
                stored_POST{branch}{idx} = [stored_POST{branch}{idx}, trapz(POST{branch}(lastadded:tt,idx))];
            end
            % Reset temp accumulators for lost spines
            temp_ampVM(branch, idx_lost) = 0;
            temp_zVM(branch, idx_lost) = 0;
            temp_fr(branch, idx_lost) = 0;
        end      

        % Turnover logic: identify spines to be added
        turnover = (W_now<wmin) & (rand(1,N)<pturnover);
        idx_added = find(turnover > 0);  % Get indices of added spines
        % Process added spines
        if ~isempty(idx_added)
            added_times{branch}(idx_added) = cellfun(@(x) [x, tt], added_times{branch}(idx_added), 'UniformOutput', false);
            prob_newspine = repmat([0, 0.3, 0.5, 0.2],sum(sum(turnover)),1); % probabilities for eye pref (ipsi, contra, binoc)
            spine_pref{branch}(1,idx_added) = sum(rand(sum(sum(turnover)),1)>cumsum(prob_newspine,2),2);
            spine_pref{branch}(2,idx_added) = lin_angles(randi(n_ori,size(spine_pref{branch}(2,turnover>0))));
            ampVM{branch}(time_ind+1,idx_added) = amp_min + (amp_max-amp_min)*rand(1,sum(sum(turnover)));
            zVM{branch}(time_ind+1,idx_added) = z_min + (z_max-z_min)*rand(1,sum(sum(turnover)));
            W{branch}(tt,idx_added) = 1.5*wmin;
            eye_pref{branch}(tt,:) = spine_pref{branch}(1,:);
            angle_pref{branch}(tt,:) = spine_pref{branch}(2,:);
        end
    end
    % bAP 
    A(tt) = sum(A_branch);
    if A(tt)>A_th(tt)
        A_th(tt) = A_th(tt)+1;
        bAP(tt) = 1;
    end
    A_th(tt+1) = A_thexp*A_th(tt);
end
for branch=1:n_branch
    for idx = 1:N
        if (numel(added_times{branch}{idx}) - numel(lost_times{branch}{idx}))==1
            lost_times{branch}{idx}(end+1) = T;
            surv_time = lost_times{branch}{idx}(end)-added_times{branch}{idx}(end);
            surv_times{branch}{idx} = [surv_times{branch}{idx} surv_time];
            stored_ampVM{branch}{idx} = [stored_ampVM{branch}{idx} temp_ampVM(branch,idx)./surv_time];
            stored_zVM{branch}{idx} = [stored_zVM{branch}{idx} temp_zVM(branch,idx)./surv_time];
            stored_fr{branch}{idx} = [stored_fr{branch}{idx} temp_fr(branch,idx)./surv_time];
            stored_prefangle{branch}{idx} = [stored_prefangle{branch}{idx} spine_pref{branch}(2,idx)];
            stored_prefeye{branch}{idx} = [stored_prefeye{branch}{idx} spine_pref{branch}(1,idx)];
            stored_PRE{branch}{idx} = [stored_PRE{branch}{idx} trapz(PRE{branch}(added_times{branch}{idx}(end):lost_times{branch}{idx}(end),idx))];
            stored_POST{branch}{idx} = [stored_POST{branch}{idx} trapz(POST{branch}(added_times{branch}{idx}(end):lost_times{branch}{idx}(end),idx))];
        end
    end
end

%% Saving other variables 
% Number of time points to save
n_points = 1e+3;
time_indices1 = round(linspace(1, tt, n_points));
time_indices2 = round(linspace(1, time_ind, n_points));
% Initialize cells for each point
saving_W = cell(1, n_branch);
saving_ampVM = cell(1, n_branch);
saving_zVM = cell(1, n_branch);
saving_fr = cell(1, n_branch);
saving_eye = cell(1, n_branch);
saving_angles = cell(1, n_branch);
saving_PRE = cell(1, n_branch);
saving_POST = cell(1, n_branch);
saving_pos = pos;
saving_dMat = dMat;
saving_dSoma = dSoma;   
saving_A_th = A_th(time_indices1);
for branch = 1:n_branch
    saving_fr{branch} = r{branch}(time_indices2, :);
    saving_W{branch} = W{branch}(time_indices1, :);
    saving_PRE{branch} = PRE{branch}(time_indices1, :);
    saving_POST{branch} = POST{branch}(time_indices1, :);
    saving_zVM{branch} = zVM{branch}(time_indices2, :);
    saving_ampVM{branch} = ampVM{branch}(time_indices2, :);
    saving_eye{branch} = eye_pref{branch}(time_indices1, :);
    saving_angles{branch} = angle_pref{branch}(time_indices1, :);
end
function [x_L,x_R,fr_soma] = testing(eye,angles,W,view,amp,zvm,dMat,dSoma,bAP_switch,rndVM)

keep_angle = 2000; dt=1; taupre = 60; taupost = 30; % ms
N = size(W,2);
n_ori = 8; 
lin_angles = 0:360/n_ori:315; % Angles in degrees
phi = 2; % ms 
W_test = W;
amp_test = amp; amp_test(amp_test<0) = 0.001; % control
zVM_test = zvm;
angles_pref = angles;
dSoma_L = dSoma{1};  dSoma_R = dSoma{2};
Waug_L = repmat(W_test(1,:),[N,1]).*dMat{1}; 
Waug_R = repmat(W_test(2,:),[N,1]).*dMat{2}; 
view_idxL = [find(view(1,:) == eye) find(view(1,:) == 3)];
view_idxR = [find(view(2,:) == eye) find(view(2,:) == 3)];
x_L = zeros(n_ori*keep_angle,N);
x_R = zeros(n_ori*keep_angle,N);
if eye==3
    view_idxL = 1:N;
    view_idxR = 1:N;
end
rndtest_L = rndVM{1}(view_idxL);
rndtest_R = rndVM{2}(view_idxR);
fr_soma = zeros(1, n_ori);
r_L = zeros(1, N);
r_R = zeros(1, N);
% Tuning curve inputs
VM_curve = @(A, theta_shown, theta_pref, width, eps) (A.*exp(cos(((theta_shown - theta_pref))/180*pi).*width)./(2.*pi.*besseli(0,width)) +...
    eps.* A.*exp(cos(((theta_shown - 180 - theta_pref))/180*pi).*width)./(2.*pi.*besseli(0,width)) ) ./(1+eps);
for ii = 1:n_ori
    A_th = .8;
    angle_test = lin_angles(ii);
    xt_L = zeros(keep_angle,N); % inputs 
    xt_R = zeros(keep_angle,N); % inputs 
    r_L(view_idxL) = VM_curve(amp_test(1,view_idxL), angle_test, angles_pref(1,view_idxL), zVM_test(1,view_idxL), rndtest_L);
    r_R(view_idxR) = VM_curve(amp_test(2,view_idxR), angle_test, angles_pref(2,view_idxR), zVM_test(2,view_idxR), rndtest_R);
    A_test = zeros(keep_angle,1);
    bAP = zeros(keep_angle,1);
    PRE_L = zeros(keep_angle,N);
    PRE_R = zeros(keep_angle,N);
    POST_test = zeros(keep_angle,2,N);
    for tt = 2:keep_angle
        xt_L(tt,rand(1,N)<r_L*dt)=1;
        xt_R(tt,rand(1,N)<r_R*dt)=1;
        Xnew_L = xt_L(tt,view_idxL);
        Xnew_R = xt_R(tt,view_idxR);
        PRE_L(tt,view_idxL) = PRE_L(tt-1,view_idxL) + (1./taupre)*(-PRE_L(tt-1,view_idxL) + phi*Xnew_L);
        PRE_R(tt,view_idxR) = PRE_R(tt-1,view_idxR) + (1./taupre)*(-PRE_R(tt-1,view_idxR) + phi*Xnew_R);    
        POST_now = reshape(POST_test(tt-1,:,:),2,N);
        POST_now(1,view_idxL) = POST_now(1,view_idxL) + (1./taupost)*(-POST_now(1,view_idxL)+(Waug_L(view_idxL,view_idxL)*Xnew_L')'+dSoma_L(view_idxL)*bAP_switch*bAP(tt-1));
        POST_now(2,view_idxR) = POST_now(2,view_idxR)+(1./taupost)*(-POST_now(2,view_idxR)+(Waug_R(view_idxR,view_idxR)*Xnew_R')'+dSoma_R(view_idxR)*bAP_switch*bAP(tt-1));
        POST_test(tt,:,:) = POST_now;
        POST_now = reshape(POST_test(tt,:,:),2,N);
        A_test(tt) = W_test(1,view_idxL)*(POST_now(1,view_idxL)') + W_test(2,view_idxR)*(POST_now(2,view_idxR)');
        if A_test(tt)>A_th
            A_th = A_th+.5;
            bAP(tt) = 1;
        end
        A_th = A_th*.9;
    end
    fr_soma(ii) =  sum(bAP);
    x_L((keep_angle*(ii-1)+1):keep_angle*ii,:) = PRE_L;
    x_R((keep_angle*(ii-1)+1):keep_angle*ii,:) = PRE_R;
end
end
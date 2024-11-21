function [xtest_complete_L,xtest_complete_R,fr_soma] = test(eye,angles,W,view,amp,kvm,dMat,dSoma,bap_switch,rndVM)

keep_angle_test = 2000; dt=1; taupre = 60; taupost = 30; % ms
n_ori = 8;    lin_angles = linspace(0,360,n_ori+1); lin_angles(end) = []; 
    
phi = 2; % ms 
N = size(W,2);
% BEFORE
angles_test = angles;
dSoma_left = dSoma{1}; 
dSoma_right = dSoma{2};
W_test = W;
view_test = view;

amp_test = amp;
amp_test(amp_test<0) = 0.001;
k_VM_test = kvm;

Waug_left = repmat(W_test(1,:),[N,1]).*dMat{1}; 
Waug_right = repmat(W_test(2,:),[N,1]).*dMat{2}; 

view_indexes_left = [find(view_test(1,:) == eye) find(view_test(1,:) == 3)];
view_indexes_right = [find(view_test(2,:) == eye) find(view_test(2,:) == 3)];% find(view_test(W_test>1.5*wth) == eye);
xtest_complete_L = zeros(n_ori*keep_angle_test,N);
xtest_complete_R = zeros(n_ori*keep_angle_test,N);

if eye==3
    view_indexes_left = 1:N;
    view_indexes_right = 1:N;
end
rndtest_left = rndVM{1}(view_indexes_left);
rndtest_right = rndVM{2}(view_indexes_right);

fr_soma = zeros(1, n_ori);
r_test_left = zeros(1, N);
r_test_right = zeros(1, N);

for ii = 1:n_ori
    A_th = .8;
    alpha_angle = lin_angles(ii);
    xt_testleft = zeros(keep_angle_test,N); % inputs 
    xt_testright = zeros(keep_angle_test,N); % inputs 
    r_test_left(view_indexes_left) = (amp_test(1,view_indexes_left).*exp(cos(((alpha_angle- angles_test(1,view_indexes_left)))/180*pi).*k_VM_test(1,view_indexes_left))./(2.*pi.*besseli(0,k_VM_test(1,view_indexes_left))) +...
    rndtest_left.* amp_test(1,view_indexes_left).*exp(cos(((alpha_angle- 180 -angles_test(1,view_indexes_left)))/180*pi).*k_VM_test(1,view_indexes_left))./(2.*pi.*besseli(0,k_VM_test(1,view_indexes_left))) ) ./(1+rndtest_left);
    r_test_right(view_indexes_right) = (amp_test(2,view_indexes_right).*exp(cos(((alpha_angle- angles_test(2,view_indexes_right)))/180*pi).*k_VM_test(2,view_indexes_right))./(2.*pi.*besseli(0,k_VM_test(2,view_indexes_right))) +...
           rndtest_right.*amp_test(2,view_indexes_right).*exp(cos(((alpha_angle- 180-angles_test(2,view_indexes_right)))/180*pi).*k_VM_test(2,view_indexes_right))./(2.*pi.*besseli(0,k_VM_test(2,view_indexes_right))) ) ./(1+rndtest_right);
    A_test = zeros(keep_angle_test,1);
    bAP = zeros(keep_angle_test,1);
    PRE_test_left = zeros(keep_angle_test,N);
    PRE_test_right = zeros(keep_angle_test,N);
    POST_test = zeros(keep_angle_test,2,N);
    for tt = 2:keep_angle_test
        xt_testleft(tt,rand(1,N)<r_test_left*dt)=1;
        xt_testright(tt,rand(1,N)<r_test_right*dt)=1;
        Xnew_left = xt_testleft(tt,view_indexes_left);
        Xnew_right = xt_testright(tt,view_indexes_right);
        PRE_test_left(tt,view_indexes_left) = PRE_test_left(tt-1,view_indexes_left) + (1./taupre)*(-PRE_test_left(tt-1,view_indexes_left) + phi*Xnew_left);
        PRE_test_right(tt,view_indexes_right) = PRE_test_right(tt-1,view_indexes_right) + (1./taupre)*(-PRE_test_right(tt-1,view_indexes_right) + phi*Xnew_right);
      
        POST_now = reshape(POST_test(tt-1,:,:),2,N);
        POST_now(1,view_indexes_left) = POST_now(1,view_indexes_left) + (1./taupost)*(-POST_now(1,view_indexes_left)+(Waug_left(view_indexes_left,view_indexes_left)*Xnew_left')'+dSoma_left(view_indexes_left)*bap_switch*bAP(tt-1));
        POST_now(2,view_indexes_right) = POST_now(2,view_indexes_right)+(1./taupost)*(-POST_now(2,view_indexes_right)+(Waug_right(view_indexes_right,view_indexes_right)*Xnew_right')'+dSoma_right(view_indexes_right)*bap_switch*bAP(tt-1));
        POST_test(tt,:,:) = POST_now;
        POST_now = reshape(POST_test(tt,:,:),2,N);
        A_test(tt) = W_test(1,view_indexes_left)*(POST_now(1,view_indexes_left)') + W_test(2,view_indexes_right)*(POST_now(2,view_indexes_right)');
        if A_test(tt)>A_th
            A_th = A_th+.5;
            bAP(tt) = 1;
        end
        A_th = A_th*.9;
    end
    fr_soma(ii) =  sum(bAP);
    xtest_complete_L((2000*(ii-1)+1):2000*ii,:) = PRE_test_left;
    xtest_complete_R((2000*(ii-1)+1):2000*ii,:) = PRE_test_right;
end
end

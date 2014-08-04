clear all
% close all
addpath(genpath('~/James_scripts'));
cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-7 7;-7 7];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;

Pix2Deg = 0.018837;
% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]

RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

n_iter = 3;
deltax_range = 12;
deltay_range = 12;
smooth_sigma = 0.015;
error_sigma = 0.3;

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
% load fixation_stim_righteye
load fixation_data_v3
% load fixation_imagXe_patches
load fixation_image_patches_d2_nt_origcal
disparity = all_sac_locs_right_nt - all_sac_locs_left_nt;

%%
sm_win = 10;
sm_hor = smooth(disparity(:,1),sm_win,'rlowess');
sm_ver = smooth(disparity(:,2),sm_win,'rlowess');

disp_error(:,1) = disparity(:,1)-sm_hor;
disp_error(:,2) = disparity(:,2)-sm_ver;

max_hor_error = 0.2;
max_ver_error = 0.15;
bad_disp = find(abs(disp_error(:,1)) > max_hor_error | abs(disp_error(:,2)) > max_ver_error);
good_disp = setdiff(1:size(disp_error,1),bad_disp);

% to_use = find(ismember(used_fixs,good_disp));
% X_left = X_left(to_use,:,:);
% X_right = X_right(to_use,:,:);
% X_avg = X_avg(to_use,:,:);
% used_fixs = used_fixs(to_use);

used_left_pos = all_sac_locs_left_nt(used_fixs,:);
used_right_pos = all_sac_locs_right_nt(used_fixs,:);
used_disp = used_right_pos - used_left_pos;
avg_disp = mean(used_disp);
used_left_pos = bsxfun(@plus,used_left_pos,0.5*avg_disp);
used_right_pos = bsxfun(@minus,used_right_pos,0.5*avg_disp);
used_avg_pos = 0.5*used_left_pos + 0.5*used_right_pos;
used_left_err = used_left_pos - used_avg_pos;
used_right_err = used_right_pos - used_avg_pos;
used_disp = used_right_pos - used_left_pos;

%%
X = X_avg;
NT = size(X,1);
SDIM = size(X,2);
kern_len = SDIM-1;
[Xg,Yg] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);

cellids = [1 2 3 4 5 6 7 8 10];
hor_posd = [4.5 4.72 4.5 4.65 4.65 4.65 4.88 4.57 4.57];
vert_posd = [-2.61 -2 -1.93 -2.08 -1.7 -1.85 -1.7 -1.85 -1.63];
orientations = [0 0.31 0.31 0 0 1.9 pi/2 1.9 0];
lambdas = [1.0 0.35 0.35 0.5 0.67 0.35 0.67 0.35 0.67]*Fsd;
phases = [0 pi/2];
% vert_loc = [12 20 21 19 24 22 24 22 25];
% hor_loc = [14 17 14 16 16 16 19 15 15];
for i = 1:length(cellids)
    [~,vert_ind(i)] = min(abs(yax(ypatch_inds) - vert_posd(i)));
    [~,hor_ind(i)] = min(abs(xax(xpatch_inds) - hor_posd(i)));
end
vert_pos = Yg(vert_ind,1);
hor_pos = Xg(1,hor_ind);

muaids = [1 2 4 5 6 11 12 13 14];
mhor_posd = [4.5 4.57 nan 4.65 4.65 4.6 nan nan 4.9 4.72 4.72 4.7 4.7 4.88];
mvert_posd = [-2.0 -1.78 nan -2.0 -2.0 -1.85 nan nan -2 -1.78 -1.85 -1.9 -1.6 -1.4];
morientations = [0.39 0.39 nan 0.39 0 0 nan nan 0 1.57 0.79 0 0 0];
mlambdas = [0.5 0.5 nan 0.36 0.36 0.5 nan nan 0.5 0.36 0.67 0.5 0.67 1]*Fsd;
mhor_posd = mhor_posd(muaids); mvert_posd = mvert_posd(muaids);
morientations = morientations(muaids); mlambdas = mlambdas(muaids);
phases = [0 pi/2];
for i = 1:length(muaids)
    [~,mvert_ind(i)] = min(abs(yax(ypatch_inds) - mvert_posd(i)));
    [~,mhor_ind(i)] = min(abs(xax(xpatch_inds) - mhor_posd(i)));
end
mvert_pos = Yg(mvert_ind,1);
mhor_pos = Xg(1,mhor_ind);

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];
hor_pos = [hor_pos mhor_pos];
vert_pos = [vert_pos mvert_pos];
orientations = [orientations morientations];
lambdas = [lambdas mlambdas];

n_used_cells = size(spk_cnts,2);

%% COMPUTE GABOR TEMPLATES
NT = size(X,1);
SDIM = size(X,2);
kern_len = SDIM-1;
[Xg,Yg] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);

Xmat = reshape(X,NT,SDIM^2);

for n = 1:n_used_cells;
    fprintf('Fitting cell %d of %d\n',n,n_used_cells);
    init_params(1) = hor_pos(n); %x0
    init_params(2) = vert_pos(n); %y0
    init_params(3) = orientations(n); %theta
    init_params(4) = lambdas(n); %lambda
    init_params(5) = 0.5*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %weight of lin term 0 phase
    init_params(9) = 0; %weight of lin term pi/2 phase
    init_params(10) = 0; %const offset
    ip(n,:) = init_params(1:6);
    hold_const = [0 0 0 0 0 0 0 0 0 0];
    LB = [-5 -5 0 6 2 0.5 -Inf -Inf -Inf -Inf];
    UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf Inf];
    [gabor_params{1}(n,:),LL{1}(n)] = fit_gabor_params(init_params,Xmat,spk_cnts(used_fixs,n),[SDIM SDIM],hold_const,LB,UB);
    
    pref_phase(n) = atan2(gabor_params{1}(n,9),gabor_params{1}(n,8));
    cur_mask1 = get_pgabor_mask(gabor_params{1}(n,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params{1}(n,1:6),pi/2,[SDIM SDIM]);
    mask1_out = Xmat*cur_mask1(:);
    mask2_out = Xmat*cur_mask2(:);
    e_out = sqrt(mask1_out.^2 + mask2_out.^2);
    gmat = gabor_params{1}(n,7)*e_out + gabor_params{1}(n,8)*mask1_out + gabor_params{1}(n,9)*mask2_out;
    [temp_scale{1}(n,:),temp_offset{1}(n,:),temp_LL{1}(n)] = fit_gnm_gabormodel_blocks_v3(gmat,spk_cnts(used_fixs,n),blockids(used_fixs));
end
%%
% for i = 1:n_used_cells
%     subplot(5,4,i)
%     cur_mask_fin = get_pgabor_mask(gabor_params{1}(i,1:6),0,[SDIM SDIM]);
%     imagesc(cur_mask_fin);
% end

%% COMPUTE GABOR FILTER BANKS
delta_x = -deltax_range:deltax_range;
delta_y = -deltay_range:deltay_range;
use_x = find(ismember(Xg(1,:),delta_x)); use_x = use_x(:);
use_y = find(ismember(Yg(:,1),delta_y)); use_y = use_y(:);
[dX,dY] = meshgrid(delta_x,delta_y);
dX = dX(:);
dY = dY(:);
K = length(dX);

gabor_bank = zeros(n_used_cells,2,SDIM,SDIM);
for n = 1:n_used_cells
    gabor_bank(n,1,:,:) = get_pgabor_mask(gabor_params{1}(n,1:6),0,[SDIM SDIM]);
    gabor_bank(n,2,:,:) = get_pgabor_mask(gabor_params{1}(n,1:6),pi/2,[SDIM SDIM]);
end

%% CONVOLVE GABOR TEMPLATES WITH IMAGE PATCH
X_gabor_output = zeros([n_used_cells 2 NT length(delta_y) length(delta_x)]);
resh_X = reshape(X,NT,SDIM^2);
for n = 1:n_used_cells
    fprintf('Computing Gabor output matricies, Cell %d of %d\n',n,n_used_cells);
    cur_mask1 = squeeze(gabor_bank(n,1,:,:));
    cur_mask2 = squeeze(gabor_bank(n,2,:,:));
    for i = 1:length(delta_x)
        for j = 1:length(delta_y)
            smask1 = dist_shift2d(cur_mask1,delta_x(i),2,0);
            smask1 = dist_shift2d(smask1,delta_y(j),1,0);
            smask2 = dist_shift2d(cur_mask2,delta_x(i),2,0);
            smask2 = dist_shift2d(smask2,delta_y(j),1,0);
            smask1 = smask1(:);
            smask2 = smask2(:);
            X_gabor_output(n,1,:,j,i) = resh_X*smask1;
            X_gabor_output(n,2,:,j,i) = resh_X*smask2;
        end
    end
end

%%
X_gabor_output = reshape(X_gabor_output,[n_used_cells 2 NT K]);

zero_state = find(dX==0 & dY==0);
cur_state_seq = ones(NT,1)*zero_state;
cur_seq_dX = dX(cur_state_seq);
cur_seq_dY = dY(cur_state_seq);

dist_mat = squareform(pdist([dX dY]));
dist_mat = dist_mat/Fsd;
A = exp(-dist_mat.^2/(2*smooth_sigma.^2));
for i = 1:K
    A(i,:) = A(i,:)/sum(A(i,:));
end
state_dists_left = bsxfun(@plus,used_left_err(:,1),dX'/Fsd).^2 + bsxfun(@plus,used_left_err(:,2),dY'/Fsd).^2;
state_priors_left = exp(-state_dists_left.^2/(2*error_sigma^2));
state_dists_right = bsxfun(@plus,used_right_err(:,1),dX'/Fsd).^2 + bsxfun(@plus,used_right_err(:,2),dY'/Fsd).^2;
state_priors_right = exp(-state_dists_right.^2/(2*error_sigma^2));

%%
LLs = zeros(n_used_cells,NT,K);
for n = 1:n_used_cells
    LLs(n,:,:)= get_pgabor_tempmod_LLs_v2(gabor_params{1}(n,:),temp_scale{1}(n,:),temp_offset{1}(n,:),spk_cnts(used_fixs,n),...
        squeeze(X_gabor_output(n,:,:,:,:)),blockids(used_fixs));
end

clear cur_X_seq cur_Y_seq loglik cur_model cur_r2 post_entropy post_stdX post_stdY
for n = 2:n_iter
    fprintf('Iteration %d of %d\n',n,n_iter);
    
    B = squeeze(exp(sum(LLs,1))); %convert to overall likes
    
    Pi = ones(K,1);
    Pi = Pi'/sum(Pi);
    
    %initialize forward and backward messages
    alpha=zeros(NT,K);
    beta=zeros(NT,K);
    gamma=zeros(NT,K);
    
    scale=zeros(NT,1); %initialize rescaling parameters
    
    disp('Computing forward and backward messages');
    %compute rescaled forward messages
    alpha(1,:)=Pi.*B(1,:);
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)/scale(1);
    for t=2:NT
        alpha(t,:)=(alpha(t-1,:)*A).*B(t,:);
        scale(t) = sum(alpha(t,:));
        alpha(t,:)=alpha(t,:)/scale(t);
    end
    
    %compute rescaled backward messages
    beta(NT,:)=ones(1,K)/scale(NT);
    for t=NT-1:-1:1
        beta(t,:)=(beta(t+1,:).*B(t+1,:))*(A')/scale(t);
    end
    
    %compute posteriors over hidden states
    gamma=(alpha.*beta);
    gamma = gamma.*state_priors_left.*state_priors_right; %multiply by state prior matrices
    gamma = gamma./repmat(sum(gamma,2),1,K);
    
    lscale = log(scale);
    loglik(n) = sum(lscale); %rabiner eq(103), scale is defined as inverse here though
    post_entropy(n,:) = sum(-gamma.*log2(gamma),2);
    post_meanX = sum(bsxfun(@times,gamma,dX'),2);
    post_meanY = sum(bsxfun(@times,gamma,dY'),2);
    post_diffX = bsxfun(@minus,post_meanX,dX');
    post_diffY = bsxfun(@minus,post_meanY,dY');
    post_stdX(n,:) = sqrt(sum(gamma.*post_diffX.^2,2));
    post_stdY(n,:) = sqrt(sum(gamma.*post_diffY.^2,2));
    
    %NOW FOR VITERBI SEQUENCE ESTIMATION
    fprintf('Computing viterbi sequence\n');
    %initialize variables
    delta=zeros(NT,K); %delta is the maximized log probability as a function of time
    psi=zeros(NT,K); %Psi stores the most likely preceeding state
    %initialization
    delta(1,:) = log(Pi)+log(B(1,:));    % Eq. 105(a) Rabiner
    for t=2:NT
        temp = repmat(delta(t-1,:),K,1) + log(A);
        [delta(t,:),psi(t,:)] = max(temp,[],2);
        delta(t,:) = delta(t,:) + log(B(t,:)) + log(state_priors_left(t,:)) + log(state_priors_right(t,:));
    end
    
    % Backtracking Viterbi
    state_seq = zeros(NT,1);
    [llik_best,state_seq(NT)] = max(delta(NT,:));
    for t=NT-1:-1:1,
        state_seq(t) = psi(t+1,state_seq(t+1));
    end
    
    cur_X_seq(n,:) = dX(state_seq);
    cur_Y_seq(n,:) = dY(state_seq);
    
    %% UPDATE MODEL
    disp('Updating model');
    SDIM = size(X,2);
    Xmat_resh = reshape(X,size(X,1),SDIM,SDIM);
    for i = 1:size(X,1)
        Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-cur_X_seq(end,i),2,0);
        Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-cur_Y_seq(end,i),1,0);
    end
    Xmatr = reshape(Xmat_resh,NT,SDIM^2);
    for c = 1:n_used_cells;
        fprintf('Fitting cell %d of %d\n',c,n_used_cells);
        hold_const = [0 0 0 0 0 0 0 0 0 0];
        [gabor_params{n}(c,:),LL{n}(c)] = fit_gabor_params(gabor_params{n-1}(c,:),Xmatr,spk_cnts(used_fixs,c),[SDIM SDIM],hold_const,LB,UB);
        
        pref_phase(c) = atan2(gabor_params{n}(c,9),gabor_params{n}(c,8));
        cur_mask1 = get_pgabor_mask(gabor_params{n}(c,1:6),0,[SDIM SDIM]);
        cur_mask2 = get_pgabor_mask(gabor_params{n}(c,1:6),pi/2,[SDIM SDIM]);
        mask1_out = Xmat*cur_mask1(:);
        mask2_out = Xmat*cur_mask2(:);
        e_out = sqrt(mask1_out.^2 + mask2_out.^2);
        gmat = gabor_params{n}(c,7)*e_out + gabor_params{n}(c,8)*mask1_out + gabor_params{n}(c,9)*mask2_out;
        [temp_scale{n}(c,:),temp_offset{n}(c,:),temp_LL{n}(c)] = fit_gnm_gabormodel_blocks_v3(gmat,spk_cnts(used_fixs,c),blockids(used_fixs));
    end
    
    if n < n_iter
        gabor_bank = zeros(n_used_cells,2,SDIM,SDIM);
        for c = 1:n_used_cells
            gabor_bank(c,1,:,:) = get_pgabor_mask(gabor_params{n}(c,1:6),0,[SDIM SDIM]);
            gabor_bank(c,2,:,:) = get_pgabor_mask(gabor_params{n}(c,1:6),pi/2,[SDIM SDIM]);
        end
        
        X_gabor_output = zeros([n_used_cells 2 NT length(delta_y) length(delta_x)]);
        for c = 1:n_used_cells
            fprintf('Computing Filter Outputs Cell %d of %d\n',c,n_used_cells);
            cur_mask1 = squeeze(gabor_bank(c,1,:,:));
            cur_mask2 = squeeze(gabor_bank(c,2,:,:));
            for i = 1:length(delta_x)
                for j = 1:length(delta_y)
                    smask1 = dist_shift2d(cur_mask1,delta_x(i),2,0);
                    smask1 = dist_shift2d(smask1,delta_y(j),1,0);
                    smask2 = dist_shift2d(cur_mask2,delta_x(i),2,0);
                    smask2 = dist_shift2d(smask2,delta_y(j),1,0);
                    smask1 = smask1(:);
                    smask2 = smask2(:);
                    X_gabor_output(c,1,:,j,i) = resh_X*smask1;
                    X_gabor_output(c,2,:,j,i) = resh_X*smask2;
                end
            end
        end
        X_gabor_output = reshape(X_gabor_output,[n_used_cells 2 NT K]);
        
        LLs = zeros(n_used_cells,NT,K);
        for c = 1:n_used_cells
            LLs(c,:,:)= get_pgabor_tempmod_LLs_v2(gabor_params{n}(c,:),temp_scale{n}(c,:),temp_offset{n}(c,:),spk_cnts(used_fixs,c),...
                squeeze(X_gabor_output(c,:,:,:,:)),blockids(used_fixs));
        end
    end
end

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
save sing_eye_pgabortrack_d2 used_fixs loglik cur_*_seq dX dY post_std* post_entropy gabor_params LL temp_scale temp_offset temp_LL

%%
use_it = n_iter;
X_mean = cur_X_seq(use_it,:)/Fsd;
X_std = post_stdX(use_it,:)/Fsd;
Y_mean = cur_Y_seq(use_it,:)/Fsd;
Y_std = post_stdY(use_it,:)/Fsd;

est_X_pos = used_avg_pos(:,1)+X_mean';
est_Y_pos = used_avg_pos(:,2)+Y_mean';
orig_errors = used_disp/2;
est_error_left = [(used_left_pos(:,1)-est_X_pos) (used_left_pos(:,2)-est_Y_pos)];
est_error_right = [(used_right_pos(:,1)-est_X_pos) (used_right_pos(:,2)-est_Y_pos)];

figure
hold on
plot(used_fixs,orig_errors(:,1),'b')
plot(used_fixs,-orig_errors(:,1),'k')
h = shadedErrorBar(used_fixs,X_mean,X_std,'r');

figure
hold on
plot(used_fixs,orig_errors(:,2),'b')
plot(used_fixs,-orig_errors(:,2),'k')
h = shadedErrorBar(used_fixs,Y_mean,Y_std,'r');

figure
hold on
h = shadedErrorBar(used_fixs,-est_error_left(:,2),Y_std,'r');
h = shadedErrorBar(used_fixs,-est_error_right(:,2),Y_std,'b');

%%
% load sing_eye_pgabortrack_d2
% lam_mat = zeros(n_used_cells,n_iter);
% sig_mat = zeros(n_used_cells,n_iter);
% ecc_mat = zeros(n_used_cells,n_iter);
% wc_mat = zeros(n_used_cells,n_iter);
% ws_mat = zeros(n_used_cells,n_iter);
% for i = 1:n_iter
%     lam_mat(:,i) = gabor_params{i}(:,4);
%     sig_mat(:,i) = gabor_params{i}(:,5);
%     ecc_mat(:,i) = gabor_params{i}(:,6);
%     wc_mat(:,i) = gabor_params{i}(:,7);
%     ws_mat(:,i) = max(abs(gabor_params{i}(:,8:9)),[],2);
% end
% figure
% subplot(2,2,1)
% plot(mean(lam_mat))
% 
% subplot(2,2,2)
% plot(mean(sig_mat))
% 
% subplot(2,2,3)
% plot(mean(wc_mat))
% subplot(2,2,4)
% plot(mean(ws_mat))

%%

gabor_bank = zeros(n_used_cells,2,SDIM,SDIM);
for c = 1:n_used_cells
    gabor_bank(c,1,:,:) = get_pgabor_mask(gabor_params{end}(c,1:6),0,[SDIM SDIM]);
    gabor_bank(c,2,:,:) = get_pgabor_mask(gabor_params{end}(c,1:6),pi/2,[SDIM SDIM]);
end

X_gabor_output = zeros([n_used_cells 2 NT length(delta_y) length(delta_x)]);
for c = 1:n_used_cells
    fprintf('Computing Filter Outputs Cell %d of %d\n',c,n_used_cells);
    cur_mask1 = squeeze(gabor_bank(c,1,:,:));
    cur_mask2 = squeeze(gabor_bank(c,2,:,:));
    for i = 1:length(delta_x)
        for j = 1:length(delta_y)
            smask1 = dist_shift2d(cur_mask1,delta_x(i),2,0);
            smask1 = dist_shift2d(smask1,delta_y(j),1,0);
            smask2 = dist_shift2d(cur_mask2,delta_x(i),2,0);
            smask2 = dist_shift2d(smask2,delta_y(j),1,0);
            smask1 = smask1(:);
            smask2 = smask2(:);
            X_gabor_output(c,1,:,j,i) = resh_X*smask1;
            X_gabor_output(c,2,:,j,i) = resh_X*smask2;
        end
    end
end
X_gabor_output = reshape(X_gabor_output,[n_used_cells 2 NT K]);

LLs = zeros(n_used_cells,NT,K);
for c = 1:n_used_cells
    LLs(c,:,:)= get_pgabor_tempmod_LLs_v2(gabor_params{end}(c,:),temp_scale{end}(c,:),temp_offset{end}(c,:),spk_cnts(used_fixs,c),...
        squeeze(X_gabor_output(c,:,:,:,:)),blockids(used_fixs));
end

cur_state_err_X = bsxfun(@minus,X_mean',dX'/Fsd);
cur_state_err_Y = bsxfun(@minus,Y_mean',dY'/Fsd);
error_sigma = 0.075;
state_priors = exp(-(cur_state_err_X.^2+cur_state_err_Y.^2)/(2*error_sigma^2));

B = squeeze(exp(sum(LLs,1))); %convert to overall likes
posteriors = B.*state_priors;
posteriors = posteriors./repmat(sum(posteriors,2),1,K);

[~,post_max] = max(posteriors,[],2);
post_maxX = dX(post_max)/Fsd;
post_maxY = dY(post_max)/Fsd;
post_meanX = sum(bsxfun(@times,posteriors,dX'),2);
post_meanY = sum(bsxfun(@times,posteriors,dY'),2);
post_diffX = bsxfun(@minus,post_meanX,dX');
post_diffY = bsxfun(@minus,post_meanY,dY');
post_stdX = sqrt(sum(posteriors.*post_diffX.^2,2))/Fsd;
post_stdY = sqrt(sum(posteriors.*post_diffY.^2,2))/Fsd;

% figure
% shadedErrorBar(used_fixs,post_maxX,post_stdX);
% 
% figure
% shadedErrorBar(used_fixs,post_maxY,post_stdY,'r');
tot_std_X = 1/sqrt(2)*(X_std' + post_stdX);
tot_std_Y = 1/sqrt(2)*(Y_std' + post_stdY);

tot_est_X = used_avg_pos(:,1) + post_maxX;
tot_est_Y = used_avg_pos(:,2) + post_maxY;

orig_errors = used_disp/2;
est_error_left = [(used_left_pos(:,1)-tot_est_X) (used_left_pos(:,2)-tot_est_Y)];
est_error_right = [(used_right_pos(:,1)-tot_est_X) (used_right_pos(:,2)-tot_est_Y)];

figure
subplot(2,1,1)
shadedErrorBar(used_fixs,post_maxX,tot_std_X,'b');
hold on
plot(used_fixs,X_mean,'r','linewidth',2)
axis tight
xlabel('Fixation number','fontsize',14)
ylabel('Estimated error (degrees)','fontsize',14)
subplot(2,1,2)
plot(used_fixs,used_disp(:,1),'k')
hold on
plot(used_fixs,est_error_left(:,1),'b')
plot(used_fixs,est_error_right(:,1),'r')
xlabel('Fixation number','fontsize',14)
ylabel('Measured disparity (degrees)','fontsize',14)
axis tight
legend('Observed Disparity','Left error','Right error')

figure
subplot(2,1,1)
shadedErrorBar(used_fixs,post_maxY,tot_std_Y,'b');
hold on
plot(used_fixs,Y_mean,'r','linewidth',2)
axis tight
xlabel('Fixation number','fontsize',14)
ylabel('Estimated error (degrees)','fontsize',14)
subplot(2,1,2)
plot(used_fixs,used_disp(:,2),'k')
hold on
plot(used_fixs,est_error_left(:,2),'b')
plot(used_fixs,est_error_right(:,2),'r')
xlabel('Fixation number','fontsize',14)
ylabel('Measured disparity (degrees)','fontsize',14)
axis tight
legend('Observed Disparity','Left error','Right error')

%%
close all

figure
set(gca,'fontsize',14,'fontname','arial')
subplot(2,1,1)
shadedErrorBar(used_fixs,post_maxY,tot_std_Y,'k');
hold on
h = shadedErrorBar(used_fixs,Y_mean,Y_std,'r');
ylim([-0.5 0.5])
xlim([0 used_fixs(end)])
xlabel('Fixation Number','fontsize',14)
ylabel('Estimated Error (degrees)','fontsize',14)
title('Total Vertical Error Estimate','fontsize',16)
subplot(2,1,2)
h = shadedErrorBar(used_fixs,post_maxY-Y_mean',post_stdY,'k');
xlim([0 used_fixs(end)])
ylim([-0.5 0.5])
xlabel('Fixation Number','fontsize',14)
ylabel('Estimated Error (degrees)','fontsize',14)
title('High-frequency Vertical Error Estimate','fontsize',16)

hax = linspace(-0.5,0.5,25);
temp1 = hist(post_maxY,hax);
temp2 = hist(post_maxY-Y_mean',hax);
figure
set(gca,'fontsize',14,'fontname','arial')
bar(hax,temp1)
xlim([-0.5 0.5])
xlabel('Vertical Error (degrees)','fontsize',14)
figure
set(gca,'fontsize',14,'fontname','arial')
bar(hax,temp2)
xlim([-0.5 0.5])
xlabel('Vertical Error (degrees)','fontsize',14)


figure
set(gca,'fontsize',14,'fontname','arial')
subplot(2,1,1)
shadedErrorBar(used_fixs,post_maxX,tot_std_X,'k');
hold on
h = shadedErrorBar(used_fixs,X_mean,X_std,'b');
ylim([-0.5 0.5])
xlim([0 used_fixs(end)])
xlabel('Fixation Number','fontsize',14)
ylabel('Estimated Error (degrees)','fontsize',14)
title('Total Horizontal Error Estimate','fontsize',16)
subplot(2,1,2)
h = shadedErrorBar(used_fixs,post_maxX-X_mean',post_stdX,'k');
xlim([0 used_fixs(end)])
ylim([-0.5 0.5])
xlabel('Fixation Number','fontsize',14)
ylabel('Estimated Error (degrees)','fontsize',14)
title('High-frequency Horizontal Error Estimate','fontsize',16)


hax = linspace(-0.5,0.5,25);
temp1 = hist(post_maxX,hax);
temp2 = hist(post_maxX-X_mean',hax);
figure
set(gca,'fontsize',14,'fontname','arial')
bar(hax,temp1)
xlim([-0.5 0.5])
xlabel('Horizontal Error (degrees)','fontsize',14)
figure
set(gca,'fontsize',14,'fontname','arial')
bar(hax,temp2)
xlim([-0.5 0.5])
xlabel('Horizontal Error (degrees)','fontsize',14)

%%
close all

figure
set(gca,'fontsize',14,'fontname','arial')
plot(used_fixs,used_disp(:,2),'k')
hold on
plot(used_fixs,est_error_left(:,2),'r')
plot(used_fixs,est_error_right(:,2),'b')
ylim([-0.4 0.6])
xlim([0 used_fixs(end)])
xlabel('Fixation Number','fontsize',14)
ylabel('Estimated Error (degrees)','fontsize',14)
legend('Overall Error','Left Eye','Right Eye')
title('Vertical Error','fontsize',16)

hax = linspace(-0.6,0.6,50);
temp1 = hist(used_disp(:,2),hax);
temp2 = hist(est_error_right(:,2),hax);
temp3 = hist(est_error_left(:,2),hax);
figure
set(gca,'fontsize',14,'fontname','arial')
stairs(hax,temp1,'k')
hold on
stairs(hax,temp2,'b')
stairs(hax,temp3,'r')
xlim([-0.4 0.6])


figure
set(gca,'fontsize',14,'fontname','arial')
plot(used_fixs,used_disp(:,1),'k')
hold on
plot(used_fixs,est_error_left(:,1),'r')
plot(used_fixs,est_error_right(:,1),'b')
ylim([-0.6 0.6])
xlim([0 used_fixs(end)])
xlabel('Fixation Number','fontsize',14)
ylabel('Estimated Error (degrees)','fontsize',14)
legend('Overall Error','Left Eye','Right Eye')
title('Horizontal Error','fontsize',16)

hax = linspace(-0.6,0.6,50);
temp1 = hist(used_disp(:,1),hax);
temp2 = hist(est_error_right(:,1),hax);
temp3 = hist(est_error_left(:,1),hax);
figure
set(gca,'fontsize',14,'fontname','arial')
stairs(hax,temp1,'k')
hold on
stairs(hax,temp2,'b')
stairs(hax,temp3,'r')
xlim([-0.6 0.6])

%%
disp('Final model update');
SDIM = size(X,2);
Xmat_resh = reshape(X,size(X,1),SDIM,SDIM);
for i = 1:size(X,1)
    Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-dX(post_max(i)),2,0);
    Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-dY(post_max(i)),1,0);
end
Xmatr = reshape(Xmat_resh,NT,SDIM^2);
for c = 1:n_used_cells;
    fprintf('Fitting cell %d of %d\n',c,n_used_cells);
    hold_const = [0 0 0 0 0 0 0 0 0 0];
    [gabor_params_fin(c,:),LL_fin(c)] = fit_gabor_params(gabor_params{end}(c,:),Xmatr,spk_cnts(used_fixs,c),[SDIM SDIM],hold_const,LB,UB);
    
    cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
    mask1_out = Xmat*cur_mask1(:);
    mask2_out = Xmat*cur_mask2(:);
    e_out = sqrt(mask1_out.^2 + mask2_out.^2);
    gmat = gabor_params_fin(c,7)*e_out + gabor_params_fin(c,8)*mask1_out + gabor_params_fin(c,9)*mask2_out;
    [fin_scale(c,:),fin_offset(c,:),fin_LL(c)] = fit_gnm_gabormodel_blocks_v3(gmat,spk_cnts(used_fixs,c),blockids(used_fixs));
end

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
save sing_eye_pgabortrack_fin_est_d2 used_fixs dX dY gabor_params_fin LL_fin tot_* temp_scale fin_*

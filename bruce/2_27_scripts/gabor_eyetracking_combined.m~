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
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
% load fixation_stim_righteye
load fixation_data
% load fixation_image_patches
load fixation_image_patches
% X = X_right;
% X = X_left;
X_eye{1} = X_left;
X_eye{2} = X_right;

%% COMPUTE GABOR TEMPLATES
SDIM = size(X_left,2);

kern_len = SDIM-1;
[Xg,Yg] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);

cellids = [1 2 3 4 5 6 7 8 10];
hor_posd = [4.5 4.72 4.5 4.65 4.65 4.65 4.88 4.57 4.57];
vert_posd = [-2.61 -2 -1.93 -2.08 -1.7 -1.85 -1.7 -1.85 -1.63];
orientations = [0 0.31 0.31 0 0 1.9 pi/2 1.9 0];
lambdas = [1.0 0.35 0.35 0.5 0.67 0.35 0.67 0.35 0.67]*Fsd;
phases = [0 pi/2];
for i = 1:length(cellids)
    [~,vert_ind(i)] = min(abs(yax(ypatch_inds) - vert_posd(i)));
    [~,hor_ind(i)] = min(abs(xax(xpatch_inds) - hor_posd(i)));
end
vert_pos = Yg(vert_ind,1);
hor_pos = Xg(1,hor_ind);
n_used_sus = length(cellids);
gabor_bank = zeros(n_used_sus,length(phases),size(Xg,1),size(Xg,2));
for n = 1:n_used_sus
    for k = 1:length(phases)
        gabor_bank(n,k,:,:) = get_gabor_template(Xg,Yg,hor_pos(n),vert_pos(n),orientations(n),lambdas(n),phases(k));
    end
end
n_used_cells = n_used_sus;
% stationary_units = [1 2 3 6 7 8 9 10];
% eliminate = find(~ismember(cellids,stationary_units));
% gabor_bank(eliminate,:,:,:) = [];
% cellids(eliminate) = [];

%% COMPUTE GABOR TEMPLATES
% muaids = [1 2 4 5 6 9 10 11 12 13 14];
muaids = [1 2 4 5 6 9 12 13];

hor_posd = [4.5 4.57 nan 4.65 4.65 2.91 nan nan 3.51 4.72 4.72 3.0 3.0 4.88];
vert_posd = [-2.0 -1.78 nan -2.0 -2.0 -1.85 nan nan -1.78 -1.78 -1.85 -1.78 -1.78 -1.4];
orientations = [0.39 0.39 nan 0.39 0 1.57 nan nan 1.57 1.57 0.79 1.57 1.57 0];
lambdas = [0.5 0.5 nan 0.36 0.36 1 nan nan 0.67 0.36 0.67 0.67 0.67 1]*Fsd;
hor_posd = hor_posd(muaids); vert_posd = vert_posd(muaids);
orientations = orientations(muaids); lambdas = lambdas(muaids);

phases = [0 pi/2];
for i = 1:length(muaids)
    [~,vert_ind(i)] = min(abs(yax(ypatch_inds) - vert_posd(i)));
    [~,hor_ind(i)] = min(abs(xax(xpatch_inds) - hor_posd(i)));
end
vert_pos = Yg(vert_ind,1);
hor_pos = Xg(1,hor_ind);
n_used_mus = length(muaids);
gabor_bank_mua = zeros(n_used_mus,length(phases),size(Xg,1),size(Xg,2));
for n = 1:n_used_mus
    for k = 1:length(phases)
        gabor_bank_mua(n,k,:,:) = get_gabor_template(Xg,Yg,hor_pos(n),vert_pos(n),orientations(n),lambdas(n),phases(k));
    end
end

% stationary_mus = [1 4 5 8 10];
% eliminate = find(~ismember(muaids,stationary_mus));
% gabor_bank_mua(eliminate,:,:,:) = [];
% muaids(eliminate) = [];

%% treat MUs as SUs
gabor_bank = cat(1,gabor_bank,gabor_bank_mua);
spk_cnts = cat(2,spk_cnts(:,cellids),mua_cnts(:,muaids));
n_used_cells = size(spk_cnts,2);


for eye = 1:2
    %%
    X = X_eye{eye};
    %% CONVOLVE GABOR TEMPLATES WITH IMAGE PATCH
    NT = size(X,1);
    SDIM = size(X,2);
    deltax_range = 10;
    deltay_range = 10;
    delta_x = -deltax_range:deltax_range;
    delta_y = -deltay_range:deltay_range;
    use_x = find(ismember(Xg(1,:),delta_x)); use_x = use_x(:);
    use_y = find(ismember(Yg(:,1),delta_y)); use_y = use_y(:);
    
    resh_X = reshape(X,NT,SDIM^2);
    X_gabor_output = zeros(n_used_cells,2,NT,length(delta_x),length(delta_y));
    for n = 1:n_used_cells
        fprintf('Cell %d of %d\n',n,n_used_cells);
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
    
    %% FIT INITIAL MODELS
    [dX,dY] = meshgrid(delta_x,delta_y);
    dX = dX(:);
    dY = dY(:);
    K = length(dX);
    X_gabor_output = reshape(X_gabor_output,[n_used_cells length(phases) NT K]);
    
    zero_state = find(dX==0 & dY==0);
    cur_state_seq = ones(NT,1)*zero_state;
    cur_seq_dX = dX(cur_state_seq);
    cur_seq_dY = dY(cur_state_seq);
    
    % MODEL OUTPUT FOR EACH PERTERBATION
    dist_mat = squareform(pdist([dX dY]));
    sigma = 0.25;
    A = exp(-dist_mat.^2/(2*sigma.^2));
    for i = 1:K
        A(i,:) = A(i,:)/sum(A(i,:));
    end
    
    %%
    n_iter = 5;
    
    simple_outs = zeros(n_used_cells,length(phases),NT);
    for i = 1:NT
        simple_outs(:,:,i) = squeeze(X_gabor_output(:,:,i,cur_state_seq(i)));
    end
    complex_outs = squeeze(sqrt(simple_outs(:,1,:).^2+simple_outs(:,2,:).^2));
    
    for n = 1:n_used_cells;
        Robs = spk_cnts(used_fixs,n);
        simple_out = squeeze(simple_outs(n,:,:))';
        complex_out = complex_outs(n,:)';
        %     gabor_model(1,n) = fit_gnm_gabormodel([simple_out complex_out],Robs);
        gabor_model(1,n) = fit_gnm_gabormodel_blocks([simple_out complex_out],Robs,blockids(used_fixs));
        
    end
    
    LLs = zeros(n_used_cells,NT,K);
    for n = 1:n_used_cells
        LLs(n,:,:) = get_gabor_blocks_LLs(gabor_model(1,n),spk_cnts(used_fixs,n),squeeze(X_gabor_output(n,:,:,:)),...
            squeeze(sqrt(X_gabor_output(n,1,:,:).^2 + X_gabor_output(n,2,:,:).^2)),blockids(used_fixs));
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
            delta(t,:) = delta(t,:) + log(B(t,:));
        end
        
        % Backtracking Viterbi
        state_seq = zeros(NT,1);
        [llik_best,state_seq(NT)] = max(delta(NT,:));
        for t=NT-1:-1:1,
            state_seq(t) = psi(t+1,state_seq(t+1));
        end
        
        %%
        disp('Updating model');
        cur_X_seq(n,:) = dX(state_seq);
        cur_Y_seq(n,:) = dY(state_seq);
        
        %% UPDATE MODEL
        simple_outs = zeros(n_used_cells,length(phases),NT);
        for i = 1:NT
            simple_outs(:,:,i) = squeeze(X_gabor_output(:,:,i,state_seq(i)));
        end
        complex_outs = squeeze(sqrt(simple_outs(:,1,:).^2+simple_outs(:,2,:).^2));
        
        for c = 1:n_used_cells;
            Robs = spk_cnts(used_fixs,c);
            simple_out = squeeze(simple_outs(c,:,:))';
            complex_out = complex_outs(c,:)';
            %         gabor_model(n,c) = fit_gnm_gabormodel([simple_out complex_out],Robs);
            gabor_model(n,c) = fit_gnm_gabormodel_blocks([simple_out complex_out],Robs,blockids(used_fixs));
        end
        
        LLs = zeros(n_used_cells,NT,K);
        for c = 1:n_used_cells
            LLs(c,:,:) = get_gabor_blocks_LLs(gabor_model(n,c),spk_cnts(used_fixs,c),squeeze(X_gabor_output(c,:,:,:)),...
                squeeze(sqrt(X_gabor_output(c,1,:,:).^2 + X_gabor_output(c,2,:,:).^2)),blockids(used_fixs));
        end
        
    end
    
    disp('DONE!');
    %%
    cd /Users/James/Data/bruce/2_27_12/stimrecon
    if eye == 1
        save left_eye_gabortrack_combined_v2 loglik cur_*_seq dX dY post_std* post_entropy gabor_model
    else
        save right_eye_gabortrack_combined_v2 loglik cur_*_seq dX dY post_std* post_entropy gabor_model
    end
end
%%
cd /Users/James/Data/bruce/2_27_12/stimrecon

load left_eye_gabortrack_combined_v2
left_X_mean = cur_X_seq(end,:)/Fsd;
left_X_std = post_stdX(end,:)/Fsd;
left_Y_mean = cur_Y_seq(end,:)/Fsd;
left_Y_std = post_stdY(end,:)/Fsd;
% left_loc_i = best_loc_i;
% left_loc_j = best_loc_j;
load right_eye_gabortrack_combined_v2
right_X_mean = cur_X_seq(end,:)/Fsd;
right_X_std = post_stdX(end,:)/Fsd;
right_Y_mean = cur_Y_seq(end,:)/Fsd;
right_Y_std = post_stdY(end,:)/Fsd;
% right_loc_i = best_loc_i;
% right_loc_j = best_loc_j;

figure
H=shadedErrorBar(1:length(left_X_mean),left_X_mean,left_X_std,'b',1);
hold on
H=shadedErrorBar(1:length(right_X_mean),right_X_mean,right_X_std,'r',1);
axis tight
ylim([-10 10]/Fsd)
ylabel('Horizontal error (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)

figure
H=shadedErrorBar(1:length(left_Y_mean),left_Y_mean,left_Y_std,'b',1);
hold on
H=shadedErrorBar(1:length(right_Y_mean),right_Y_mean,right_Y_std,'r',1);
axis tight
ylim([-10 10]/Fsd)
ylabel('Vertical error (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)
hold on

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
observed_ydisp = all_sac_locs_right(used_fixs,2)-all_sac_locs_left(used_fixs,2);
smoothed_ydisp = smooth(observed_ydisp,20,'rlowess');

% inferred_ydisp = right_Y_mean - left_Y_mean+avg_ydisp;
inferred_ydisp = right_Y_mean - left_Y_mean;
inferred_ystd = 1/sqrt(2)*(right_Y_std+left_Y_std);

figure
hold on
plot(-observed_ydisp,'r','linewidth',0.5)
% plot(-smoothed_ydisp,'k','linewidth',2)
H=shadedErrorBar(1:length(left_Y_mean),inferred_ydisp,inferred_ystd,'b',0);
axis tight
ylim([-30 30]/Fsd)
xl = xlim();
line(xl,[0 0],'color','k')
ylabel('Vertical disparity (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)

observed_xdisp = all_sac_locs_right(used_fixs,1)-all_sac_locs_left(used_fixs,1);
smoothed_xdisp = smooth(observed_xdisp,20,'rlowess');

% inferred_xdisp = right_X_mean - left_X_mean+avg_xdisp;
inferred_xdisp = right_X_mean - left_X_mean;
inferred_xstd = 1/sqrt(2)*(right_X_std+left_X_std);

figure
hold on
plot(-observed_xdisp,'r','linewidth',0.5)
% plot(-smoothed_xdisp,'k','linewidth',2)
H=shadedErrorBar(1:length(left_X_mean),inferred_xdisp,inferred_xstd,'b',0);
axis tight
ylim([-30 30]/Fsd)
xl = xlim();
line(xl,[0 0],'color','k')
ylabel('Horizontal disparity (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)


%%

simple_outs = zeros(n_used_cells,length(phases),NT);
for i = 1:NT
    simple_outs(:,:,i) = squeeze(X_gabor_output(:,:,i,state_seq(i)));
end
complex_outs = squeeze(sqrt(simple_outs(:,1,:).^2+simple_outs(:,2,:).^2));

for c = 1:n_used_cells;
    Robs = spk_cnts(used_fixs,c);
    simple_out = squeeze(simple_outs(c,:,:))';
    complex_out = complex_outs(c,:)';
    %         gabor_model(n,c) = fit_gnm_gabormodel([simple_out complex_out],Robs);
    gabor_model(n,c) = fit_gnm_gabormodel_blocks([simple_out complex_out],Robs,blockids(used_fixs));
end

LLs = zeros(n_used_cells,NT,K);
for c = 1:n_used_cells
    LLs(c,:,:) = get_gabor_blocks_LLs(gabor_model(n,c),spk_cnts(used_fixs,c),squeeze(X_gabor_output(c,:,:,:)),...
        squeeze(sqrt(X_gabor_output(c,1,:,:).^2 + X_gabor_output(c,2,:,:).^2)),blockids(used_fixs));
end

%%


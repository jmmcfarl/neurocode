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

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
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
% load fixation_imagXe_patches
load fixation_image_patches_d2
% X = X_right;
% X = X_left;
X_eye{1} = X_left;
X_eye{2} = X_right;

for eye = 1:2
    X = X_eye{eye};
    %% COMPUTE GABOR TEMPLATES
    NT = size(X,1);
    SDIM = size(X,2);
    kern_len = SDIM-1;
    [Xg,Yg] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);
    
    Xmat = reshape(X,NT,SDIM^2);
    
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
    n_used_cells = length(cellids);
    
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
        [gabor_params{1}(n,:),LL{1}(n)] = fit_gabor_params(init_params,Xmat,spk_cnts(used_fixs,cellids(n)),[SDIM SDIM],hold_const);
    end
    
    %%
%     n = 6;
% %     close all
%     cur_mask_init = get_pgabor_mask(ip(n,1:6),0,[SDIM SDIM]);
%     cur_mask_fin = get_pgabor_mask(gabor_params{1}(n,1:6),0,[SDIM SDIM]);
%     subplot(2,1,1)
%     imagesc(cur_mask_init)
%     subplot(2,1,2)
%     imagesc(cur_mask_fin);
%     
    
    %% MODEL OUTPUT FOR EACH PERTURBATION
    deltax_range = 20;
    deltay_range = 20;
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
    
    %%
    X_gabor_output = reshape(X_gabor_output,[n_used_cells 2 NT K]);
    
    zero_state = find(dX==0 & dY==0);
    cur_state_seq = ones(NT,1)*zero_state;
    cur_seq_dX = dX(cur_state_seq);
    cur_seq_dY = dY(cur_state_seq);
    
    dist_mat = squareform(pdist([dX dY]));
    sigma = 0.25;
    A = exp(-dist_mat.^2/(2*sigma.^2));
    for i = 1:K
        A(i,:) = A(i,:)/sum(A(i,:));
    end
    
    %%
    n_iter = 30;
    
    LLs = zeros(n_used_cells,NT,K);
    for n = 1:n_used_cells
        LLs(n,:,:) = get_pgabor_LLs(gabor_params{1}(n,:),spk_cnts(used_fixs,cellids(n)),squeeze(X_gabor_output(n,:,:,:,:)));
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
            %         hold_const = [1 1 1 1 1 0 0 1 0];
            hold_const = [0 0 0 0 0 0 0 0 0 0];
            [gabor_params{n}(c,:),LL{n}(c)] = fit_gabor_params(gabor_params{n-1}(c,:),Xmatr,spk_cnts(used_fixs,cellids(c)),[SDIM SDIM],hold_const);
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
                LLs(c,:,:) = get_pgabor_LLs(gabor_params{n}(c,:),spk_cnts(used_fixs,cellids(c)),squeeze(X_gabor_output(c,:,:,:,:)));
            end
        end
    end
    
    %%
    cd /Users/James/Data/bruce/2_27_12/stimrecon
    
    if eye == 1
        save left_eye_pgabortrack_d2v1 loglik cur_*_seq dX dY post_std* post_entropy gabor_params LL
    else
        save right_eye_pgabortrack_d2v1 loglik cur_*_seq dX dY post_std* post_entropy gabor_params LL
        
    end
end
%%
    cd /Users/James/Data/bruce/2_27_12/stimrecon
close all
load left_eye_pgabortrack_d2v1
left_X_mean = cur_X_seq(end,:)/Fsd;
left_X_std = post_stdX(end,:)/Fsd;
left_Y_mean = cur_Y_seq(end,:)/Fsd;
left_Y_std = post_stdY(end,:)/Fsd;
left_x_pos = gabor_params{4}(:,1);
left_y_pos = gabor_params{4}(:,2);

load right_eye_pgabortrack_d2v1
n_iter = size(cur_X_seq,1);
right_X_mean = cur_X_seq(end,:)/Fsd;
right_X_std = post_stdX(end,:)/Fsd;
right_Y_mean = cur_Y_seq(end,:)/Fsd;
right_Y_std = post_stdY(end,:)/Fsd;
right_x_pos = gabor_params{n_iter}(:,1);
right_y_pos = gabor_params{n_iter}(:,2);

avg_x_pos = 0.5*left_x_pos + 0.5*right_x_pos;
left_xoffset = mean(left_x_pos - avg_x_pos);
right_xoffset = mean(right_x_pos - avg_x_pos);
% left_xoffset = -2.3;
% right_xoffset = 2.3;
left_X_mean = left_X_mean + left_xoffset/Fsd;
right_X_mean = right_X_mean + right_xoffset/Fsd;

avg_y_pos = 0.5*left_y_pos + 0.5*right_y_pos;
left_yoffset = mean(left_y_pos - avg_y_pos);
right_yoffset = mean(right_y_pos - avg_y_pos);
% left_yoffset = 3;
% right_yoffset = -3;
left_Y_mean = left_Y_mean + left_yoffset/Fsd;
right_Y_mean = right_Y_mean + right_yoffset/Fsd;

t1 = find(blockids(used_fixs) > 1,1,'first');
t2 = find(blockids(used_fixs) > 2,1,'first');
t3 = find(blockids(used_fixs) > 3,1,'first');

figure
set(gca,'fontsize',12)
H=shadedErrorBar(1:length(left_X_mean),left_X_mean,left_X_std,'b',1);
hold on
H=shadedErrorBar(1:length(right_X_mean),right_X_mean,right_X_std,'r',1);
axis tight
% ylim([-10 10]/Fsd)
ylim([-20 20]/Fsd)
yl = ylim();
line([t1 t1],yl,'color','k');
line([t2 t2],yl,'color','k');
line([t3 t3],yl,'color','k');
ylabel('Horizontal error (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)

figure
set(gca,'fontsize',12)
H=shadedErrorBar(1:length(left_Y_mean),left_Y_mean,left_Y_std,'b',1);
hold on
H=shadedErrorBar(1:length(right_Y_mean),right_Y_mean,right_Y_std,'r',1);
axis tight
% ylim([-10 10]/Fsd)
ylim([-20 20]/Fsd)
yl = ylim();
line([t1 t1],yl,'color','k');
line([t2 t2],yl,'color','k');
line([t3 t3],yl,'color','k');
ylabel('Vertical error (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)



observed_ydisp = all_sac_locs_right(used_fixs,2)-all_sac_locs_left(used_fixs,2);
smoothed_ydisp = smooth(observed_ydisp,20,'rlowess');

inferred_ydisp = right_Y_mean - left_Y_mean;
inferred_ystd = 1/sqrt(2)*(right_Y_std+left_Y_std);

figure
set(gca,'fontsize',12)
hold on
plot(-observed_ydisp,'r','linewidth',0.5)
H=shadedErrorBar(1:length(left_Y_mean),inferred_ydisp,inferred_ystd,'b',0);
axis tight
% ylim([-15 15]/Fsd)
ylim([-26 5]/Fsd)
yl = ylim();
line([t1 t1],yl,'color','k');
line([t2 t2],yl,'color','k');
line([t3 t3],yl,'color','k');
xl = xlim();
line(xl,[0 0],'color','k','linestyle','--','linewidth',2)
ylabel('Vertical disparity (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)

observed_xdisp = all_sac_locs_right(used_fixs,1)-all_sac_locs_left(used_fixs,1);
smoothed_xdisp = smooth(observed_xdisp,20,'rlowess');

inferred_xdisp = right_X_mean - left_X_mean;
inferred_xstd = 1/sqrt(2)*(right_X_std+left_X_std);

figure
set(gca,'fontsize',12)
hold on
plot(-observed_xdisp,'r','linewidth',0.5)
H=shadedErrorBar(1:length(left_X_mean),inferred_xdisp,inferred_xstd,'b',0);
axis tight
% ylim([-15 15]/Fsd)
ylim([-15 30]/Fsd)
xl = xlim();
yl = ylim();
line([t1 t1],yl,'color','k');
line([t2 t2],yl,'color','k');
line([t3 t3],yl,'color','k');
line(xl,[0 0],'color','k','linestyle','--','linewidth',2)
ylabel('Horizontal disparity (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)
set(gca,'ytick',-0.5:0.2:1);

%%
est_left_X = all_sac_locs_left(used_fixs,1) + left_X_mean';
est_left_Y = all_sac_locs_left(used_fixs,2) + left_Y_mean';

est_right_X = all_sac_locs_right(used_fixs,1) + right_X_mean';
est_right_Y = all_sac_locs_right(used_fixs,2) + right_Y_mean';

est_disparity_x = est_right_X - est_left_X;
est_disparity_y = est_right_Y - est_left_Y;



%%
close all
    cd /Users/James/Data/bruce/2_27_12/stimrecon
load left_eye_pgabortrack_d2v1
cd ~/James_scripts/bruce
lgabor_params = gabor_params;
for n = 1:n_used_cells
    
    cur_gabor_params = gabor_params{end}(n,:);
    comp_ratio(n) = log(cur_gabor_params(7)/norm(cur_gabor_params(8:9)));
    
%     cur_gabor_params(1:5) = ip(n,:);
    cur_mask1 = get_pgabor_mask(cur_gabor_params(1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(cur_gabor_params(1:6),pi/2,[SDIM SDIM]);
    
    cur_mask = cur_gabor_params(8)*cur_mask1 + cur_gabor_params(9)*cur_mask2;
    subplot(3,3,n)
    imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_mask); colormap(gray);
    axis square
    set(gca,'ydir','normal');
    set(gca,'fontsize',12)
    xlabel('Horizontal position (degrees)','fontsize',14)
    ylabel('Vertical position (degrees)','fontsize',14)
    title(sprintf('Cell%d, CI: %.1f PO: %.1f',cellids(n),comp_ratio(n),180/pi*cur_gabor_params(3)),'fontsize',14);
%     fname = sprintf('Fit_Gabor_filter_cell%d',cellids(n));
%     print('-dpdf',fname);
%     close 
end
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [11 11]);

figure
cd /Users/James/Data/bruce/2_27_12/stimrecon
load right_eye_pgabortrack_d2v1
cd ~/James_scripts/bruce
rgabor_params = gabor_params;
for n = 1:n_used_cells
    
    cur_gabor_params = gabor_params{end}(n,:);
    comp_ratio(n) = log(cur_gabor_params(7)/norm(cur_gabor_params(8:9)));
    
%     cur_gabor_params(1:5) = ip(n,:);
    cur_mask1 = get_pgabor_mask(cur_gabor_params(1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(cur_gabor_params(1:6),pi/2,[SDIM SDIM]);
    
    cur_mask = cur_gabor_params(8)*cur_mask1 + cur_gabor_params(9)*cur_mask2;
    subplot(3,3,n)
    imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_mask); colormap(gray);
    axis square
    set(gca,'ydir','normal');
    set(gca,'fontsize',12)
    xlabel('Horizontal position (degrees)','fontsize',14)
    ylabel('Vertical position (degrees)','fontsize',14)
    title(sprintf('Cell%d, CI: %.1f PO: %.1f',cellids(n),comp_ratio(n),180/pi*cur_gabor_params(3)),'fontsize',14);
%     fname = sprintf('Fit_Gabor_filter_cell%d',cellids(n));
%     print('-dpdf',fname);
%     close 
end
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [11 11]);

%%
load left_eye_pgabortrack_d2v1
llam_mat = zeros(n_used_cells,n_iter);
lsig_mat = zeros(n_used_cells,n_iter);
lecc_mat = zeros(n_used_cells,n_iter);
lwc_mat = zeros(n_used_cells,n_iter);
lws_mat = zeros(n_used_cells,n_iter);
for i = 1:n_iter
    llam_mat(:,i) = gabor_params{i}(:,4);
    lsig_mat(:,i) = gabor_params{i}(:,5);
    lecc_mat(:,i) = gabor_params{i}(:,6);
    lwc_mat(:,i) = gabor_params{i}(:,7);
    lws_mat(:,i) = max(abs(gabor_params{i}(:,8:9)),[],2);

end
load right_eye_pgabortrack_d2v1
rlam_mat = zeros(n_used_cells,n_iter);
rsig_mat = zeros(n_used_cells,n_iter);
recc_mat = zeros(n_used_cells,n_iter);
rwc_mat = zeros(n_used_cells,n_iter);
rws_mat = zeros(n_used_cells,n_iter);
for i = 1:n_iter
    rlam_mat(:,i) = gabor_params{i}(:,4);
    rsig_mat(:,i) = gabor_params{i}(:,5);
    recc_mat(:,i) = gabor_params{i}(:,6);
    rwc_mat(:,i) = gabor_params{i}(:,7);
    rws_mat(:,i) = max(abs(gabor_params{i}(:,8:9)),[],2);
end

figure
subplot(2,2,1)
plot(mean(llam_mat))
hold on
plot(mean(rlam_mat),'r')

subplot(2,2,2)
plot(mean(lsig_mat))
hold on
plot(mean(rsig_mat),'r')

subplot(2,2,3)
plot(mean(lwc_mat))
hold on
plot(mean(rwc_mat),'r')
subplot(2,2,4)
plot(mean(lws_mat))
hold on
plot(mean(rws_mat),'r')

%%
% load right_eye_pgabortrack_d2v1
% rgabor_params = gabor_params;
% load left_eye_pgabortrack_d2v1
% lgabor_params = gabor_params;
% 
% n = 1;
% %     close all
% cur_lmask_init = get_pgabor_mask(lgabor_params{1}(n,1:6),0,[SDIM SDIM]);
% cur_lmask_fin = get_pgabor_mask(lgabor_params{end}(n,1:6),0,[SDIM SDIM]);
% cur_rmask_init = get_pgabor_mask(rgabor_params{1}(n,1:6),0,[SDIM SDIM]);
% cur_rmask_fin = get_pgabor_mask(rgabor_params{end}(n,1:6),0,[SDIM SDIM]);
% subplot(2,2,1)
% imagesc(cur_lmask_init)
% subplot(2,2,2)
% imagesc(cur_lmask_fin);
% subplot(2,2,3)
% imagesc(cur_rmask_init)
% subplot(2,2,4)
% imagesc(cur_rmask_fin);


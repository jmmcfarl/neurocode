clear all
% close all

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

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_stim_lefteye

%% COMPUTE GABOR TEMPLATES
Pix2Deg = 0.018837;
% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

kern_len = 33;
[Xg,Yg] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);

cellids = [10 2 3 4 5 6 8];
n_used_cells = length(cellids);
orientations = [0 0.39 0.39 0 1.18 1.96 1.96];
% lambdas = 1./[0.5 1 1.5 2 2.5 3]*Fsd;
lambdas = 1./[1.5 2.5 2.5 2.5 1.5 2.5 2.5]*Fsd;
phases = [0 pi/2];
gabor_bank = zeros(n_used_cells,length(phases),size(Xg,1),size(Xg,2));
for n = 1:n_used_cells
    for k = 1:length(phases)
        gabor_bank(n,k,:,:) = get_gabor_template(Xg,Yg,0,0,orientations(n),lambdas(n),phases(k));
    end
end

%% CONVOLVE GABOR TEMPLATES WITH IMAGE PATCH
NT = size(X,1);
SDIM = size(X,2);

X_gabor_energy = zeros([n_used_cells size(X)]);
for i = 1:NT
    for n = 1:n_used_cells
        gabor_phase1 = conv2(squeeze(X(i,:,:)),squeeze(gabor_bank(n,1,:,:)),'same');
        gabor_phase2 = conv2(squeeze(X(i,:,:)),squeeze(gabor_bank(n,2,:,:)),'same');
        X_gabor_energy(n,i,:,:) = sqrt(gabor_phase1.^2+gabor_phase2.^2);
    end
end

X_gabor_renergy = reshape(X_gabor_energy,[n_used_cells,NT,SDIM^2]);


%% FIND RF CENTER FOR EACH CELL
Bmat = zeros(n_used_cells,SDIM,SDIM,2);
R2vals = zeros(n_used_cells,SDIM,SDIM,1);
for n = 1:n_used_cells
    for i = 1:SDIM
        for j = 1:SDIM
            [b,bint,r,rint,stats] = regress(spk_cnts(:,cellids(n)),...
                [ones(NT,1) squeeze(X_gabor_energy(n,:,i,j))']);
            Bmat(n,i,j,:) = b;
            R2vals(n,i,j) = stats(1);
        end
    end
    temp = squeeze(R2vals(n,:,:));
    [best_R2(n),best_loc(n)] = max(temp(:));
    [best_loc_i(n),best_loc_j(n)] = ind2sub(size(temp),best_loc(n));
end


%% MODEL OUTPUT FOR EACH PERTERBATION
deltax_range = 10;
deltay_range = 10;
delta_x = -deltax_range:deltax_range;
delta_y = -deltay_range:deltay_range;
[dX,dY] = meshgrid(delta_x,delta_y);
dX = dX(:);
dY = dY(:);
K = length(dX);

dist_mat = squareform(pdist([dX dY]));
sigma = 0.25;
A = exp(-dist_mat.^2/(2*sigma.^2));
for i = 1:K
    A(i,:) = A(i,:)/sum(A(i,:));
end

relative_gabor_energy = zeros(n_used_cells,NT,K);
for n = 1:n_used_cells
    for i = 1:K
        if best_loc_i(n)+dY(i) < SDIM && best_loc_j(n)+dX(i) < SDIM
            relative_gabor_energy(n,:,i) = reshape(X_gabor_energy(n,:,best_loc_i(n)+dY(i),...
                best_loc_j(n)+dX(i)),NT,1);
        end
    end
end

%%
n_iter = 5;
clear cur_X_seq cur_Y_seq loglik cur_model cur_r2 post_entropy post_stdX post_stdY
for c = 1:n_used_cells
    cur_model(1,c,:) = Bmat(c,best_loc_i(c),best_loc_j(c),:);
    cur_r2(1,c) = best_R2(c);
end
for n = 2:n_iter
    fprintf('Iteration %d of %d\n',n,n_iter);
    
    disp('Computing Model Likelihoods');
    rpreds = zeros(size(relative_gabor_energy));
    LLs = zeros(size(relative_gabor_energy));
    for c = 1:n_used_cells
        rpreds(c,:,:) = cur_model(n-1,c,1) + cur_model(n-1,c,2)*relative_gabor_energy(c,:,:);
%         LLs(c,:,:) = repmat(spk_cnts(:,cellids(c)),1,K).*log(squeeze(rpreds(c,:,:))) - squeeze(rpreds(c,:,:));
        spk_cnt_mat = repmat(spk_cnts(:,cellids(c)),1,K);
        set1 = find(spk_cnts(:,cellids(c)) > 10); %use stirling's approx for log(R!)
        set2 = find(spk_cnts(:,cellids(c)) <= 10); %use log(R!) directly
        LLs(c,set1,:) = spk_cnt_mat(set1,:).*log(squeeze(rpreds(c,set1,:))) - squeeze(rpreds(c,set1,:)) - ...
            spk_cnt_mat(set1,:).*log(spk_cnt_mat(set1,:))+spk_cnt_mat(set1,:);
        LLs(c,set2,:) = spk_cnt_mat(set2,:).*log(squeeze(rpreds(c,set2,:))) - squeeze(rpreds(c,set2,:)) - ...
            log(factorial(spk_cnt_mat(set2,:)));
    end
    
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
    
    for c = 1:n_used_cells
        new_gabor_out = zeros(NT,1);
        for i = 1:NT
            if best_loc_i(c)+cur_Y_seq(n,i) < SDIM && best_loc_j(c)+cur_X_seq(n,i) < SDIM
                new_gabor_out(i) = X_gabor_energy(c,i,best_loc_i(c)+cur_Y_seq(n,i),best_loc_j(c)+cur_X_seq(n,i));
            end
        end
        [b,bint,r,rint,stats] = regress(spk_cnts(:,cellids(c)),...
            [ones(NT,1) new_gabor_out]);
        cur_model(n,c,:) = b;
        cur_r2(n,c) = stats(1);
    end    
end
disp('DONE!');
%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
save left_eye_gabortrack_sm cur_model cur_r2 loglik cur_*_seq best_loc_* dX dY post_std* post_entropy

%%

load left_eye_gabortrack_sm
left_X_mean = cur_X_seq(end,:)/Fsd;
left_X_std = post_stdX(end,:)/Fsd;
left_Y_mean = cur_Y_seq(end,:)/Fsd;
left_Y_std = post_stdY(end,:)/Fsd;
left_loc_i = best_loc_i;
left_loc_j = best_loc_j;
load right_eye_gabortrack_sm
right_X_mean = cur_X_seq(end,:)/Fsd;
right_X_std = post_stdX(end,:)/Fsd;
right_Y_mean = cur_Y_seq(end,:)/Fsd;
right_Y_std = post_stdY(end,:)/Fsd;
right_loc_i = best_loc_i;
right_loc_j = best_loc_j;

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
load fixation_stim_righteye

avg_ydisp = mean(right_loc_i-left_loc_i)/Fsd;
avg_xdisp = mean(right_loc_j-left_loc_j)/Fsd;

observed_ydisp = all_sac_locs_right(:,2)-all_sac_locs_left(:,2);
smoothed_ydisp = smooth(observed_ydisp,20,'rlowess');

inferred_ydisp = right_Y_mean - left_Y_mean+avg_ydisp;
inferred_ystd = 1/sqrt(2)*(right_Y_std+left_Y_std);

figure
hold on
plot(-observed_ydisp,'r','linewidth',0.5)
% plot(-smoothed_ydisp,'k','linewidth',2)
H=shadedErrorBar(1:length(left_Y_mean),inferred_ydisp,inferred_ystd,'b',0);
axis tight
ylim([-15 15]/Fsd)
xl = xlim();
line(xl,[0 0],'color','k')
ylabel('Vertical disparity (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)

observed_xdisp = all_sac_locs_right(:,1)-all_sac_locs_left(:,1);
smoothed_xdisp = smooth(observed_xdisp,20,'rlowess');

inferred_xdisp = right_X_mean - left_X_mean+avg_xdisp;
inferred_xstd = 1/sqrt(2)*(right_X_std+left_X_std);

figure
hold on
plot(-observed_xdisp,'r','linewidth',0.5)
% plot(-smoothed_xdisp,'k','linewidth',2)
H=shadedErrorBar(1:length(left_X_mean),inferred_xdisp,inferred_xstd,'b',0);
axis tight
ylim([-15 15]/Fsd)
xl = xlim();
line(xl,[0 0],'color','k')
ylabel('Horizontal disparity (degrees)','fontsize',14)
xlabel('Fixation number','fontsize',14)



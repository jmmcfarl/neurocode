clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12/G029/

load ./jbeG029.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G029Expts.mat
load ./eye_calibration_data

dt = 118/1e4;
frames_per_jump = 44;
ims_per_trial = 4;
dt_per_jump = frames_per_jump*dt;
fst = 1/dt;

Pix2Deg = 0.018837;
dsfrac = 1;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
RF_patch = [-0.2 1.; -1. 0.2]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [-1 2; -2. 1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));


%%
Expt_nu = [14 22 27]; %these are the grating expts
spk_win = 0.15; %in sec
n_allunits = 96;
min_trial_dur = 1.5+spk_win;
all_xo_vec = [];
all_yo_vec = [];
all_image_vec = [];
all_angles = [];
all_jump_sizes = [];
all_expt_vec = [];
all_trial_vec = [];
all_binned_spks = [];
all_stim_num = [];

single_units = find(CellList(Expt_nu(1),:,1) > 0);
multi_units = setdiff(1:n_allunits,single_units);
for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    n_sus = length(single_units);
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    %     completed_trials = 1:length(Trial_durs); %use all trials
    completed_trials = find(Trial_durs > min_trial_dur);
    Trial_im_nums = [Expts{Expt_nu(ee)}.Trials(:).se];
    Trial_angles = [Expts{Expt_nu(ee)}.Trials(:).Fa];
    jump_sizes = [Expts{Expt_nu(ee)}.Trials(:).sz];
    
    %     Trial_angles = Trial_angles + 90;
    
    for i = 1:length(completed_trials)
        cur_images = repmat(Trial_im_nums(completed_trials(i)),ims_per_trial,1);
        
        cur_xo = (0:(ims_per_trial-1))*jump_sizes(completed_trials(i))*cos(degtorad(Trial_angles(completed_trials(i))));
        cur_yo = (0:(ims_per_trial-1))*jump_sizes(completed_trials(i))*sin(degtorad(Trial_angles(completed_trials(i))));
        
        cur_t_starts = Trial_starts(completed_trials(i)) + (0:3)*dt_per_jump;
        cur_t_stops = cur_t_starts + spk_win;
        cur_t_edges = sort([cur_t_starts cur_t_stops]);
        
        cur_binned_spks = nan(n_allunits,length(cur_t_edges));
        for j = 1:n_allunits
            cur_binned_spks(j,:) = histc(Clusters{(j)}.times,cur_t_edges);
        end
        cur_binned_spks = cur_binned_spks(:,1:2:length(cur_t_edges));
        
        all_image_vec = [all_image_vec; cur_images];
        all_xo_vec = [all_xo_vec; cur_xo'];
        all_yo_vec = [all_yo_vec; cur_xo'];
        all_expt_vec = [all_expt_vec; ones(length(cur_images),1)*Expt_nu(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_images),1)*completed_trials(i)];
        all_binned_spks = [all_binned_spks; cur_binned_spks'];
        all_angles = [all_angles; repmat(Trial_angles(completed_trials(i)),ims_per_trial,1)];
        all_jump_sizes = [all_jump_sizes; repmat(jump_sizes(completed_trials(i)),ims_per_trial,1)];
        all_stim_num = [all_stim_num; (1:ims_per_trial)'];
    end
end

%%
% % all_image_vec = all_image_vec(randperm(length(all_image_vec)));
% cd ~/James_scripts/data_processing/Images/image_set_A
% tot_images = length(unique(all_image_vec));
% used_images = unique(all_image_vec);
% tot_samps = length(all_image_vec);
% all_im_patches = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
% for i = 1:tot_images
%     filename = sprintf('%.4d.png',used_images(i));
%     IMAGEorg = imread(filename);
%     IMAGEorg = double(IMAGEorg); % convert to double format
%     IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
%     IMAGE = flipud(IMAGEorg);
%     IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
%
%     cur_samp_set = find(all_image_vec == used_images(i));
%     fprintf('Analyzing image %d, %d samps\n',i,length(cur_samp_set));
%
%     for j = 1:length(cur_samp_set)
%         ypatch_inds_adj = round(ypatch_inds - all_yo_vec(cur_samp_set(j))*Fsd);
%         xpatch_inds_adj = round(xpatch_inds - all_xo_vec(cur_samp_set(j))*Fsd);
%         all_im_patches(cur_samp_set(j),:,:) = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
%     end
% end
%%
% all_image_vec = all_image_vec(randperm(length(all_image_vec)));

cd ~/James_scripts/data_processing/Images/image_set_A
tot_images = length(unique(all_image_vec));
used_images = unique(all_image_vec);
tot_samps = length(all_image_vec);
all_im_patches_pp = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
all_im_patches_pn = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
all_im_patches_np = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
all_im_patches_nn = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
for i = 1:tot_images
    filename = sprintf('%.4d.png',used_images(i));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGE = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
%     IMAGE = flipud(IMAGE);
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(all_image_vec == used_images(i));
    fprintf('Analyzing image %d, %d samps\n',i,length(cur_samp_set));
    
    for j = 1:length(cur_samp_set)
        ypatch_inds_adj = round(ypatch_inds - all_yo_vec(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds - all_xo_vec(cur_samp_set(j))*Fsd);
        all_im_patches_nn(cur_samp_set(j),:,:) = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        
        ypatch_inds_adj = round(ypatch_inds - all_yo_vec(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds + all_xo_vec(cur_samp_set(j))*Fsd);
        all_im_patches_pn(cur_samp_set(j),:,:) = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        
        ypatch_inds_adj = round(ypatch_inds + all_yo_vec(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds - all_xo_vec(cur_samp_set(j))*Fsd);
        all_im_patches_np(cur_samp_set(j),:,:) = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        
        ypatch_inds_adj = round(ypatch_inds + all_yo_vec(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds + all_xo_vec(cur_samp_set(j))*Fsd);
        all_im_patches_pp(cur_samp_set(j),:,:) = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        
        cur_im = all_im_patches_nn(cur_samp_set(j),:,:);
        all_im_patches_nn(cur_samp_set(j),:,:) = cur_im - mean(cur_im(:));
        cur_im = all_im_patches_np(cur_samp_set(j),:,:);
        all_im_patches_np(cur_samp_set(j),:,:) = cur_im - mean(cur_im(:));
        cur_im = all_im_patches_pn(cur_samp_set(j),:,:);
        all_im_patches_pn(cur_samp_set(j),:,:) = cur_im - mean(cur_im(:));
        cur_im = all_im_patches_pp(cur_samp_set(j),:,:);
        all_im_patches_pp(cur_samp_set(j),:,:) = cur_im - mean(cur_im(:));
    end
end

%%
% fullX = reshape(all_im_patches,size(all_im_patches,1),length(ypatch_inds)*length(xpatch_inds));
fullX_nn = reshape(all_im_patches_nn,size(all_im_patches_nn,1),length(ypatch_inds)*length(xpatch_inds));
fullX_np = reshape(all_im_patches_np,size(all_im_patches_nn,1),length(ypatch_inds)*length(xpatch_inds));
fullX_pn = reshape(all_im_patches_pn,size(all_im_patches_nn,1),length(ypatch_inds)*length(xpatch_inds));
fullX_pp = reshape(all_im_patches_pp,size(all_im_patches_nn,1),length(ypatch_inds)*length(xpatch_inds));

xv_frac = 0.2;
NT = size(fullX_nn,1);
xv_NT = round(xv_frac*NT);
xv_samp = randperm(NT);
% xv_samp = 1:NT;
xv_samp(xv_NT+1:end) = [];
tr_samp = setdiff(1:NT,xv_samp);

% tr_samp(all_stim_num(tr_samp)==1) = [];
% xv_samp(all_stim_num(xv_samp)==1) = [];
%%
sdim = 64;
bandwidth = 0.25;
ar = 1.5;
kern_len = 21-1;
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
orientations = linspace(0,pi,10);
orientations(end) = [];
% orientations = 0.9;
lambdas = [1/6];
% gab_xo = [0.35];
% gab_yo = [-0.48];
% gab_xo = [0.1:0.025:0.7];
% gab_yo = [-0.7:0.025:-0.1];
gab_xo = [-0.2:0.025:1.];
gab_yo = [-1:0.025:0.2];
gabor_bank = zeros(length(orientations)*length(lambdas)*length(gab_xo)*length(gab_yo),size(XX,1),size(XX,2));
gabor_bank2 = zeros(length(orientations)*length(lambdas)*length(gab_xo)*length(gab_yo),size(XX,1),size(XX,2));
cnt = 1;
clear gabor_props
for i = 1:length(orientations)
    for j = 1:length(lambdas)
        for m = 1:length(gab_xo)
            for l = 1:length(gab_yo)
                gabor_bank(cnt,:,:) = get_gabor_template(XX,YY,gab_xo(l),gab_yo(m),orientations(i),lambdas(j),0,bandwidth,ar);
                gabor_props(cnt).yo = gab_yo(l);,fullX_x(:,2100:3000))
                gabor_props(cnt).xo = gab_xo(m);
                gabor_props(cnt).lambda = lambdas(j);
                gabor_props(cnt).ori = orientations(i);
                gabor_props(cnt).bandwidth = bandwidth;
                gabor_props(cnt).ar = ar;
                
                gabor_bank2(cnt,:,:) = get_gabor_template(XX,YY,gab_xo(l),gab_yo(m),orientations(i),lambdas(j),pi/2,bandwidth,ar);
                cnt = cnt + 1;
            end
        end
    end
end
gabor_filts = reshape(gabor_bank,size(gabor_bank,1),sdim^2);
gabor_filts2 = reshape(gabor_bank2,size(gabor_bank,1),sdim^2);

%%

% gabor_out_X1 = fullX*gabor_filts';
% gabor_out_X2 = fullX*gabor_filts2';
% gabor_out_X = gabor_out_X1.^2 + gabor_out_X2.^2;
% gabor_out_X = sqrt(gabor_out_X);
% gabor_out_X = zscore(gabor_out_X);
% n_gabor_bases = size(gabor_out_X,2);
% oris = [gabor_props(:).ori];
% gab_lambdas = [gabor_props(:).lambda];
% XO = [gabor_props(:).xo];
% YO = [gabor_props(:).yo];

gabor_out_X1 = fullX_nn*gabor_filts';
gabor_out_X2 = fullX_nn*gabor_filts2';
gabor_out_X_nn = gabor_out_X1.^2 + gabor_out_X2.^2;
gabor_out_X_nn = sqrt(gabor_out_X_nn);
gabor_out_X_nn = zscore(gabor_out_X_nn);

gabor_out_X1 = fullX_np*gabor_filts';
gabor_out_X2 = fullX_np*gabor_filts2';
gabor_out_X_np = gabor_out_X1.^2 + gabor_out_X2.^2;
gabor_out_X_np = sqrt(gabor_out_X_np);
gabor_out_X_np = zscore(gabor_out_X_np);

gabor_out_X1 = fullX_pn*gabor_filts';
gabor_out_X2 = fullX_pn*gabor_filts2';
gabor_out_X_pn = gabor_out_X1.^2 + gabor_out_X2.^2;
gabor_out_X_pn = sqrt(gabor_out_X_pn);
gabor_out_X_pn = zscore(gabor_out_X_pn);

gabor_out_X1 = fullX_pp*gabor_filts';
gabor_out_X2 = fullX_pp*gabor_filts2';
gabor_out_X_pp = gabor_out_X1.^2 + gabor_out_X2.^2;
gabor_out_X_pp = sqrt(gabor_out_X_pp);
gabor_out_X_pp = zscore(gabor_out_X_pp);

n_gabor_bases = size(gabor_out_X1,2);
oris = [gabor_props(:).ori];
gab_lambdas = [gabor_props(:).lambda];
XO = [gabor_props(:).xo];
YO = [gabor_props(:).yo];

%%


cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
% load ./onednoise_rffits
% load ./oned_noise_mapping
load ./oned_fixation_fits

disp = 0;

clear sta_gabor xv_LL*
% for t = 1:15;
for t = 1:96
    fprintf('MU %d of %d\n',t,96);
    cur_binned_spks = all_binned_spks(tr_samp,t);
    % cur_binned_spks = all_binned_spks(:,single_units{1}(t));
    unique_spk_cnts = unique(cur_binned_spks);
    spikebins = [];
    for i = 2:length(unique_spk_cnts)
        cur_set = find(cur_binned_spks == unique_spk_cnts(i));
        spikebins = [spikebins; repmat(cur_set,unique_spk_cnts(i),1)];
    end
    spikebins = sort(spikebins);
    
    % %     sta_gabor(t,:) = mean(gabor_out_X(tr_samp(spikebins),:));
    % %     cur_gabor_resh = reshape(sta_gabor(t,:),[length(gab_xo)*length(gab_yo),length(lambdas),length(orientations)]);
    %     sta_gabor_nn(t,:) = mean(gabor_out_X_nn(tr_samp(spikebins),:));
    %     cur_gabor_resh_nn = reshape(sta_gabor_nn(t,:),[length(gab_xo)*length(gab_yo),length(lambdas),length(orientations)]);
    %     sta_gabor_np(t,:) = mean(gabor_out_X_np(tr_samp(spikebins),:));
    %     cur_gabor_resh_np = reshape(sta_gabor_np(t,:),[length(gab_xo)*length(gab_yo),length(lambdas),length(orientations)]);
    %     sta_gabor_pn(t,:) = mean(gabor_out_X_pn(tr_samp(spikebins),:));
    %     cur_gabor_resh_pn = reshape(sta_gabor_pn(t,:),[length(gab_xo)*length(gab_yo),length(lambdas),length(orientations)]);
    %     sta_gabor_pp(t,:) = mean(gabor_out_X_pp(tr_samp(spikebins),:));
    %     cur_gabor_resh_pp = reshape(sta_gabor_pp(t,:),[length(gab_xo)*length(gab_yo),length(lambdas),length(orientations)]);
    
    
    %     min_amp = min(sta_gabor(t,:));
    %     max_amp = max(sta_gabor(t,:));
    
    gsigma = Fsd*0.01;
    gauss_filt = fspecial('gaussian',round(3*gsigma),gsigma);
    
    %     for ll = 1:length(lambdas)
    %         if disp == 1;figure;end;
    %         clear cur_im
    %         cnt = 1;
    %         for i = 1:3
    %             for j = 1:3
    % %                 if cnt <= size(cur_gabor_resh,3)
    %                 if cnt <= size(cur_gabor_resh_nn,3)
    %
    % %                     cur_im(cnt,:,:) = reshape(cur_gabor_resh(:,ll,cnt),length(gab_yo),length(gab_xo));
    %                     cur_im_nn(cnt,:,:) = reshape(cur_gabor_resh_nn(:,ll,cnt),length(gab_yo),length(gab_xo));
    %                     cur_im_pn(cnt,:,:) = reshape(cur_gabor_resh_pn(:,ll,cnt),length(gab_yo),length(gab_xo));
    %                     cur_im_np(cnt,:,:) = reshape(cur_gabor_resh_np(:,ll,cnt),length(gab_yo),length(gab_xo));
    %                     cur_im_pp(cnt,:,:) = reshape(cur_gabor_resh_pp(:,ll,cnt),length(gab_yo),length(gab_xo));
    %
    %                     if disp == 1
    %                         subplot(3,3,cnt)
    %                         %             cur_im(cnt,:,:) = filter2(gauss_filt,squeeze(cur_im(cnt,:,:)),'same');
    %                         imagesc(gab_xo,gab_yo,squeeze(cur_im(cnt,:,:)));
    %                         hold on
    %                         rectangle('Position',[0.1 -0.7 0.6 0.6],'edgecolor','w','linewidth',2)
    %                         set(gca,'ydir','normal')
    %                         caxis([min_amp max_amp]);
    %                     end
    %                     cnt = cnt + 1;
    %                 end1
    %             end
    %         end
    %         %         [a,b] = nanmax(cur_im(:));
    %         %         [oriloc(t),ypeakloc(t),xpeakloc(t)] = ind2sub(size(cur_im),b);
    %         [a,b] = nanmax(cur_im_nn(:));
    %         [oriloc_nn(t),ypeakloc_nn(t),xpeakloc_nn(t)] = ind2sub(size(cur_im_nn),b);
    %         [a,b] = nanmax(cur_im_np(:));
    %         [oriloc_np(t),ypeakloc_np(t),xpeakloc_np(t)] = ind2sub(size(cur_im_np),b);
    %         [a,b] = nanmax(cur_im_pn(:));
    %         [oriloc_pn(t),ypeakloc_pn(t),xpeakloc_pn(t)] = ind2sub(size(cur_im_pn),b);
    %         [a,b] = nanmax(cur_im_pp(:));
    %         [oriloc_pp(t),ypeakloc_pp(t),xpeakloc_pp(t)] = ind2sub(size(cur_im_pp),b);
    %     end
    %
    %     if disp == 1
    %         subplot(3,3,oriloc(t))
    %         plot(gab_xo(xpeakloc(t)),gab_yo(ypeakloc(t)),'wo','linewidth',2)
    %         set(gcf,'Position',[500 1000 1200 1000])
    %     end90
    
    %     init_params(1) = gab_xo(xpeakloc(t)); %x0
    %     init_params(2) = gab_yo(ypeakloc(t)); %y0
    %     init_params(3) = orientations(oriloc(t)); %theta
    %     init_params(4) = 1/6; %lambda
    %     init_params(5) = 0.3*init_params(4); %sigma
    %     init_params(6) = 1; %eccentricity of gaussian ellipse
    %     init_params(7) = 0; %weight of quad term
    %     init_params(8) = 0; %weight of lin term 0 phase
    %     init_params(9) = 0; %weight of lin term pi/2 phase
    %     init_params(10) = 0; %const offset
    %     hold_const = [0 0 0 1 0 0 0 1 1 0];
    %     LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    %     UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];
    %     [gabor_params(t,:),LL(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    %     xv_LL(t) = get_pgabor_LL_v2(gabor_params(t,:),fullX(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    %
    %     init_params(1) = mean_x(t); %x0
    %     init_params(2) = mean_y(t); %y0
    %     init_params(3) = degtorad(pref_mu_ori(t)+90); %theta
    %     init_params(4) = 1/6; %lambda
    %     init_params(5) = 0.3*init_params(4); %sigma
    %     init_params(6) = 1; %eccentricity of gaussian ellipse
    %     init_params(7) = 0; %weight of quad term
    %     init_params(8) = 0; %weight of lin term 0 phase
    %     init_params(9) = 0; %weight of lin term pi/2 phase
    %     init_params(10) = 0; %const offset
    %     hold_const = [1 1 1 1 1 1 0 1 1 0];
    %     LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    %     UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];
    %     [gabor_params_emp(t,:),LL_emp(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank_emp(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_emp(t,1:6),0);
    %     xv_LL_emp(t) = get_pgabor_LL_v2(gabor_params_emp(t,:),fullX(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    
    %     init_params(1) = gab_xo(xpeakloc_nn(t)); %x0
    %     init_params(2) = gab_yo(ypeakloc_nn(t)); %,fullX_x(:,2100:3000))y0
    %     init_params(3) = orientations(oriloc_nn(t)); %theta
    %     init_params(4) = 1/6; %lambda
    %     init_params(5) = 0.3*init_params(4); %sigma
    %     init_params(6) = 1; %eccentricity of gaussian ellipse
    %     init_params(7) = 0; %weight of quad term
    %     init_params(8) = 0; %weight of lin term 0 phase
    %     init_params(9) = 0; %weight of lin term pi/2 phase
    %     init_params(10) = 0; %const offset
    %     hold_const = [0 0 0 1 0 0 0 1 1 0];
    %     LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    %     UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];
    %     [gabor_params_nn(t,:),LL_nn(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX_nn(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank_nn(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_nn(t,1:6),0);
    %     xv_LL_nn(t) = get_pgabor_LL_v2(gabor_params_nn(t,:),fullX_nn(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %weight of lin term 0 phase
    init_params(9) = 0; %weight of lin term pi/2 phase
    init_params(10) = 0; %const offset
    hold_const = [1 1 1 1 1 1 0 1 1 0];
    LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];    
    for i = 1:4
        init_params(3) = degtorad(pref_mu_ori(t)) + (i-1)*pi/4;
        [gabor_params_emp_nn(t,:),LL_emp_nn(t,i)] = fit_gabor_params_v3(XX,YY,init_params,fullX_nn(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
        gabor_bank_emp_nn(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_emp_nn(t,1:6),0);
        xv_LL_emp_nn(t,i) = get_pgabor_LL_v2(gabor_params_emp_nn(t,:),fullX_nn(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    end
    
    %     hold_const = [0 0 0 0 0 1 0 1 1 0];
    %     [gabor_params_uemp_nn(t,:),LL_uemp_nn(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX_nn(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank_uemp_nn(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_uemp_nn(t,1:6),0);
    %     xv_LL_uemp_nn(t) = get_pgabor_LL_v2(gabor_params_uemp_nn(t,:),fullX_nn(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    
    %     init_params(1) = gab_xo(xpeakloc_np(t)); %x0
    %     init_params(2) = gab_yo(ypeakloc_np(t)); %y0
    %     init_params(3) = orientations(oriloc_np(t)); %theta
    %     init_params(4) = 1/6; %lambda
    %     init_params(5) = 0.3*init_params(4); %sigma
    %     init_params(6) = 1; %eccentricity of gaussian ellipse
    %     init_params(7) = 0; %weight of quad term
    %     init_params(8) = 0; %weight of lin term 0 phase
    %     init_params(9) = 0; %weight of lin term pi/2 phase
    %     init_params(10) = 0; %const offset
    %     hold_const = [0 0 0 1 0 0 0 1 1 0];
    %     LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    %     UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];
    %     [gabor_params_np(t,:),LL_np(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX_np(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank_np(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_np(t,1:6),0);
    %     xv_LL_np(t) = get_pgabor_LL_v2(gabor_params_np(t,:),fullX_np(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %weight of lin term 0 phase
    init_params(9) = 0; %weight of lin term pi/2 phase
    init_params(10) = 0; %const offset
    hold_const = [1 1 1 1 1 1 0 1 1 0];
    LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];
    for i = 1:4
        init_params(3) = degtorad(pref_mu_ori(t)) + (i-1)*pi/4;
        [gabor_params_emp_np(t,:),LL_emp_np(t,i)] = fit_gabor_params_v3(XX,YY,init_params,fullX_np(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
        gabor_bank_emp_np(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_emp_np(t,1:6),0);
        xv_LL_emp_np(t,i) = get_pgabor_LL_v2(gabor_params_emp_np(t,:),fullX_np(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    end
    
    %     hold_const = [0 0 0 0 0 1 0 1 1 0];
    %     [gabor_params_uemp_np(t,:),LL_uemp_np(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX_np(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank_uemp_np(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_uemp_np(t,1:6),0);
    %     xv_LL_uemp_np(t) = get_pgabor_LL_v2(gabor_params_uemp_np(t,:),fullX_np(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    
    
    %         init_params(1) = gab_xo(xpeakloc_pn(t)); %x0
    %     init_params(2) = gab_yo(ypeakloc_pn(t)); %y0
    %     init_params(3) = orientations(oriloc_pn(t)); %theta
    %     init_params(4) = 1/6; %lambda
    %     init_params(5) = 0.3*init_params(4); %sigma
    %     init_params(6) = 1; %eccentricity of gaussian ellipse
    %     init_params(7) = 0; %weight of quad term
    %     init_params(8) = 0; %weight of lin term 0 phase
    %     init_params(9) = 0; %weight of lin term pi/2 phase
    %     init_params(10) = 0; %const offset
    %     hold_const = [0 0 0 1 0 0 0 1 1 0];
    %     LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    %     UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];
    %     [gabor_params_pn(t,:),LL_pn(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX_pn(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank_pn(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_pn(t,1:6),0);
    %     xv_LL_pn(t) = get_pgabor_LL_v2(gabor_params_pn(t,:),fullX_pn(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %weight of lin term 0 phase
    init_params(9) = 0; %weight of lin term pi/2 phase
    init_params(10) = 0; %const offset
    hold_const = [1 1 1 1 1 1 0 1 1 0];
    LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];
    for i = 1:4
        init_params(3) = degtorad(pref_mu_ori(t)) + (i-1)*pi/4;
        [gabor_params_emp_pn(t,:),LL_emp_pn(t,i)] = fit_gabor_params_v3(XX,YY,init_params,fullX_pn(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
        gabor_bank_emp_pn(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_emp_pn(t,1:6),0);
        xv_LL_emp_pn(t,i) = get_pgabor_LL_v2(gabor_params_emp_pn(t,:),fullX_pn(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    end

    
    %     hold_const = [0 0 0 0 0 1 0 1 1 0];
    %     [gabor_params_uemp_pn(t,:),LL_uemp_pn(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX_pn(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank_uemp_pn(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_uemp_pn(t,1:6),0);
    %     xv_LL_uemp_pn(t) = get_pgabor_LL_v2(gabor_params_uemp_pn(t,:),fullX_pn(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    
    %     init_params(1) = gab_xo(xpeakloc_pp(t)); %x0
    %     init_params(2) = gab_yo(ypeakloc_pp(t)); %y0
    %     init_params(3) = orientations(oriloc_pp(t)); %theta
    %     init_params(4) = 1/6; %lambda
    %     init_params(5) = 0.3*init_params(4); %sigma
    %     init_params(6) = 1; %eccentricity of gaussian ellipse
    %     init_params(7) = 0; %weight of quad term
    %     init_params(8) = 0; %weight of lin term 0 phase
    %     init_params(9) = 0; %weight of lin term pi/2 phase
    %     init_params(10) = 0; %const offset
    %     hold_const = [0 0 0 1 0 0 0 1 1 0];
    %     LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    %     UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];
    %     [gabor_params_pp(t,:),LL_pp(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX_pp(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank_pp(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_pp(t,1:6),0);
    %     xv_LL_pp(t) = get_pgabor_LL_v2(gabor_params_pp(t,:),fullX_pp(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    %
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %weight of lin term 0 phase
    init_params(9) = 0; %weight of lin term pi/2 phasetemp_np = bsxfun(@minus,xv_LL_emp_np,null_xvLL');
>> temp_nn = bsxfun(@minus,xv_LL_emp_nn,null_xvLL');
>> temp_pp = bsxfun(@minus,xv_LL_emp_pp,null_xvLL');
>> temp_pn = bsxfun(@minus,xv_LL_emp_pn,null_xvLL');
    init_params(10) = 0; %const offset
    hold_const = [1 1 1 1 1 1 0 1 1 0];
    LB = [0.1 -0.7 0 0.08 0.04 0.5 -Inf -Inf -Inf -Inf];
    UB = [0.7 -0.1 pi 0.5 0.5 2 Inf Inf Inf Inf];
    for i = 1:4
        init_params(3) = degtorad(pref_mu_ori(t)) + (i-1)*pi/4;
        [gabor_params_emp_pp(t,:),LL_emp_pp(t,i)] = fit_gabor_params_v3(XX,YY,init_params,fullX_pp(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
        gabor_bank_emp_pp(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_emp_pp(t,1:6),0);
        xv_LL_emp_pp(t,i) = get_pgabor_LL_v2(gabor_params_emp_pp(t,:),fullX_pp(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    end
    
    %     hold_const = [0 0 0 0 0 1 0 1 1 0];
    %     [gabor_params_uemp_pp(t,:),LL_uemp_pp(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX_pp(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    %     gabor_bank_uemp_pp(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_uemp_pp(t,1:6),0);
    %     xv_LL_uemp_pp(t) = get_pgabor_LL_v2(gabor_params_uemp_pp(t,:),fullX_pp(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
    
    avg_rate(t) = mean(all_binned_spks(tr_samp,t));
    null_xvLL(t) = -sum(bsxfun(@minus,all_binned_spks(xv_samp,t)*log(avg_rate(t)),avg_rate(t)))/sum(all_binned_spks(xv_samp,t));
    
    if disp == 1
        figure
        imagesc(xax(xpatch_inds),yax(ypatch_inds),squeeze(gabor_bank(t,:,:)));set(gca,'ydir','normal')
        rectangle('Position',[0.1 -0.7 0.6 0.6],'edgecolor','w','linewidth',2)
        hold on
        plot(gabor_params(t,1),gabor_params(t,2),'wo','linewidth',2)
        title(sprintf('LL: %.3f',LL(t)));
        figure
        subplot(3,1,1)
        plot(gr_oris,avg_mu_ori_profile(t,:))
        subplot(3,1,2)
        plot(parallel_x(t,:),parallel_profile(t,:))
        hold on
        plot(orth_x(t,:),orth_profile(t,:),'r')
        subplot(3,1,3)
        plot(parallel_y(t,:),parallel_profile(t,:))
        hold on
        plot(orth_y(t,:),orth_profile(t,:),'r')
        set(gcf,'Position',[0 1000 400 800])
        pause
        close all
    end
    
end

% pref_ori = orientations(oriloc);
% pref_x = gab_xo(xpeakloc);
% pref_y = gab_yo(ypeakloc);
%
%%
temp_np = bsxfun(@minus,xv_LL_emp_np,null_xvLL');
temp_nn = bsxfun(@minus,xv_LL_emp_nn,null_xvLL');
temp_pp = bsxfun(@minus,xv_LL_emp_pp,null_xvLL');
temp_pn = bsxfun(@minus,xv_LL_emp_pn,null_xvLL');

%%
% cd ~/Data/bruce/7_15_12/
% cd G029/
% % cd G035/
% load ./grating_mu_data
%
%
% close all
% clear beta
% oris = [gabor_props(:).ori];
% XO = [gabor_props(:).xo];
% YO = [gabor_props(:).yo];
% for t = 1:n_allunits;
%     cur_binned_spks = all_binned_spks(:,t);
%     for j = 1:length(orientations)
%         cur_set = find(oris == orientations(j));
%         beta(j,:) = corr(cur_binned_spks,gabor_out_X(:,cur_set));
%     end
%     max_amp = max(beta(:));
%     min_amp = min(beta(:));
%     cnt = 1;temp_np = bsxfun(@minus,xv_LL_emp_np,null_xvLL');
>> temp_nn = bsxfun(@minus,xv_LL_emp_nn,null_xvLL');
>> temp_pp = bsxfun(@minus,xv_LL_emp_pp,null_xvLL');
>> temp_pn = bsxfun(@minus,xv_LL_emp_pn,null_xvLL');
%     %     for i = 1:3
%     %         for j = 1:3
%     %             subplot(3,3,cnt)
%     %             imagesc(reshape(beta(cnt,:),length(gab_xo),length(gab_yo)))
%     %             set(gca,'ydir','normal')
%     %             caxis([min_amp max_amp])
%     %             cnt = cnt + 1;
%     %         end
%     %     end
%     plot(radtodeg(orientations),beta/max(beta))
%     hold on
%     % grat_est = avg_mu_ori_profile(single_units{1}(t),:);
%     grat_est = avg_mu_ori_profile(t,:);
%     grat_est = grat_est/max(grat_est);
%     % plot(unique_oris,grat_est,'r')
%     % pause
%     %     close all
%
%     [~,temp] = max(beta);
%     cur_pref_ori(t) = orientations(temp);
% end
%
% %%
SDIM=sdim;
% pref_oris = degtorad(pref_mu_ori(single_units{1}));
% pref_oris = degtorad([90 60 20
for t = single_units;
    fprintf('Fitting cell %d of %d\n',t,n_allunits);
    init_params(1) = 0; %x0
    init_params(2) = 0; %y0
    %     init_params(3) = cur_pref_ori(t); %theta
    init_params(3) = 0; %theta
    init_params(4) = 1/6*Fsd; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1.25; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %weight of lin term 0 phase
    init_params(9) = 0; %weight of lin term pi/2 phase
    init_params(10) = 0; %const offset
    hold_const = [0 0 0 0 0 0 0 1 1 0];
    LB = [-SDIM -SDIM 0 4 2 0.5 -Inf -Inf -Inf -Inf];
    UB = [SDIM SDIM pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf Inf];
    [gabor_params(t,:),LL(t)] = fit_gabor_params(init_params,fullX,all_binned_spks(:,t),[SDIM SDIM],hold_const,LB,UB);
    gabor_bank(t,:,:) = get_pgabor_mask(gabor_params(t,1:6),0,[SDIM SDIM]);
end
%
% lin_strength = sqrt(gabor_params(:,8).^2 + gabor_params(:,9).^2);
% quad_strength = gabor_params(:,7);
%
% %%
% cd ~/Data/bruce/7_15_12/
% cd G029/
% load ./grating_mu_data
% load ./onednoise_rffits
%
% sdim = 21;
% bandwidth = 0.25;
% ar = 1.25;
% kern_len = 21-1;
% [XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
% cur_mean_x = mean_x;
% cur_mean_y = mean_y;
% cur_pref_ori = degtorad(pref_mu_ori)+90;
% cur_pref_lambda = 1/4;
% for t = 1:n_allunits
%     rf_mapping_gabor(t,:,:) = get_gabor_template(XX,YY,cur_mean_x(t),...
%         cur_mean_y(t),cur_pref_ori(t),cur_pref_lambda,0,bandwidth,ar,Fsd);
% end
%
% %%
% for t = 1:15
%     subplot(2,1,1)
%     imagesc(squeeze(gabor_bank(t,:,:)));set(gca,'ydir','normal')
%     subplot(2,1,2)
%     imagesc(squeeze(rf_mapping_gabor(t,:,:)));set(gca,'ydir','normal')
%     pause
% end

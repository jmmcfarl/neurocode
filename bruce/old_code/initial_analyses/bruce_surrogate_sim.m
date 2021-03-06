
clear all
close all
addpath(genpath('~/Data/bruce/2_27_12'))
cd ~/Data/bruce/2_27_12/stimrecon/

load Blocks.mat

% original image resolution
% Pix2Deg = 1.1279 / 60;
Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;
Nyd = 256;
Nxd = 320;

% desired temporal resolution
stimres = 0.05; %in s
new_dt = 0.002;

% block number
blockid = 1;

stimtime = Blocks{blockid}.stimtime;
stimID = Blocks{blockid}.stimids;
Nimage = length(stimID);

% % load eye-movement
% load( ['~/Data/bruce/2_27_12/saccades/lemM232.' num2str(50+blockid) '.em.sac.mat']);
%
% % left eye
% lEye = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
% % right eye
% rEye = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
%
% EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
% EyeEndT = Expt.Trials.End/10000; % time of last eye sample
% Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
% eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
%
% recon_t = []; %time axis for reconstructed stim matrix (sec)
% stim_num = []; %corresponding vector of image indices
%
% gaussfilt = fspecial('gaussian',5,0.5);
%
% ov_stimRec = [];

%%
% x0 = 150;
% y0 = 150;
% sigma = [0.25 0.5 1]*Fsd; %in pixels
% % sigma = 1.5*Fsd; %in pixels
% theta = [0 pi/8 pi/4 pi/2];
% gamma = 1.2;
% lambda = [1 1.5 2]*Fsd;
% psi = 0;
% [X,Y] = meshgrid(1:Nxd,1:Nyd);
% cur_filt_num = 1;
% for i = 1:length(lambda)
%     for j = 1:length(theta)
%         for k = 1:length(sigma)
%             Xprime = (X-x0)*cos(theta(j))+(Y-y0)*sin(theta(j));
%             Yprime = -(X-x0)*sin(theta(j))+(Y-y0)*cos(theta(j));
%             filt{cur_filt_num} = exp(-(Xprime.^2+gamma^2*Yprime.^2)/(2*sigma(k)^2)) .* cos(2*pi*Xprime/lambda(i)+psi);
%             filt_r(:,cur_filt_num) = filt{cur_filt_num}(:)/1.7046e3;
%             
%             filt_props(cur_filt_num).theta = theta(j);
%             filt_props(cur_filt_num).gamma = gamma;
%             filt_props(cur_filt_num).lambda = lambda(i);
%             filt_props(cur_filt_num).psi = psi;
%             filt_props(cur_filt_num).x0 = x0;
%             filt_props(cur_filt_num).y0 = y0;
%             filt_props(cur_filt_num).sigma = sigma(k);
%             
%             cur_filt_num = cur_filt_num + 1;
%         end
%     end
% end
% figure
% imagesc((1:Nyd)/Fsd,(1:Nxd)/Fsd,filt);

[X,Y] = meshgrid(1:Nxd,1:Nyd);

%[x0 y0 sigma theta lambda gamma psi]
% filt_props(1,:) = [150 150 0.5*Fsd pi/4 0.5*Fsd 1.2 0];
% filt_props(2,:) = [150 150 0.5*Fsd pi/4 1*Fsd 1.2 0];
% filt_props(3,:) = [150 150 0.5*Fsd pi/4 1.5*Fsd 1.2 0];
% filt_props(4,:) = [150 150 0.5*Fsd pi/4 2*Fsd 1.2 0];
% filt_props(5,:) = [150 150 0.5*Fsd pi/2 1*Fsd 1.2 0];
% filt_props(6,:) = [150 150 0.5*Fsd pi/8 1*Fsd 1.2 0];
filt_props(1,:) = [150 150 0.5*Fsd pi/2 0.5*Fsd 1.2 0];
filt_props(2,:) = [150 150 0.5*Fsd pi/2 1*Fsd 1.2 0];
filt_props(3,:) = [150 150 0.5*Fsd pi/2 1.5*Fsd 1.2 0];
filt_props(4,:) = [150 150 0.5*Fsd 0 1*Fsd 1.2 0];
filt_props(5,:) = [150 150 0.5*Fsd pi/4 1*Fsd 1.2 0];
for i = 1:size(filt_props,1)
    Xprime = (X-filt_props(i,1))*cos(filt_props(i,4))+(Y-filt_props(i,2))*sin(filt_props(i,4));
    Yprime = -(X-filt_props(i,1))*sin(filt_props(i,4))+(Y-filt_props(i,2))*cos(filt_props(i,4));
    filt{i} = exp(-(Xprime.^2+filt_props(i,6)^2*Yprime.^2)/(2*filt_props(i,3)^2)) .* ...
        cos(2*pi*Xprime/filt_props(i,5) + filt_props(i,7));
    filt_r(:,i) = filt{i}(:)/1.7046e3;
end



%%
beta = 1;
offset = 0;
scale = 2;
cd ~/Data/bruce/2_27_12/stimrecon/
recon_t = [];
gen_funs = [];
spike_probs = [];
new_recon_t = [];
% all_spike_bins = [];

% ov_avg_stim = zeros(Nyd,Nxd);
% ov_spktrg_avg = zeros(Nyd,Nxd);
% n_tot_spks = 0;
n_recon_samps = zeros(Nimage,1);
stim_num = [];

delta_bin = 0;

% load images
for i=1:Nimage
    % for i = 1:5
    i
    load(['BLOCK' num2str(blockid) 'IMAGE' num2str(i) 'r2_p5.mat'], 'STIMrec');
    n_recon_samps(i) = size(STIMrec,1);
    
    STIMmat = reshape(STIMrec,[n_recon_samps(i) Nxd*Nyd]);
    
    cur_gen_funs = STIMmat*filt_r;
    
    cur_gen_funs = cur_gen_funs.^2; %squaring NL
    
    cur_spike_prob = scale/beta*log(1+exp(beta*(cur_gen_funs-offset)));
    
    gen_funs = [gen_funs; cur_gen_funs];
    spike_probs = [spike_probs; cur_spike_prob];
    
    %         cur_spikes_binned = poissrnd(cur_spike_prob);
    %         un_spk_cnts = unique(cur_spikes_binned);
    %         spk_bins = [];
    %         for tt = 2:length(un_spk_cnts)
    %             spk_bins = [spk_bins; repmat(find(cur_spikes_binned==un_spk_cnts(tt)),un_spk_cnts(tt),1)];
    %         end
    %         spk_bins = spk_bins + delta_bin;
    %
    %         all_spike_bins = [all_spike_bins; cur_spikes_binned];
    %
    %         cur_spkavg_stim = squeeze(nansum(STIMrec(spk_bins,:,:),1));
    %         ov_spktrg_avg = squeeze(ov_spktrg_avg) + cur_spkavg_stim;
    %
    %         n_tot_spks = n_tot_spks + length(spk_bins);
    %         cur_avg_stim = squeeze(nanmean(STIMrec));
    %
    %         ov_avg_stim = ov_avg_stim + cur_avg_stim*n_recon_samps(i);
    
    
    stim_num = [stim_num i*ones(1,n_recon_samps(i))];
    
    recon_t = [recon_t stimtime(i):stimres:(stimtime(i)+(n_recon_samps(i)-1)*stimres)];
    new_recon_t = [new_recon_t stimtime(i):new_dt:recon_t(end)];
end
% ov_avg_stim = ov_avg_stim/sum(n_used_samps);
% ov_spktrg_avg = bsxfun(@rdivide,ov_spktrg_avg,n_tot_spks);
% ov_spktrg_avg_cor = ov_spktrg_avg - ov_avg_stim;


%%
norm_gen_funs = zscore(gen_funs);
spike_prob_post = scale/beta*log(1+exp(beta*(norm_gen_funs-offset)));

% interp_axis = recon_t(1):new_dt:recon_t(end);
spike_prob_int = interp1(recon_t',spike_prob_post,new_recon_t)*new_dt/stimres;

spikes = poissrnd(spike_prob_int);
for j = 1:size(filt_r,2)
    un_spk_cnts = unique(spikes(:,j));
    spikebins = [];
    for i = 1:length(un_spk_cnts)
        cur_set = find(spikes(:,j) == un_spk_cnts(i));
        spikebins = [spikebins; repmat(cur_set(:),un_spk_cnts(i),1)];
    end
    % spikebins = unique(spikebins); %max 1 spike per bin
    spiketimes{j} = sort(new_recon_t(spikebins));
    fprintf('Cell %d Nspks: %d\n',j,length(spiketimes{j}));
end
%%
% save surrogate_spikes_CC2_40hz spiketimes all_spike_bins *recon_t stim_num all_spike_bins beta offset scale theta sigma gamma lambda x0 y0
save surrogate_spikes_CCset2 spiketimes *recon_t stim_num beta offset scale filt_props


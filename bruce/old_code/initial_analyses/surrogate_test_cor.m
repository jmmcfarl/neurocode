
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

% stimtime = Blocks{blockid}.stimtime;
% stimID = Blocks{blockid}.stimids;
% Nimage = length(stimID);

%%
Nyp = 1024;
Nxp = 1280;
Ny = Nyp/dsfrac;
Nx = Nxp/dsfrac;
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;

RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

x0 = find(xax > 4.5,1,'first');
y0 = find(yax > -2.5,1,'first');

xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));
SDIM = length(xpatch_inds);
[X,Y] = meshgrid(1:Nxd,1:Nyd);

%[x0 y0 sigma theta lambda gamma psi]
% filt_props(1,:) = [150 150 0.5*Fsd pi/4 0.5*Fsd 1.2 0];
% filt_props(2,:) = [150 150 0.5*Fsd pi/4 1*Fsd 1.2 0];
% filt_props(3,:) = [150 150 0.5*Fsd pi/4 1.5*Fsd 1.2 0];
% filt_props(4,:) = [150 150 0.5*Fsd pi/4 2*Fsd 1.2 0];
% filt_props(5,:) = [150 150 0.5*Fsd pi/2 1*Fsd 1.2 0];
% filt_props(6,:) = [150 150 0.5*Fsd pi/8 1*Fsd 1.2 0];
filt_props(1,:) = [x0 y0 0.3*Fsd pi/4 0.5*Fsd 1.2 0];
filt_props(2,:) = [x0 y0 0.3*Fsd pi/4 1*Fsd 1.2 0];
filt_props(3,:) = [x0 y0 0.3*Fsd pi/4 1.5*Fsd 1.2 0];
% filt_props(4,:) = [150 150 0.5*Fsd 0 1*Fsd 1.2 0];
% filt_props(5,:) = [150 150 0.5*Fsd pi/4 1*Fsd 1.2 0];
for i = 1:size(filt_props,1)
    Xprime = (X-filt_props(i,1))*cos(filt_props(i,4))+(Y-filt_props(i,2))*sin(filt_props(i,4));
    Yprime = -(X-filt_props(i,1))*sin(filt_props(i,4))+(Y-filt_props(i,2))*cos(filt_props(i,4));
    filt{i} = exp(-(Xprime.^2+filt_props(i,6)^2*Yprime.^2)/(2*filt_props(i,3)^2)) .* ...
        cos(2*pi*Xprime/filt_props(i,5) + filt_props(i,7));
    filt_patch{i} = filt{i}(ypatch_inds,xpatch_inds);
    filt_patch{i} = filt_patch{i} - mean(filt_patch{i}(:));
    filt_r(:,i) = filt_patch{i}(:)/1.7046e3;
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
% n_recon_samps = zeros(Nimage,1);
% stim_num = [];

delta_bin = -1;

for blockid = 1:4
    fprintf('Block %d of %d\n',blockid,4);
    sname = sprintf('image_patch_block%d_leye_50_dsf4',blockid);
    load(sname);
    
    gen_funs = reshape(ov_im_patch,length(recon_t),SDIM^2)*filt_r;
    gen_funs(isnan(gen_funs)) = 0;
    gen_funs = [zeros(-delta_bin,size(gen_funs,2)); gen_funs(1:end+delta_bin,:)];
    
    
    % create vector determining whether current stimulus was high-pass filtered
    filt_stims = mod(0:500,4) + 1 >2;
    is_stim_filtered = filt_stims(stim_num);
    
    %%
    % norm_gen_funs = zscore(gen_funs);
    norm_gen_funs = zscore(gen_funs);
    spike_prob_post = scale/beta*log(1+exp(beta*(norm_gen_funs-offset)));
    
    spikes = poissrnd(spike_prob_post);
    for j = 1:size(filt_r,2)
        un_spk_cnts = unique(spikes(:,j));
        spikebins = [];
        for i = 1:length(un_spk_cnts)
            cur_set = find(spikes(:,j) == un_spk_cnts(i));
            spikebins = [spikebins; repmat(cur_set(:),un_spk_cnts(i),1)];
        end
        % spikebins = unique(spikebins); %max 1 spike per bin
        spiketimes{blockid,j} = sort(recon_t(spikebins));
        fprintf('Cell %d Nspks: %d\n',j,length(spiketimes{blockid,j}));
    end
    
end
%%
% save surrogate_spikes_CC2_40hz spiketimes all_spike_bins *recon_t stim_num all_spike_bins beta offset scale theta sigma gamma lambda x0 y0
save surrogate_spikes_CCset_cor spiketimes filt_props


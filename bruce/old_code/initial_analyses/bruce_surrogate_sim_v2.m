
clear all
% close all
addpath(genpath('~/Data/bruce/2_27_12'))
cd ~/Data/bruce/2_27_12/stimrecon/

load Blocks.mat

% original image resolution
Pix2Deg = 1.1279 / 60;

% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;
Nyd = 256;
Nxd = 320;

% desired temporal resolution
stimres = 0.025; %in s
new_dt = 0.002;

%%
x0 = 150;
y0 = 150;
sigma = 0.45*Fsd; %in pixels
theta =pi/4;
gamma = 1.3;
lambda = 1.6*Fsd;
psi = 0;
[X,Y] = meshgrid(1:Nxd,1:Nyd);
Xprime = (X-x0)*cos(theta)+(Y-y0)*sin(theta);
Yprime = -(X-x0)*sin(theta)+(Y-y0)*cos(theta);
filt = exp(-(Xprime.^2+gamma^2*Yprime.^2)/(2*sigma^2)) .* cos(2*pi*Xprime/lambda+psi);
filt_r = filt(:);
filt_r = filt_r/1.7046e3; %this kludge normalizes gen fun to unit variance for 'on-the-fly- decoding'

beta = 1;
theta = 2;
scale = 2;

delta_bin = 0;

%%
% block number
cd ~/Data/bruce/2_27_12/stimrecon/

ov_avg_stim = zeros(Nyd,Nxd);
ov_spktrg_avg = zeros(Nyd,Nxd);
n_tot_spks = 0;
n_used_samps = 0;

for blockid = 1:4
    
    stimtime = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    Nimage = length(stimID);
    
    recon_t{blockid} = [];
    gen_fun{blockid} = [];
    spike_prob{blockid} = [];
    new_recon_t{blockid} = [];
    stim_num{blockid} = [];
    
    % load images
    for i=1:Nimage
        i
        load(['BLOCK' num2str(blockid) 'IMAGE' num2str(i) 'r2.mat'], 'STIMrec');
        n_recon_samps = size(STIMrec,1);
        STIMmat = reshape(STIMrec,[n_recon_samps Nxd*Nyd]);
        cur_gen_fun = STIMmat*filt_r;
        gen_fun{blockid} = [gen_fun{blockid}; cur_gen_fun];
        cur_spike_prob = scale/beta*log(1+exp(beta*cur_gen_fun-theta));
        spike_prob{blockid} = [spike_prob{blockid}; cur_spike_prob];
        
        cur_spikes_binned = poissrnd(cur_spike_prob);
        un_spk_cnts = unique(cur_spikes_binned);
        spk_bins = [];
        for tt = 2:length(un_spk_cnts)
            spk_bins = [spk_bins; repmat(find(cur_spikes_binned==un_spk_cnts(tt)),un_spk_cnts(tt),1)];
        end
        spk_bins = spk_bins + delta_bin;
        
        %         all_spike_bins{i} = spk_bins;
        
        cur_spkavg_stim = squeeze(nansum(STIMrec(spk_bins,:,:),1));
        ov_spktrg_avg = squeeze(ov_spktrg_avg) + cur_spkavg_stim;
        
        n_tot_spks = n_tot_spks + length(spk_bins);
        cur_avg_stim = squeeze(nanmean(STIMrec));
        
        cur_nsamps = size(STIMrec,1);
        n_used_samps = n_used_samps + cur_nsamps;
        ov_avg_stim = ov_avg_stim + cur_avg_stim*cur_nsamps;
        
        stim_num{blockid} = [stim_num{blockid} i*ones(1,n_recon_samps)];
        
        recon_t{blockid} = [recon_t{blockid} stimtime(i):stimres:(stimtime(i)+(n_recon_samps-1)*stimres)];
        new_recon_t{blockid} = [new_recon_t{blockid} stimtime(i):new_dt:recon_t{blockid}(end)];
    end
    
end

ov_avg_stim = ov_avg_stim/sum(n_used_samps);
ov_spktrg_avg = bsxfun(@rdivide,ov_spktrg_avg,n_tot_spks);
ov_spktrg_avg_cor = ov_spktrg_avg - ov_avg_stim;

%%
for blockid = 1:4
norm_gen_fun = gen_fun{blockid};
% norm_gen_fun = zscore(gen_fun);
spike_prob_post = scale/beta*log(1+exp(beta*(norm_gen_fun-theta)));

% interp_axis = recon_t(1):new_dt:recon_t(end);
spike_prob_int = interp1(recon_t{blockid}',spike_prob_post,new_recon_t{blockid})*new_dt/stimres;

spikes = poissrnd(spike_prob_int);
un_spk_cnts = unique(spikes);
spikebins = [];
for i = 1:length(un_spk_cnts)
    cur_set = find(spikes == un_spk_cnts(i));
    spikebins = [spikebins; repmat(cur_set(:),un_spk_cnts(i),1)];
end
% spikebins = unique(spikebins); %max 1 spike per bin
spiketimes{blockid} = sort(new_recon_t{blockid}(spikebins));
fprintf('Nspks: %d\n',length(spikebins));

end
%%
% save surrogate_spikes_fourblock spiketimes all_spike_bins *recon_t stim_num
save surrogate_spikes_fourblock spiketimes



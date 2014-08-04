clear all
% close all
addpath(genpath('~/James_scripts'));

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
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

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3


SDIM = length(xpatch_inds);
dt = .005;

cellids = [1:10];
muaids = [1:14];

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
backwin = 0.3;
forwardwin = 1;
all_model_time_axis = [];
all_model_blockids = [];
all_model_fixids = [];
spikes_binned = [];
time_since_stim = [];
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    cur_stim_times = Blocks{blockid}.stimtime;
    
    for i = 1:length(cur_stim_times)
        start_T = cur_stim_times(i)-backwin;
        end_T = start_T + forwardwin;
        cur_tedges = start_T:dt:end_T;
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        all_model_time_axis = [all_model_time_axis cur_tcents];
        all_model_blockids = [all_model_blockids blockid*ones(1,length(cur_tcents))];
        temp_binned = zeros(10,length(cur_tcents));
        temp_modouts = zeros(10,length(cur_tcents));
        for c = 1:n_used_cells
            if c <= 10
                temp = histc(Blocks{blockid}.spktimes{cellids(c)},cur_tedges);
            else
                temp = histc(Blocks{blockid}.mutimes{muaids(c-10)},cur_tedges);
            end
            temp_binned(c,:) = temp(1:end-1);
        end
        spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_stim = [time_since_stim cur_tcents-start_T-backwin];
        
    end
end

%% Compute TBR time-since fix onset
max_tsf = 1; nbins = 50;
used_tsf = time_since_stim;
used_tsf(used_tsf > max_tsf) = [];
tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));

Tmat = tbrep(time_since_stim,tax);


%% COMPUTE SACCADE TRIG AVG FIRING RATES
zspikes = bsxfun(@rdivide,spikes_binned,mean(spikes_binned));
nrate = mean(zspikes,2);

stim_trg_rate = zeros(n_used_cells,length(tax));
stim_trg_nrate = zeros(1,length(tax));
stim_trg_var = zeros(n_used_cells,length(tax));
taxb = [tax(1)-0.05 tax];
for i = 1:length(tax)
    curset = find(time_since_stim >= taxb(i) & time_since_stim < taxb(i+1));
    stim_trg_nrate(i) = mean(nrate(curset));
    for c = 1:n_used_cells
        stim_trg_rate(c,i) = mean(spikes_binned(curset,c));
%         stim_trg_rate(c,i) = mean(zspikes(curset,c));
        stim_trg_var(c,i) = var(spikes_binned(curset,c));
    end
    n_occ(i) = length(curset);
end
avg_rates = mean(spikes_binned)/dt;
fano_fac = stim_trg_var./stim_trg_rate;
stim_trg_rate = stim_trg_rate/dt;
stim_trg_nrate = stim_trg_nrate;
stim_trg_std = sqrt(stim_trg_var)/dt;
norm_str = bsxfun(@rdivide,stim_trg_rate,avg_rates');
% norm_str = bsxfun(@rdivide,stim_trg_rate,max(stim_trg_rate,[],2));

%%

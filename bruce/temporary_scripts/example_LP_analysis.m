clear all
close all

cd ~/Data/bruce/M232
addpath('~/Matlab/James/code_iSTAC/')
load ./lemM232Expts.mat
%%
bar_expts = [37 38 39 43 46 47];
rf_cent = [4.73 -2.2]; %deg
axis_or = 50*pi/180; %in radians

Fs = 1000; %reported sample rate of Bruce-processed LFP (1kHz on the money)
niqf = Fs/2;

[bb,aa] = butter(2,[1 100]/niqf); %band-pass filtering of LFP

%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 35;
scales = logspace(log10(10),log10(500),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fs);

%temporal parameters of model
dt = 0.01; %stimulus dt
usfrac = 2;
new_dt = dt/usfrac;
tbspace = 2;%tent-basis spacing
nLags = 15; %number of tent bases

n_bar_pos = 15;

stim_params = NIMcreate_stim_params([nLags n_bar_pos],new_dt,1,tbspace);

%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
use_units = 1:24;
%%
load ./un_bar_pos.mat
un_bar_pos(1:2) = [];

full_Xmat = [];
full_spkbinned = [];
full_expt_vec = [];
full_trial_vec = [];
full_taxis = [];
full_bar_pos = [];
full_t_since_tstart = [];
for ee = 1:length(bar_expts);
    fprintf('Processing expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
%     trial_durs{ee} = [Expts{bar_expts(ee)}.Trials(:).dur];
    n_trials(ee) = length(Expts{bar_expts(ee)}.Trials);
    all_t_axis = [];
    %all_used_inds = [];
    all_expt_vec = [];
    all_trial_vec = [];
    all_bar_Op = [];
    all_bar_Xmat = [];
    all_binned_spikes = [];
    
    all_t_since_trial_start = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
        
        cur_bar_Op_new = repmat(cur_bar_Op,1,usfrac)';
        cur_bar_Op_new = cur_bar_Op_new(:);
        
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_t_edges_new = Expts{bar_expts(ee)}.Trials(tt).Start/1e4:new_dt:Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4;
        cur_t_axis_new = 0.5*cur_t_edges_new(1:end-1) + 0.5*cur_t_edges_new(2:end);
        
        cur_bar_Op_new(length(cur_t_axis_new)+1:end) = [];
        
        %bin spikes for this trial
        cur_binned_spks = nan(length(cur_t_axis_new),length(use_units));
        for cc = 1:length(use_units)
            cur_hist = histc(Clusters{use_units(cc)}.times,cur_t_edges_new);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
                
        %compute matrix of bar-positions (1 if the bar is at that position,
        %0 otherwise)
        cur_bar_mat = zeros(length(cur_t_axis_new),length(un_bar_pos));
        for b = 1:length(un_bar_pos)
            cur_set = find(cur_bar_Op_new==un_bar_pos(b));
            cur_bar_mat(cur_set,b) = 1;
        end
        bar_Xmat = create_time_embedding(cur_bar_mat,stim_params);
                
        cur_t_since_tstart = cur_t_axis_new - cur_t_axis(1);
                
        all_t_axis = [all_t_axis; cur_t_axis_new(:)];
        all_t_since_trial_start = [all_t_since_trial_start; cur_t_since_tstart(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis_new),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis_new),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op_new(:)];
        all_bar_Xmat = [all_bar_Xmat; bar_Xmat];
    end
        
    full_Xmat = [full_Xmat; all_bar_Xmat];
    full_spkbinned = [full_spkbinned; all_binned_spikes];
    full_expt_vec = [full_expt_vec; ones(length(all_t_axis),1)*ee];
    full_taxis = [full_taxis; all_t_axis];
    full_t_since_tstart = [full_t_since_tstart; all_t_since_trial_start];
    full_trial_vec = [full_trial_vec; all_trial_vec];
    full_bar_pos = [full_bar_pos; all_bar_Op];
        
end

%%
full_lfps = [];
full_phasegrams = [];
full_ampgrams = [];
n_lfp_trials = nan(length(bar_expts),1);
for ee = 1:length(bar_expts);
    fprintf('Expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('lemM232A.%d.lfp.mat',bar_expts(ee));
    load(fname);
    
%     Fs = 1/LFP.Header.CRsamplerate;
    
%compile LFP data for current expt
    n_lfp_trials(ee) = length(LFP.Trials);
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    expt_lfp_t_axis = [];
    expt_lfps = [];
    for tt = 1:n_lfp_trials(ee)
        
        cur_npts = size(LFP.Trials(tt).LFP,1);
        cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
        cur_t_axis = lfp_trial_starts(tt):1/Fs:cur_t_end(tt);
        
        %handle overlap of LFP data
        if ~isempty(expt_lfp_t_axis)
            cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
        else
            cur_sp = 1;
        end
        cur_t_axis = cur_t_axis(cur_sp:end);
        
        cur_LFP = [LFP.Trials(tt).LFP];
        cur_LFP = cur_LFP(cur_sp:end,:);
        
        cur_LFP = filtfilt(bb,aa,cur_LFP);               
        
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
        expt_lfps = [expt_lfps; cur_LFP];
    end
    
%     %compute CSD if you want
%     vars.Fs = Fs;
%     vars.BrainBound = 1;
%     vars.ChanSep = 0.05;
%      vars.diam = 2; %0.5
%     CSD = PettersenCSD(expt_lfps','spline',vars)';
%     expt_lfps = CSD; %switch LFP to CSD

    cur_set = find(full_expt_vec==ee); %stimulus samples for the current expt block
    
    %interpolate LFP samples onto stim t-axis
    interp_lfps = interp1(expt_lfp_t_axis,expt_lfps,full_taxis(cur_set));
    
    %wavelet analysis
    cur_phasegram = nan(length(expt_lfp_t_axis),length(wfreqs),24);
    cur_ampgram = nan(length(expt_lfp_t_axis),length(wfreqs),24);
    for ll = 1:24
        fprintf('Wavelet transform on channel %d of %d\n',ll,24);
        temp = cwt(expt_lfps(:,ll),scales,'cmor1-1');
        cur_phasegram(:,:,ll) = angle(temp)';
        cur_ampgram(:,:,ll) = abs(temp)';
    end
    
    %interpolate wavelet transform
    unwr_phasegram = unwrap(cur_phasegram); %unwrap phase data before linear interpolation
    interp_phasegrams = interp1(expt_lfp_t_axis,unwr_phasegram,full_taxis(cur_set));
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi; %re-wrap
    interp_ampgrams = interp1(expt_lfp_t_axis,cur_ampgram,full_taxis(cur_set));
    
    full_lfps = [full_lfps; interp_lfps];
    full_phasegrams = cat(1,full_phasegrams, interp_phasegrams);
    full_ampgrams = cat(1,full_ampgrams, interp_ampgrams);
end

% normalize amplitude spectra to have unit variance
full_ampgrams = bsxfun(@rdivide,full_ampgrams,std(full_ampgrams(used_inds,:,:)));


%% Parse into training and XV sets
[c,ia,ic] = unique([full_expt_vec full_trial_vec],'rows'); %assigns unique trial indices
n_trials = length(ia);

rp = randperm(n_trials);
rand_trial_vec = rp(ic);
[~,ind_shuff] = sort(rand_trial_vec);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);

%create set of cross-validation trials
xv_tset = randperm(n_trials);
xv_tset(n_xv_trials+1:end) = [];

tr_tset = find(~ismember(1:n_trials,xv_tset));%training trials

xv_inds = find(ismember(ic,xv_tset));
tr_inds = find(ismember(ic,tr_tset))';

%eliminate data near the beginnings and ends of trials
buffer_time = 0.15; %in sec
trial_dur = 2;
used_inds = true(length(ic),1);
used_inds(full_t_since_tstart < buffer_time | full_t_since_tstart > (trial_dur-buffer_time)) = false;
xv_inds(~used_inds(xv_inds)) = [];
tr_inds(~used_inds(tr_inds)) = [];


%% Fit STRFs
reg_params = NIMcreate_reg_params('lambda_d2XT',200);
for cc = 1:24
    
    Robs = full_spkbinned(tr_inds,cc);
    
    mod_signs = [1]; %determines whether input is exc or sup (doesn't matter in the linear case)
    NL_types = {'lin'}; %define subunit as linear
    silent = 1; %display model optimization iterations
    
    fit0(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
    fit0(cc) = NIMfit_filters(fit0(cc),Robs,full_Xmat(tr_inds,:),[],[],silent); %fit stimulus filters
    
end

%%

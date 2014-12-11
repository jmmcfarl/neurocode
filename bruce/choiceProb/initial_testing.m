%% LOAD PROCESSED DATA
cd('/home/james/Data/bruce/ChoiceProb/')

load('M239/lemM239.image.ORBW.LFP.mat')

%% Assemble stim, Robs, LFP into trial structure
Ntrials = length(AllExpt.Expt.Trials);
LFP_Fs = 1.0 / AllExpt.Expt.Header.LFPsamplerate;  % this is 1 kHz 
dt = 0.01;  % assuming 10 ms frame time exactly
trialDur = 2; %in sec.
NT = 200; %number of video frames (@100Hz) per trial

% SU and MU data
Nunits = length(AllExpt.Header);
SUindx = []; MUindx = [];  SUchans = []; MUchans = [];  SUcellnumbers = [];
for nn = 1:Nunits
  if AllExpt.Header(nn).cellnumber > 0
    SUcellnumbers(end+1) = AllExpt.Header(nn).cellnumber;
    SUindx(end+1) = nn;
    SUchans(end+1) = AllExpt.Header(nn).probe;
  else
    MUindx(end+1) = nn;
    MUchans(end+1) = AllExpt.Header(nn).probe;
  end    
end

%%
trialNums = [AllExpt.Expt.Trials(:).Trial];
trialOR = [AllExpt.Expt.Trials(:).or];
trialOB = [AllExpt.Expt.Trials(:).ob];
trialrwDir = [AllExpt.Expt.Trials(:).rwdir];
trialrespDir = [AllExpt.Expt.Trials(:).RespDir];
trialSe = [AllExpt.Expt.Trials(:).se];

%% Get Binned Spikes
fullRobs = nan(NT,Ntrials,length(SUindx));
for tr = 1:Ntrials
  for nn = 1:length(SUindx)
    trindx = find( AllExpt.Spikes{SUindx(nn)}.Trial == trialNums(tr));
    if ~isempty(trindx)
      spks = double( AllExpt.Spikes{SUindx(nn)}.Spikes{trindx} ) * 1e-4;  % no adjustment from start of trial
      cur_Robs = histc( spks, (0:(NT))*dt);
      fullRobs(:,tr,nn) = cur_Robs(1:end-1);
    end
  end 
end
totSpkCnts = squeeze(nansum(fullRobs));
%%
LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes*1e-4;
R_trial_taxis = (0:(NT-1))*dt + dt/2;

nprobes = 24;

trial_LFP_state = nan(NT,Ntrials,nprobes);

lcf = 2; hcf = 4; %cut-off freqs for LFP filter
[filt_bb,filt_aa] = butter(2,[lcf hcf]/(LFP_Fs/2));
for tr = 1:Ntrials
    cur_LFPs = double(AllExpt.Expt.Trials(tr).LFP);
    cur_LFPs(isnan(cur_LFPs)) = 0;
    cur_LFPs = filtfilt(filt_bb,filt_aa,cur_LFPs);
    cur_phase = unwrap(angle(hilbert(cur_LFPs)));
    cur_phase_interp = mod(interp1(LFP_trial_taxis,cur_phase,R_trial_taxis),2*pi);
    trial_LFP_state(:,tr,:) = cur_phase_interp;
end

%%
clear all; 
%%%%%%% INITIAL PROCESSING OF Tim Data
dt   = 19.9920031987/1000; 
sdim = 64; 

rstimfiles = {'af2-64p-50h-2m-mn125-ct015.s0', 'wn2-64p-50h-2m-mn125-ct015.s1', 'pt2-64p-50h-2m-mn125-ct015.s2',...
    'ps2-64p-50h-2m-mn125-ct015.s3', 'ns2-64p-50h-2m-mn125-ct015.s4', 'pn2-64p-50h-2m-mn125-ct015.s5',...
    'af2-64p-50h-2m-mn125-ct025.s6', 'wn2-64p-50h-2m-mn125-ct025.s7', 'ns2-64p-50h-2m-mn125-ct025.s8',...
    'pn2-64p-50h-2m-mn125-ct025.s9', 'ps2-64p-50h-2m-mn125-ct025.s10', 'pt2-64p-50h-2m-mn125-ct025.s11',...
    'pn2-64p-50h-2m-mn125-ct035.s12', 'ps2-64p-50h-2m-mn125-ct035.s13', 'ns2-64p-50h-2m-mn125-ct035.s14',...
    'pt2-64p-50h-2m-mn125-ct035.s15', 'af2-64p-50h-2m-mn125-ct035.s16', 'wn2-64p-50h-2m-mn125-ct035.s17',...
    'pt2-64p-50h-2m-mn125-ct045.s18', 'pn2-64p-50h-2m-mn125-ct045.s19', 'ns2-64p-50h-2m-mn125-ct045.s20',...
    'af2-64p-50h-2m-mn125-ct045.s21', 'ps2-64p-50h-2m-mn125-ct045.s22', 'wn2-64p-50h-2m-mn125-ct045.s23'};

stimfiles  = cellfun(@(x)x(1:26),rstimfiles,'UniformOutput',0); 
nstims = length(stimfiles);

datdir = '~/Data/blanche/rec_75';
%stimfiles = dir([datdir,'/stim_data/*.stim']); 
cd(datdir)
save stimfiles75.mat stimfiles;  

nfiles=length(stimfiles); allstims =[]; rstimparts= cell(nfiles,1); 
for istim =1:nfiles; 
	tstimfile = stimfiles{istim}; %current dataset
    %read in current raw data
	fid    = fopen([datdir,'/stim_data/',tstimfile],'r');
	  rawdat = fread(fid,inf,'uchar'); 
	fclose(fid);
    
    %compute stats on data
	[min(rawdat),max(rawdat),mean(rawdat),var(rawdat)]
	%rawdat = (rawdat-mean(rawdat))/sqrt(var(rawdat)); 
	rawdat = 2*((rawdat/255)-0.5); %not sure why this transformation 
	temp(istim) = length(rawdat);
    rstimparts{istim} = rawdat; 
	allstims=[allstims;rawdat]; 
end; 


%% reshape, normalize & downsample the stimulus [-1,1]
dsr=2; %spatial downsampling factor
%stims75   = reshape(allstims, sdim*sdim,nfiles*6000)';  
%dstims75  = dsamplebstim(stims75,dsr); 
stimps75  = cellfun(@(x)reshape(x,sdim*sdim,6000)',rstimparts,'UniformOutput',0); %reshape the two minutes of data 
dstimps75 = cellfun(@(x)dsamplebstim(x,dsr),stimps75,'UniformOutput',0); 




%% Read timing data
fid  = fopen([datdir,'/stim_data/75_-_track_7c_rf_mapping-2.din'],'r');
frt  = fread(fid,inf,'uint64'); fclose(fid);
ntot=length(frt); 
rfrt = reshape(frt,2,ntot/2)';
rfrt(1:10,:)
ts0  = frt(1,1); %initial time stamp

tims      = (rfrt(:,1)-ts0)/1000; %zero initial time stamp
fonsets   = reshape(tims,4,size(tims,1)/4)'; 
seqstarts = [0; fonsets( (1:23)*6000+1 , 1)]/1000 %start times of each 2min sequence (24 of them)



%% read the spikes
spikefiles = dir([datdir,'/spk_data/*.spk']); 
ncells = length(spikefiles); 
rspks75 = cell(ncells,1); 
% Read spike times
for ispfile = 1:ncells
  fid = fopen([datdir,'/spk_data/',spikefiles(ispfile).name],'r');
  rawts = fread(fid,inf,'uint64'); fclose(fid);
  ts = (rawts - ts0)/1000; %subtract off initial time stamp
  rspks75{ispfile} = ts(ts > 0)/1000;
end

%parse spike times into separate 2m recording sequences
seqlen = dt*6000; %number of time samples per sequence 
spksegs75 = cell(ncells,nfiles);
for icell=1:ncells;
	tspks = rspks75{icell}; 
	for iep =1:nfiles; 
		tstart = seqstarts(iep); tend=tstart+seqlen; %start and stop timestamps of current sequence
		spksegs75{icell,iep} = tspks( tspks>tstart & tspks < tend) - tstart; %spikes from current cell in current sequence
	end; 
end; 

%cell array containing all spike times (concatenated) for each cell
spks75 = cell(ncells,1);
for icell=1:ncells; 
	spks =[]; 
	for iseq=1:nfiles; spks=[spks;spksegs75{icell,iseq}+(iseq-1)*seqlen]; end; 
    spks75{icell}=spks; 
end; 


figure; 
ecdf(spks75{1}); hold on; for icell =2:26; ecdf(spks75{icell}); end; hold off; 
vline((1:24)*6000*dt)

ecdf(rspks75{1}); hold on; for icell =2:26; ecdf(rspks75{icell}); end; hold off; 
vline(seqstarts); axis tight

plot(rspks75{1},1:length(rspks75{1})); hold on; 
   for icell =2:26; plot(rspks75{icell},1:length(rspks75{icell})); end; hold off;



% spksRec75b = []; 
% for tcell = 1:ncells; 
% 	tspks = [];
% 	for ifile =1:nfiles; tebo=ebounds{ifile};
% 		tspks=[tspks;spksegs{tcell,ifile}+tebo(1)];
% 	end;
% 	spksRec75b{tcell} = tspks;
% end; 


nspks = cellfun(@length,spks75)

%save stmsRec75.mat stims75
%save dstmsRec75.mat dstims75

save spksRec75.mat spikefiles spks75 dt %save spike file names, and the concatenated arrays of spike times for each cell

save dstimpsRec75.mat dstimps75 %save cell array of spatially down-sampled stimulus sequences
save spksegsRec75.mat spksegs75; %save segmented spike times for each cell


figure; 

%%
cnames   = [stimfiles]; cellfun(@disp,cnames)
nconds   = length(cnames);
keys = {'af','wn','pn','pt','ps','ns'}
rids = cellfun(@(tkey)...
    find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
    keys,'UniformOutput',0);

taf        = rids{1}; %cellfun(@disp,cnames(tpn))
twn        = rids{2};
tpn        = rids{3}; %cellfun(@disp,cnames(tpn))
tpt        = rids{4};
tps        = rids{5}; %cellfun(@disp,cnames(tps))
tns        = rids{6}; %cellfun(@disp,cnames(tns))
usets = [rids{3}(3:4) rids{4}(3:4) rids{5}(3:4) rids{6}(3:4)]; ctype='all'; 
nspks = cellfun(@(x)length(x),spksegs75(ccell,:));

ccell =5;
stimlen = 6000;
longspks = [];
longstim = [];
for cond = 1:length(usets)
   cur_cond = usets(cond);
   longstim = [longstim; dstimps75{cur_cond}];
   longspks = [longspks; (cond-1)*stimlen*dt+spksegs75{ccell,cur_cond}];
    
end

NT = size(longstim,1);
spikebin_ids = 1+floor(longspks/dt);
spikebins = hist(spikebin_ids,1:NT)';
iisp = find(spikebins > 0);
splen = length(iisp);
nsp = sum(spikebins(iisp));
SS = makeStimRows(longstim,flen,iisp);
STA = (spikebins(iisp)'*SS)'/nsp;


clear STA
flen = 6;
for cond = 1:24
        cur_stim = dstimps75{cond};
    cur_spike_times = spksegs75{ccell,cond};
    if length(cur_spike_times) > 10
    stimfiles{cond}
        spikebin_ids = 1+floor(cur_spike_times/dt);
%         spikebin_ids = ceil(6000*rand(size(cur_spike_times)));
        spikebins = hist(spikebin_ids,1:6000)';
%         spikebins = hist(spikebin_ids,1:250)';
        
        iisp = find(spikebins>0);
        splen = length(iisp);
        nsp = sum(spikebins(iisp));
%         Msz = splen*size(cur_stim,2)*flen;  % size of "full" spike-triggered stimulus matrix
        %     SO = makeStimRows(cur_stim,flen,1);  % Convert stimulus to matrix where each row is one stim
        %     ov_mean = mean(SO);
        SS = makeStimRows(cur_stim, flen, iisp);
        rowlen = size(SS,2);
        STA(:,cond) = (spikebins(iisp)'*SS)'/nsp;
        plotfilterbank(STA(:,cond),32,1:1024)
        pause
        close all
    end
end

% avg_sta = STA(:,usets)*nspks(usets)';
% figure
% plotfilterbank(avg_sta,32,1:1024);
% 
% 


%% make a little stimulus movie
% tstim = dstimps75{1}; tsdim = sqrt(size(tstim,2)) 
% aviobj = avifile('natstimmovie2.avi','fps',15);
% for tslice =5800+(1:100);
% 	tslice
% 	plot2drfmat(reshape(tstim(tslice,:),tsdim,tsdim),[-2.,2.]);
% 	frame = getframe(gca); aviobj = addframe(aviobj,frame);
% end;
% aviobj = close(aviobj);

%%
close all
tstim = dstimps30{5}; tsdim = sqrt(size(tstim,2))
cstd = std(tstim(:));
for tslice =1:1000;
	tslice
	plot2drfmat(reshape(tstim(tslice,:),tsdim,tsdim),[-3*cstd,3*cstd]);
    pause(.1)
end;



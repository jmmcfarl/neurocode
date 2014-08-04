

clear all; pars = load('~/Data/rust/infos/StandardParametersRust.mat'); 
datdir = '~/Data/rust/stcbar/DATA/'; 

cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles); 
allcells = cell(ncells,1); 

tcellid=17

for tcellid = 1:ncells;
	disp(sprintf('*** working on cell %i',tcellid))
	fname = cfiles(tcellid).name;
	eval(['load ',['~/Data/rust/stcbar/Data/',fname]])
    
	stimlen = size(stim,1); SDIM = size(stim,2); 
	
	spikes_per_frm = spikes_per_frm(:); %binned spike counts per stim frame
	nspikes        = sum(spikes_per_frm); %total N spikes
	rbins          = (find(spikes_per_frm>0.5)); % bins with responses
	flen           = pars.flen; %number of filter temporal lags
	offset         = flen + pars.hilen; %?? 
	nsp            = spikes_per_frm(rbins); %spikes per response frame

    %create a vector containing the bin number for each individual spike
	spikebins =[];
	for isp =1:length(rbins); spikebins= [spikebins;repmat(rbins(isp),nsp(isp),1)]; end;
   
	csbins = spikebins(spikebins>flen)-flen+1; %??

    %compute sta and stcs, and stimulus projection onto stc subspace
%     [sta, stcs, fstim, evs] = getSTCfilters(stim,spikebins,flen,npos,nneg)
    [sta, stcs, fstim] = getSTCfilters(stim,spikebins,flen,8,8);	

    %plot resultant sta and stc filters
	figure; plot1dfilterbank([sta',stcs],SDIM); 
	
    %find localized e and i filters and plot
	suse               = findLocalSubunits(stcs(:,1:8),7,SDIM);
	susi               = findLocalSubunits(stcs(:,9:16),7,SDIM);
	plot1dfilterbanks({stcs,[suse,susi]},SDIM); drawnow; 
	
	tcell.name = fname; 
	tcell.sdim = SDIM; 
	tcell.stlen = stimlen; 
	tcell.nspks = nspikes; 
	tcell.sta  = sta'; 
	tcell.stcs = stcs; 
	tcell.suse = suse; 
	tcell.susi = susi; 
	allcells{tcellid} = tcell; 
end; 

save('RustSTCsallcells.mat','allcells');


for tcellid = 1:ncells;
	disp(sprintf('*** working on cell %i',tcellid))
	fname = cfiles(tcellid).name;
	eval(['load ',['~/myrepo/workspace/DataRepository/rust/stcbar/Data/',fname]])

	stimlen = size(stim,1); SDIM = size(stim,2); 
	
	spikes_per_frm = spikes_per_frm(:);
	nspikes        = sum(spikes_per_frm); 
	rbins          = (find(spikes_per_frm>0.5)); % bins with responses
	flen           = pars.flen;
	offset         = flen + pars.hilen; 
	nsp            = spikes_per_frm(rbins);

	spikebins =[];
	for isp =1:length(rbins); spikebins= [spikebins;repmat(rbins(isp),nsp(isp),1)]; end;
	csbins = spikebins(spikebins>flen)-flen+1;

	[sta, stcs, fstim] = getSTCfilters2(stim,spikebins,flen,8,8);	
	figure; plot1dfilterbank([sta',stcs],SDIM); 
	
		
	tcell.name = fname; 
	tcell.sdim = SDIM; 
	tcell.stlen = stimlen; 
	tcell.nspks = nspikes; 
	tcell.sta  = sta'; 
	tcell.stcs = stcs; 
	allcells{tcellid} = tcell;
end; 

save('RustSTCsallcellsProjOut.mat','allcells');


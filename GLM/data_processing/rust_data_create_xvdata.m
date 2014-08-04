clear all; pars = load('~/Data/rust/infos/StandardParametersRust.mat');
datdir = '~/Data/rust/stcbar/DATA/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);


for stcfileno = 1:ncells;
    fname = cfiles(stcfileno).name;
    uid   = [fname(2:3),'-',fname(6:7)];
    % eval(['load ',['~/myrepo/workspace/DataRepository/rust/stcbar/Data/',fname]])
    eval(['load ',['~/Data/rust/stcbar/Data/',fname]])
    
    ustim = stim; stimlen = size(ustim,1); SDIM = size(ustim,2); pars.SDIM=SDIM;
    
    spikes_per_frm = spikes_per_frm(:); %binned spike counts
    nspikes        = sum(spikes_per_frm)
    rbins          = (find(spikes_per_frm>0.5)); % bins with responses
    offset = pars.flen + pars.hilen; flen = pars.flen; npos=6; nneg=6;
    nsp            = spikes_per_frm(rbins);
    
    %create spike bin vector
    spikebins =[];
    for isp =1:length(rbins); spikebins= [spikebins;repmat(rbins(isp),nsp(isp),1)]; end;
    csbins = spikebins(spikebins>flen)-flen+1;
    
    
    %% specify cross validation sequences
    nparts   = 9;
    partlen  = floor(stimlen/nparts); %length of each XV part
    
    %??
    xvlen    = 8000; %length of chunks for XV
    xvdelay  = floor((partlen-xvlen)/2); %length of training chunks
    
    %boundaries of parts
    pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';
    enparts  = size(pbounds,1);
    
    %boundaries of parts for XV
    xvbounds = [pbounds(:,1)+xvdelay,pbounds(:,1)+xvdelay+xvlen-1];
    
    %boundaries of parts for training (two sets before and after each XV
    %sequence)
    trbound1 = [pbounds(:,1),xvbounds(:,1)-1];
    trbound2 = [xvbounds(:,2)+1,pbounds(:,2)];
    
    
    %% stimuli + downsampled stimuli
    [xvstims, trstims] = deal(cell(enparts,1));
    [axvstims,atrstims] = deal([]);
    
    for ipart = 1:enparts;
        xvb            = xvbounds(ipart,:);
        fstim          = ustim(xvb(1):xvb(2),:); %current stimulus chunk for XV
        xvstims{ipart} = fstim; %separated stim matrices
        axvstims       = [axvstims;fstim];  %concatenated stim matrices
        
        trb1    = trbound1(ipart,:);   trstim1 = ustim(trb1(1):trb1(2),:);
        trb2    = trbound2(ipart,:);   trstim2 = ustim(trb2(1):trb2(2),:);
        trstims{ipart} = [trstim1; trstim2];
        atrstims       = [atrstims; trstims{ipart}];
    end;
    
    
    %% spikes ** just for this one cell here.
    ncells            = 1; cid=1;
    [trsbins,xvsbins]   = deal(cell(ncells,nparts));
    [atrsbins,axvsbins] = deal(cell(ncells,1));
    
    sbins =spikebins;
    ctrspks = [];
    for part = 1:enparts;
        xvb     = xvbounds(part,:)-[1,1];
        trb1    = trbound1(part,:)-[1,1];
        trb2    = trbound2(part,:)-[1,1];
        
        xvsbins{part} = isInside(sbins,xvb)  - xvb(1);
        axvsbins{cid}  = [axvsbins{cid}; (part-1)*xvlen+ xvsbins{cid,part}];
        
        trsbins{cid,part} = [          isInside(sbins,trb1) - trb1(1); ...
            xvdelay + isInside(sbins,trb2) - trb2(1)];
        
        atrsbins{cid}     = [atrsbins{cid} ; ...
            (part-1)*(partlen-xvlen) + trsbins{cid,part}];
    end;
    
    save(['XVTRData-',uid,'.mat'],'trstims','atrstims','trsbins','atrsbins',...
        'xvstims','axvstims','xvsbins','axvsbins','fname',...
        'xvlen','xvdelay','pbounds','xvbounds','trbound1','trbound2');
end;

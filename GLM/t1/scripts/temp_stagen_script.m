dt   = 19.9920031987/1000;
cd '/Users/James/James_scripts/GLM/t1/rec74_modfits/'
flen = 6;
nconds   = length(dstimps74)
cnames   = [stimfiles]; cellfun(@disp,cnames)

keys = {'ar','wr','pr','nr'}
rids = cellfun(@(tkey)...
	find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
	keys,'UniformOutput',0); 

tns        = rids{4}; cellfun(@disp,cnames(tns))
tpn        = rids{3}; cellfun(@disp,cnames(tpn))

tconds = [tpn,tns]; ctype='all'; 
for tcell = 1:25
    segstas = [];
    allstims = [];
    allspks = [];
    for ifile =1:16;
        ifile
        spks    = spksegs74{tcell,ifile};
        if length(spks) > 10
            tstim   = repmat(dstimps74{ifile},24,1);
            allstims = [allstims; tstim - repmat(mean(tstim),6000,1)];
            spkbs   = 1+floor(spks/dt); spkbs=spkbs(spkbs<size(tstim,1));
            allspks = [allspks; spkbs+(ifile-1)*6000];
            longsta = getSTA(tstim,spkbs,flen)';
        else
            longsta = zeros(1024*flen,1);
        end
        segstas = [segstas,longsta];
    end;
    nspks = cellfun(@(x)length(x),spksegs74(tcell,:));
    avg_sta = mean(segstas(:,tconds),2);
    figure; plotfilterbank([segstas(:,1:16) avg_sta],32,1:1024)
    for ii = 1:16
        subplot(17,flen,(ii-1)*flen+1)
        title(sprintf('Spk: %d',nspks(ii)));
        yl = sprintf('%s-%s',stimfiles{ii}(1:2),stimfiles{ii}(22:26));
        ylabel(yl);
    end
         subplot(17,flen,(17-1)*flen+1)
        title(sprintf('Spk: %d',sum(nspks(tconds))));
        ylabel('Avg');
       
        f1 = gcf;
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 60]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    pname = sprintf('cell%d_allstas',tcell);
    print('-dpng',pname);
    close all
    
end

%%
dt   = 19.9920031987/1000;
cd '/Users/James/James_scripts/GLM/t1/rec30_modfits/'

for tcell = 1:25
    segstas = [];
    for ifile =1:16;
        ifile
        spks    = spksegs30{tcell,ifile};
        if length(spks) > 10
            tstim   = dstimps30{ifile};
            
            spkbs   = 1+floor(spks/dt); spkbs=spkbs(spkbs<size(tstim,1));
            longsta = getSTA(tstim,spkbs,flen)';
        else
            longsta = zeros(1024*flen,1);
        end
        segstas = [segstas,longsta];
    end;
    nspks = cellfun(@(x)length(x),spksegs30(tcell,:));
    figure; plotfilterbank(segstas(:,1:16),32,1:1024)
    for ii = 1:16
        subplot(16,flen,(ii-1)*flen+1)
        title(sprintf('Spk: %d',nspks(ii)));
        yl = sprintf('%s-%s',stimfiles{ii}(1:2),stimfiles{ii}(22:26));
        ylabel(yl);
    end
    f1 = gcf;
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 60]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    pname = sprintf('cell%d_allstas',tcell);
    print('-dpng',pname);
    close all
    
end


%%
dt   = 19.9920031987/1000;
cd '/Users/James/James_scripts/GLM/t1/rec75_modfits/'
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
tconds = [tpn tpt tps tns];
tconds(tconds > 16) = [];
for tcell = 1:26
    segstas = [];
    allstims = [];
    allspks = [];
    for ifile =1:16;
        ifile
        spks    = spksegs75{tcell,ifile};
        if length(spks) > 10
            tstim   = dstimps75{ifile};
            allstims = [allstims; tstim - repmat(mean(tstim),6000,1)];          
            spkbs   = 1+floor(spks/dt); spkbs=spkbs(spkbs<size(tstim,1));
               allspks = [allspks; spkbs+(ifile-1)*6000];
         longsta = getSTA(tstim,spkbs,flen)';
        else
            longsta = zeros(1024*flen,1);
        end
        segstas = [segstas,longsta];
    end;
    nspks = cellfun(@(x)length(x),spksegs75(tcell,:));
    avg_sta = mean(segstas(:,tconds),2);
    figure; plotfilterbank([segstas(:,1:16) avg_sta],32,1:1024)
    for ii = 1:16
        subplot(17,flen,(ii-1)*flen+1)
        title(sprintf('Spk: %d',nspks(ii)));
        yl = sprintf('%s-%s',stimfiles{ii}(1:2),stimfiles{ii}(22:26));
        ylabel(yl);
    end
    subplot(17,flen,(17-1)*flen+1)
    title(sprintf('Spk: %d',sum(nspks(tconds))));
    ylabel('Avg');
    f1 = gcf;
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 60]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    pname = sprintf('cell%d_allstas',tcell);
    print('-dpng',pname);
    close all
    
end
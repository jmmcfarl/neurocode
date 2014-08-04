function MeanSpike = PlotMeanSpike(DATA, varargin)
    recalc = 0;
    clnum = DATA.currentcluster;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'cluster',4)
            j = j+1;
            clnum = varargin{j};
        elseif strncmpi(varargin{j},'recalc',4)
            recalc = 1;
        end
        j = j+1;
    end
    id = find(DATA.clst == clnum+1);
    nid = find(DATA.clst == 1);
    MeanSpike.addmean = DATA.addmean;
    if DATA.interactive >= 0
        SetFigure(DATA.tag.meanspike, DATA);
        subplot(1,2,1);
        hold off;
        AllV = GetAllV(DATA);
        MeanV = getappdata(DATA.toplevel,'MeanV');
    else
        AllV = DATA.AllV;
        MeanV = DATA.meanV;
    end
    if DATA.keepsmooth && ~isempty(id)
        lfpts = -20000:10:20000;
        sampleid = DATA.trigtimes{1}(id);
        samplenid = DATA.trigtimes{1}(nid);
        Vall = mygetappdata(DATA,'Vall');
        allid = RectAdd(lfpts, sampleid, 'range', [1 length(Vall.Vsmooth)]);
        MeanSpike.lfpmean = mean(Vall.Vsmooth(allid)+double(Vall.V(allid)),1);
        [st, tid] = ShuffleTimes(Vall.t(sampleid),DATA.Expt);
        sid  = FullVTimes2id(Vall,st);
        allid = RectAdd(lfpts, sid, 'range', [1 length(Vall.Vsmooth)]);
        MeanSpike.shufflemean = mean(Vall.Vsmooth(allid)+double(Vall.V(allid)),1);
        muid = RectAdd(lfpts, samplenid, 'range', [1 length(Vall.Vsmooth)]);
        MeanSpike.lfpumean = mean(Vall.Vsmooth(muid)+double(Vall.V(muid)),1);
    end
    if ~isfield(DATA,'MeanSpike') || recalc
        ms = mean(AllV(:,:,id),3);
        mu = mean(AllV(:,:,nid),3);
        if DATA.addmean == 0 %if 1, already added to AllV
            mum =mean(MeanV(:,nid)');
            msm = mean(MeanV(:,id)');
            for j = 1:DATA.nprobes
                ms(j,:) = ms(j,:) + msm;
                mu(j,:) = mu(j,:) + mum;
            end
        end
        MeanSpike.ms = ms;
        MeanSpike.mu = mu;
        xc = corrcoef(ms(:),mu(:));
        MeanSpike.muxc = xc(1,2);
    else
        if isfield(DATA.cluster,'MeanSpike')
            ms = DATA.cluster.MeanSpike.ms;
            mu = DATA.cluster.MeanSpike.mu;
        else
        ms = DATA.MeanSpike.ms;
        mu = DATA.MeanSpike.mu;
        end
    end

    chspk = DATA.probe(1)+ [-1:1]; 
    chspk = DATA.chspk;
    dj = 0;
    %do csd first, then can do dvdt to CSD
    if DATA.plotcsd
        ms = (diff(ms,2,1));
        mu = (diff(mu,2,1));
        csd = diff(AllV,DATA.plotcsd,1);
        chspk = chspk-1;
        dj = 1;
    end
    if DATA.plotdvdt
        ms = (diff(ms,1,2));
        mu = (diff(mu,1,2));
    end
    if DATA.interactive >= 0

        if DATA.plot.comparemean > 0
            C = mygetappdata(DATA,'Clusters');
            imagesc([1 size(ms,2)],[1 size(ms,1)],ms);
            title(sprintf('P%d Cl %d',DATA.probe(1),clnum));
            subplot(1,2,2);
            hold off;
            msa = C{DATA.plot.comparemean}.MeanSpike.ms;
            imagesc([1 size(msa,2)],[1 size(msa,1)],msa);
            title(sprintf('P%d Cl %d',DATA.plot.comparemean,1));
            return;
        elseif strcmp(DATA.plot.meantype,'sidebyside')
            imagesc([1 size(ms,2)],[1 size(ms,1)-2],ms(1:2:end,:));
            subplot(1,2,2);
            hold off;
            imagesc([1 size(ms,2)],[2 size(ms,1)-1],ms(2:2:end,:));
            return;
        else
            imagesc(ms);
        end
        for j = 1:length(DATA.triggerchan)
            p = DATA.triggerchan(j);
            xl = get(gca,'xlim');
            line([xl(1) xl(2)/5],[p p ],'color','k');
        end
    title(sprintf('P%d Cl %d',DATA.probe(1),clnum));
    if size(AllV,1)  == 1
        if DATA.keepsmooth && isfield(MeanSpike,'lfpmean')
            subplot(2,1,2);
            hold off;
            plot(MeanSpike.lfpmean);
            hold on;
            plot(get(gca,'xlim'),[0 0],'k:');
            plot(MeanSpike.lfpumean,'g');
            plot(MeanSpike.shufflemean,'r');
            subplot(2,1,1);
        else
            subplot(1,1,1);
        end
    else
        subplot(1,2,2);
        if 0  %might add this one day
            hold off;
            imagesc(MeanSpike.mu);
            return;
        end
    end
    
    

    hold off;
    end
    chspk = chspk(chspk >0 & chspk <= size(ms,1));
    voff = CalcVoffset(ms,DATA.chspk,0);
    voff = voff-voff(DATA.probe(1));
    voff = DATA.voffset-DATA.voffset(DATA.probe(1));
    for j = chspk
        if DATA.plotdvdt && DATA.plotcsd
                dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(diff(csd(j,:,id),1,2),[],3) var(diff(csd(j,:,nid),1,2),[],3)]));
        elseif DATA.plotdvdt
        dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(diff(AllV(j,:,nid),1,2),[],3) var(diff(AllV(j,:,id),1,2),[],3)]));
        elseif DATA.plotcsd
                dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(csd(j,:,id),[],3) var(csd(j,:,nid),[],3)]));
        else
        dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(AllV(j,:,nid),[],3) var(AllV(j,:,id),[],3)]));
        end
        dpscale = max(ms(:))./max(dp);
        if DATA.interactive >= 0
            plot(mu(j,:)+voff(j)/5,'color',[0.5 0.5 0.5]);
            hold on;
            plot(ms(j,:)+voff(j)/5,'r');
            if j == DATA.probe(1)
                plot(abs(dp).*dpscale,'m');
            else
                plot(abs(dp).*dpscale+voff(j)/5,'g');
            end
            text(size(ms,2),voff(j)/5,sprintf('%d',j+dj),'fontsize',DATA.gui.fontsize(1));
        end
        peaks = find(diff(sign(diff(abs(dp)))) > 0);
        [a,b] = sort(abs(dp(peaks)),'descend');
        MeanSpike.dpmax(j,1:length(b)) = peaks(b);
        MeanSpike.dp(j,:) = dp;
        if recalc
        MeanSpike.vdprime(j,:) = dp;
        end
    end
    if strcmp(DATA.plot.meantype,'dprimeimage')
        imagesc(ms-mu);
        subplot(1,2,2);
        hold off;
        if size(MeanSpike.dp,1) < size(ms,1)
            MeanSpike.dp(size(ms,1),:,:) = 0;
        end
        imagesc(abs(MeanSpike.dp));
    end
    if DATA.interactive >= 0
    set(gca,'xlim',[1 size(ms,2)+1]);
    end
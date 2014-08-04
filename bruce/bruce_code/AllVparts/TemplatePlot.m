function [out, TemplateUsed, DprimeUsed] = TemplatePlot(DATA, varargin)

    projectout = 0;
    calcdips = 0;
    usemean = 0;
    usestd = 0;
    %normalized correlation problem is htat get lots close to 1 (or
    %std of template, so its compresed at that end. Not good for GM fits
    normalize = 0;
    plottype = 1;
    if DATA.profiling > 1
        ts = now;
        profile on;
    end
    j = 1;
    while  j <= length(varargin)
        if strncmpi(varargin{j},'calcdip',7)
            calcdips = 1;
        elseif strncmpi(varargin{j},'nodip',4)
            calcdips = 0;
        elseif strncmpi(varargin{j},'noplot',6)
            plottype = 0;
        elseif strncmpi(varargin{j},'projectout',7)
            projectout = 1;
        elseif strncmpi(varargin{j},'recalc',6) %redo cal with 'templateused field
            usemean = 2;
        elseif strncmpi(varargin{j},'usemean',7)
            usemean = 1;
        elseif strncmpi(varargin{j},'stdtemplate',7)
            usestd = 1;
        end
        j = j+1;
    end

    AllV = GetAllV(DATA);
  
    DataClusters = mygetappdata(DATA,'Clusters');
    %default in case there is no template to use
    if length(DATA.cluster.spts) > size(DATA.StdTemplate,2)
        DATA.StdTemplate(end,length(DATA.cluster.spts)) = 0;
    end
    nprobes = size(AllV,1);
    othermeans{1} = repmat(DATA.StdTemplate(1,1:length(DATA.cluster.spts)),nprobes,1);
    othermeans{2} = repmat(DATA.StdTemplate(2,1:length(DATA.cluster.spts)),nprobes,1);
    othermeans{3} = repmat(DATA.StdTemplate(2,1:length(DATA.cluster.spts)),nprobes,1);
    if DATA.currentcluster > 1 && isfield(DATA.cluster,'next') && length(DATA.cluster.next) >= DATA.currentcluster-1
        C = DATA.cluster.next{DATA.currentcluster-1};
        C.cluster = DATA.currentcluster;
        if isfield(DATA.cluster,'MeanSpike')
            othermeans{1} = DATA.cluster.MeanSpike.ms;
        end
    else
        C = DATA.cluster;
        C.cluster = 1;
        if length(DATA.cluster.next) && isfield(DATA.cluster.next{1},'MeanSpike')
        othermeans{1} = DATA.cluster.next{1}.MeanSpike.ms;
        end
    end
    if C.probe(1) >1 && length(DataClusters) >= C.probe(1) && isfield(DataClusters{C.probe(1)-1},'MeanSpike')
        othermeans{2} = DataClusters{C.probe(1)-1}.MeanSpike.ms;
    end
    if ~isfield(C,'quick')
        C.quick = 0;
    end
    if isfield(C,'neednewtemplate') && C.neednewtemplate == 0
        SetFigureName(DATA.toplevel,sprintf('Calculating Template Scores for C%d...',DATA.currentcluster));
    else
        SetFigureName(DATA.toplevel,sprintf('Calculating New Templates for C%d...',DATA.currentcluster));
    end
    DATA.templatecluster = DATA.currentcluster;
    
    if C.quick || ~isfield(C,'MeanSpike')
        C.MeanSpike = PlotMeanSpike(DATA,'recalc');
        if C.cluster > 1
            DATA.cluster.next{C.cluster-1}= C;
        else
            DATA.cluster = C;
        end
    end
    DATA.usestdtemplates = usestd;
    id = find(DATA.clst == DATA.currentcluster+1);
    nid = find(DATA.clst == 1);
    nevents = DATA.nevents;
    %If a new cluster has stolen points from an old one, want to keep teh
    %old template scores for that old cluster (until its boudnary is
    %chaned)
    if isfield(C,'neednewtemplate') && C.neednewtemplate == 0
        SetFigureName(DATA.toplevel,sprintf('Calculating Template Scores for C%d...',DATA.currentcluster));
        ms = C.TemplateUsed;
        mu = C.MeanSpike.mu;
        if length(DATA.chspk) == size(C.mumeanUsed,1)
            mu(DATA.chspk,:) = C.mumeanUsed;
%        elseif max(DATA.chspk) == size(C.mumeanUsed,1)
%            mu(DATA.chspk,:) = C.mumeanUsed;
        end
    elseif usemean == 1 || isempty(id)
        SetFigureName(DATA.toplevel,sprintf('Calculating Templates from mean C%d...',DATA.currentcluster));
        ms = C.MeanSpike.ms;
        mu = C.MeanSpike.mu;
    else
        SetFigureName(DATA.toplevel,sprintf('Calculating New Templates for C%d...',DATA.currentcluster));
        for j = size(AllV,1):-1:1
            ms(j,:) = mean(AllV(j,:,id),3);
            mu(j,:) = mean(AllV(j,:,nid),3);
        end
    end
    if size(mu,1) < size(ms,1)
        a = size(mu,1);
        b = size(ms,1);
        mu(a+1:b,:) = 0;
    end
    if C.cluster > 1
        DATA.cluster.next{C.cluster-1}.neednewtemplate = 0;
    else
        DATA.cluster.neednewtemplate = 0;
    end

   for j = DATA.nprobes:-1:1
            mu(j,:) = mu(j,:)-mean(mu(j,:));
            ms(j,:) = ms(j,:)-mean(ms(j,:));
            %        xc(j,1,:) = squeeze(TemplateScores(j,1,:))./(DATA.spkvar(j,:)' .* std(ms(j,:)));
            %        xc(j,2,:) = squeeze(TemplateScores(j,2,:))./(DATA.spkvar(j,:)' .* std(mu(j,:)));
   end

   if isnan(sum(ms(:)))
       fprintf('Cannot calculate template for empty cluster\n');
       DATA.MeanSpike.ms = zeros(size(ms));
       TemplateUsed = DATA.MeanSpike.ms;
       DprimeUsed = [];
       set(DATA.toplevel,'UserData',DATA);
       return;
   end
    ispk = find(DATA.chspk == DATA.probe(1));
    if usestd
        j = DATA.probe(1);
        meanV = repmat(mean(AllV(j,:,:),2),[1  size(AllV,2) 1]); %mean for each spike
        TemplateScores(ispk,1,:) = DATA.StdTemplate(1,:) * squeeze(AllV(j,:,:) - meanV);
        TemplateScores(ispk,2,:) = DATA.StdTemplate(2,:) * squeeze(AllV(j,:,:) - meanV);
        TemplateScores(ispk,3,:) = squeeze(diff(DATA.StdTemplate(1,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
        TemplateScores(ispk,7,:) = squeeze(diff(DATA.StdTemplate(2,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
        TemplateScores(ispk,8,:) = squeeze(diff(DATA.StdTemplate(2,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
        TemplateUsed = DATA.StdTemplate;
        DprimeUsed = [];
    else
%TemplateScores(p,1) is meanspike
%TemplateScores(p,2) is mu
%TemplateScores(p,3) is sprimt weighted
        for k = length(DATA.chspk):-1:1
            j = DATA.chspk(k);
            meanV = repmat(mean(AllV(j,:,:),2),[1  size(AllV,2) 1]);
            TemplateScores(k,1,:) = squeeze(ms(j,:)) * squeeze(AllV(j,:,:) - meanV);
            TemplateScores(k,2,:) = squeeze(mu(j,:)) * squeeze(AllV(j,:,:) - meanV);
            for nt = 1:3
            if size(othermeans{nt},2) ==  size(AllV,2)
                TemplateScores(k,9+nt,:) = squeeze(othermeans{nt}(j,:)) * squeeze(AllV(j,:,:) - meanV);
            elseif size(othermeans{nt},2)>  size(AllV,2)
                TemplateScores(k,9+nt,:) = squeeze(othermeans{nt}(j,1:size(AllV,2))) * squeeze(AllV(j,:,:) - meanV);
            end
            end
            if length(id)
                dp(j,:) = (mean(AllV(j,:,id),3)-mean(AllV(j,:,nid),3))./sqrt(mean([var(AllV(j,:,nid),[],3); var(AllV(j,:,id),[],3)]));
            else
                dp(j,:) = DATA.cluster.MeanSpike.dp(j,:);
            end
            if normalize
                TemplateScores(k,1,:) = TemplateScores(k,1,:) ./ std(AllV(j,:,:));
                TemplateScores(k,2,:) = TemplateScores(k,1,:) ./ std(AllV(j,:,:));
            end
            TemplateScores(k,9,:) = sum(abs(repmat(ms(j,:)',1,nevents) -squeeze(AllV(j,:,:)-meanV)));
        end
        TemplateUsed = ms;
        DprimeUsed = dp(DATA.chspk,:);
        clear meanV;
        dpcrit = 1;
        for k = length(DATA.chspk):-1:1
            j = DATA.chspk(k);
            %        dp * squeeze(AllV(j,:,:) - repmat(mu(j,:),[1 1 DATA.nevents]));
            TemplateScores(k,3,:) = squeeze(dp(j,:)) * squeeze(AllV(j,:,:) - repmat(mu(j,:),[1 1 DATA.nevents]));
            id = find(abs(dp(j,:)) > dpcrit);
            if length(id)
                TemplateScores(k,4,:) = squeeze(dp(j,id)) * squeeze(AllV(j,id,:));
                TemplateScores(k,5,:) = squeeze(sign(dp(j,id))) * squeeze(AllV(j,id,:));
            else
                [a,b] = sort(dp(j,:),'descend');
                id = b(1:2);
                TemplateScores(k,4,:) = squeeze(dp(j,id)) * squeeze(AllV(j,id,:));
                TemplateScores(k,5,:) = squeeze(sign(dp(j,id))) * squeeze(AllV(j,id,:));
            end
            TemplateScores(k,6,:) = squeeze(diff(ms(j,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
            TemplateScores(k,19,:) = squeeze(diff(othermeans{1}(j,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
            if normalize
                TemplateScores(k,6,:) = TemplateScores(k,6,:) ./ std(diff(AllV(j,:,:),1,2));
            end

        end
        if size(AllV,1) > 2 %look at spatial derivatives
            if max(DATA.chspk) < nprobes
                dyspk = [DATA.chspk DATA.chspk(end)+1];
            elseif min(DATA.chspk) <= 1
                dyspk = [1:length(DATA.chspk)+1];
            else
                dyspk = [nprobes-length(DATA.chspk):nprobes];
            end
            dyspk = dyspk(dyspk > 0);
            dv = diff(AllV(dyspk,:,:),1,1);  %diff along length
            mdv = diff(ms(dyspk,:),1,1);
            for k = 1:size(dv,1)
                TemplateScores(k,7,:) = squeeze(mdv(k,:)) * squeeze(dv(k,:,:));
                if normalize
                    TemplateScores(k,7,:) = TemplateScores(k,7,:) ./ std(squeeze(mdv(k,:)) * squeeze(dv(k,:,:)));
                end

            end
            clear dv;
            clear mdv;

            if min(dyspk) > 1
                csdspk = [dyspk(1)-1 dyspk];
            else
                csdspk = [1:length(dyspk)+1];
            end
            csdspk = csdspk(csdspk > 0 & csdspk <= nprobes);

            csd = diff(AllV(csdspk,:,:),2,1);
            mcsd = diff(ms(csdspk,:),2,1);
            for j = 1:size(csd,1)
                TemplateScores(j,8,:) = squeeze(mcsd(j,:)) * squeeze(csd(j,:,:));
            end
            clear csd;
            clear mcsd;
        else %Temporary - fill 7 adn 8 with somethins
            if size(AllV,2) >  size(DATA.StdTemplate,2)
                vpts = 1:size(DATA.StdTemplate,2);
            else
                vpts = 1:size(AllV,2);
            end
            meanV = repmat(mean(AllV(DATA.chspk,:,:),2),[1  size(AllV,2) 1]); %mean for each spike
            TemplateScores(1,7,:) = squeeze(diff(mu(j,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2)); %mu dvdt
%            TemplateScores(1,7,:) = DATA.StdTemplate(1,vpts) * squeeze(AllV(DATA.chspk,vpts,:) - meanV);
            TemplateScores(1,8,:) = DATA.StdTemplate(2,vpts) * squeeze(AllV(DATA.chspk,vpts,:) - meanV);
            TemplateScores(1,12,:) = squeeze(diff(DATA.StdTemplate(1,vpts),1,2)) * squeeze(diff(AllV(j,vpts,:),1,2));
%13-18 also free            
        end
    end
    DATA.DprimeUsed = DprimeUsed;
    DATA.TemplateUsed = TemplateUsed;
    DATA.mumeanUsed = mu(DATA.chspk,:);
    DATA.Template.othermeans = othermeans;
    %These shouls be copied into cluster ONLY if this space is used for
    %classification, so don't it here
    if 0 
    if DATA.currentcluster > 1 && isfield(DATA.cluster,'next')
        DATA.cluster.next{DATA.currentcluster-1}.TemplateUsed = TemplateUsed;
        DATA.cluster.next{DATA.currentcluster-1}.DprimeUsed = DprimeUsed;
        DATA.cluster.next{DATA.currentcluster-1}.mumeanUsed = mu(DATA.chspk,:);
    elseif DATA.currentcluster == 1
        DATA.cluster.TemplateUsed = TemplateUsed;
        DATA.cluster.DprimeUsed = DprimeUsed;
        DATA.cluster.mumeanUsed = mu(DATA.chspk,:);
    end
    end
    if length(DATA.chspk) > 2
    chspk = DATA.chspk;
    else
    chspk = DATA.probe(1)+ [-1:1];
    end
    chspk = chspk(chspk >0 & chspk <= nprobes);
    if min(chspk) > 1
    csdspk = chspk-1;
    else
        csdspk = chspk;
    end
    
    if max(chspk)  == nprobes
        xspk = min(chspk)-1;
    else
        xspk = max(chspk)+1;
    end
        
    if projectout  %doesn't seem much use...., and uses too much memory
%project out the templates, then redo the pca
    mg = sum(ms.^2,2);
    G = sum(AllV.*repmat(ms,[1 1 DATA.nevents]),2)./repmat(mg,[1 1 DATA.nevents]);
    nv = AllV - repmat(ms,[1 1 DATA.nevents]) .* repmat(G,[1 size(AllV,2) 1]);
    TV = nv(chspk(1),:,:);
    for j = 2:length(chspk)
        TV = cat(2,TV,nv(chspk(j),:,:));
    end
    TV = squeeze(TV)';
    [pc, E] = eig(cov(TV));
    pc = fliplr(pc); %put largest first;
    pcs = TV*pc;
    end
    ispk = find(DATA.chspk == DATA.probe(1));

    if projectout == 2
        TMPL.pcs(:,1) = sum(TemplateScores(:,1,:));
        TMPL.pcs(:,2:9) = pcs(:,1:8);
    elseif usestd
        TMPL.pcs(:,1) = TemplateScores(ispk,1,:);  %1r
        TMPL.pcs(:,2) = TemplateScores(ispk,2,:); %
        TMPL.pcs(:,3) = TemplateScores(ispk,3,:);
        TMPL.pcs(:,4) = TemplateScores(ispk,7,:); 
        TMPL.pcs(:,5) = TemplateScores(ispk,1,:);
        TMPL.pcs(:,6) = TemplateScores(ispk,2,:);
        TMPL.pcs(:,7) = TemplateScores(ispk,3,:);
        TMPL.pcs(:,8) = TemplateScores(ispk,7,:); %2dt 
        TMPL.pcs(:,9) = TemplateScores(ispk,1,:);
        TMPL.pcs(:,10) = TemplateScores(ispk,7,:); %2dt
        TMPL.pcs(:,11) = TemplateScores(ispk,3,:);
        TMPL.pcs(:,12) = TemplateScores(ispk,7,:); 
    elseif length(DATA.chspk) == 1
        TMPL.pcs(:,1) = TemplateScores(ispk,1,:);
        TMPL.pcs(:,2) = TemplateScores(1,2,:);
        TMPL.pcs(:,8) = TemplateScores(ispk,6,:); %dvdt
        TMPL.pcs(:,3) = TemplateScores(ispk,7,:); %std 1
        TMPL.pcs(:,4) = TemplateScores(ispk,8,:); %std 2
        TMPL.pcs(:,6) = TemplateScores(ispk,3,:); %dprime weighted
        TMPL.pcs(:,7) = TemplateScores(1,3,:); %dprime weighted
        TMPL.pcs(:,3) = TemplateScores(1,7,:); %mu dvdt
        TMPL.pcs(:,4) = TemplateScores(1,19,:); %dvdt for other cluster
        TMPL.pcs(:,10) = TemplateScores(1,3,:); %dprime
        TMPL.pcs(:,11) =DATA.rV;
        TMPL.pcs(:,12) = (TemplateScores(1,12,:)); %std1 dvdt
        TMPL.pcs(:,13) = (TemplateScores(1,9,:)); %sum abs diffs
        TMPL.pcs(:,14) = (TemplateScores(1,2,:)); %sum mu score
        TMPL.pcs(:,15) = (TemplateScores(1,10,:)); %sum score for other template
        TMPL.pcs(:,9) = TemplateScores(1,10,:); %Other Cluster
        TMPL.pcs(:,5) = (TemplateScores(1,10,:)) - TemplateScores(1,1,:); %diff in sums
        TMPL.pcs(:,16) = (TemplateScores(1,10,:)) - TemplateScores(1,2,:); %diff in sums
        TMPL.pcs(:,17) = (TemplateScores(1,11,:)); %TemplateScore for cluster above
        TMPL.pcs(:,18) = (TemplateScores(1,12,:)); %TemplateScore for cluster below
    else
        TMPL.pcs(:,1) = TemplateScores(ispk,1,:);
        TMPL.pcs(:,2) = sum(TemplateScores(:,1,:)); %sum of r
        TMPL.pcs(:,8) = TemplateScores(ispk,6,:); %dvdt
        TMPL.pcs(:,9) = TemplateScores(ispk,7,:); %dvdy
        TMPL.pcs(:,5) = TemplateScores(ispk,8,:); %csd
        TMPL.pcs(:,6) = TemplateScores(ispk,3,:); %dprime weighted
        TMPL.pcs(:,7) = TemplateScores(1,3,:); %dprime weighted
%This needs work if length(chspk) > 3
        if length(chspk) > 1
            TMPL.pcs(:,3) = TemplateScores(1,1,:);
        end
        if usestd
            TMPL.pcs(:,4) = TemplateScores(ispk,1,:);
        elseif length(chspk) > 2
            if ispk == 3
            TMPL.pcs(:,4) = TemplateScores(2,1,:);
            else
            TMPL.pcs(:,4) = TemplateScores(length(chspk),1,:);
            end
        else
            TMPL.pcs(:,4) = TemplateScores(end,1,:);
        end
        TMPL.pcs(:,10) = sum(TemplateScores(:,6,:));
        if DATA.trigdt == 4
            TMPL.pcs(:,11) =DATA.rV;
        else
            TMPL.pcs(:,11) = sum(TemplateScores(:,3,:)); %dprime weighted
        end
        TMPL.pcs(:,12) = sum(TemplateScores(:,7,:)); %sum dy
        TMPL.pcs(:,13) = sum(TemplateScores(:,9,:)); %sum abs diffs
        TMPL.pcs(:,14) = sum(TemplateScores(:,2,:)); %sum mu score
        TMPL.pcs(:,15) = sum(TemplateScores(:,10,:)); %sum score for other template
        TMPL.pcs(:,16) = sum(TemplateScores(:,10,:)) - sum(TemplateScores(:,2,:)); %diff in sums
        TMPL.pcs(:,17) = sum(TemplateScores(:,11,:)); %TemplateScore for cluster above
        TMPL.pcs(:,18) = sum(TemplateScores(:,12,:)); %TemplateScore for cluster below
        if projectout
            TMPL.pcs(:,13:15) = pcs(:,1:3);
        end
    end
    DATA.TemplateLabels = TemplateLabels(DATA, usestd);
    TMPL.dvdt = 0;
    TMPL.csd = 0;
    TMPL.clid = id;
    TMPL.nid = nid;
    TMPL.toplevel = DATA.toplevel;
    TMPL.pcplots = DATA.pcplots;
    TMPL.clplot = DATA.clplot;
    TMPL.plottype = 1;
    DATA.TemplateScores = TMPL.pcs;
    for j = 1:size(DATA.TemplateScores,2)
        sds(j) = std(DATA.TemplateScores(:,j));
        if sds(j) > 0
        DATA.TemplateScores(:,j) = DATA.TemplateScores(:,j)./sds(j);
        DATA.TemplateScaling(j) = sds(j);
        end
    end
    
    if calcdips
        for j = 1:5
            TMPL.dipvals(j) = HartigansDipTest(sort(TMPL.pcs(:,j)));
        end
        DATA.tmpdips = CalculateTemplateDips(DATA);
        theta = 0:pi/36:pi * 35/36;
        for j = 1:length(theta)
            for k = 2:5;
                xy = xyrotate(TMPL.pcs(:,1),TMPL.pcs(:,k),theta(j));
                rdip(j,k) = HartigansDipTest(xy(:,1));
            end
        end
    else
        DATA.tmpdips = zeros(1,8);
    end

    if ~ismember(DATA.plottype,[3 4])
    DATA.plottype = 3;
    end

    if DATA.watchplots  && plottype == 1
    SetFigure(DATA.tag.tmplscore, DATA);
    subplot(1,1,1);
    PlotTemplateScores(TMPL,TemplateScores, DATA.chspk);
%   plot(DATA.TemplateScores(:,1),DATA.spkvar(:,1));

    DATA = ReplotPCs(DATA,[]);
    end
%currently sets DATA.TemplateScores, DATA.TemplateLabels, and DATA.tmpdips
%if output is requestsed, don't set the figure userdata
    if DATA.profiling > 1
       fprintf('Templates tookd %.2f\n',mytoc(ts));
       profile viewer;
    end
    if nargout
        set(DATA.toplevel,'Name',get(DATA.toplevel,'Tag'));
        out = DATA.TemplateScores;
        return;
    end
    set(DATA.toplevel,'UserData',DATA);
    SetFigureName(DATA.toplevel,get(DATA.toplevel,'Tag'));

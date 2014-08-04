function [distance, obj, xy, details] = BestSpace(DATA, varargin)

newtemplate = 0;
nloops = 0; %to test multiple fits for consistency
pconly = 0;
nr = 1;
ntr = 1;
nc = 2;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'newtemplate',7);
        newtemplate = 1;
        nr=5;
    elseif strncmpi(varargin{j},'template',5)
        ntr=5; %redoing a template fit. Make sure its right.
    elseif strncmpi(varargin{j},'pconly',5)
        pconly = 1; %redoing a template fit. Make sure its right.
    elseif strncmpi(varargin{j},'ncells',5)
        j = j+1;
        nc = varargin{j};
    end
    j = j+1;
end
p= DATA.probe(1);
uid = DATA.uid;
%when doing this from scratch (remaking template from other cut), use
%several starts to make sure the PC cut is best.
%Now uses GMfit that set a sensible start point.

[P, distance(1), a] = GMfit(DATA.pcs(uid,DATA.pcspace),nc,1);
details.pctook = a.took;
objs{1} = P;
details.Converged(1) = P.Converged;
ll(1) = P.NlogL;
if pconly
    best = 1;
    obj = objs{best};
    [xy, details.cid] = ProjectND(DATA, best, objs{best});
    e(1) = mean(DATA.energy(1,details.cid == 1));
    e(2) = mean(DATA.energy(1,details.cid == 2));
    if e(1) > e(2) && nc < 3
        details.cid = 3 - details.cid;
    end
    return;
end
%The Var/Energy plot is particularly sensitive to electdoe movement
%artifacts, so exclude spikes after trial end when making the cut
AllV = GetAllV(DATA);
iid = setdiff(uid, DATA.postevents);
xy(:,1) = DATA.energy(1,iid);
xy(:,2) = DATA.spkvar(p,iid)./DATA.energy(1,iid);

[E, distance(2)] = GMfit(xy,nc,1);
details.Converged(2) = E.Converged;
objs{2} = E;
distance(2) = gmdistance(E);
ll(2) = E.NlogL;

try
    [V, distance(3), a] = GMfit(squeeze(AllV(p,DATA.vspace,uid))',nc,1);
catch
    cprintf('errors','ERROR!!! E%dP%d uid %d-%d of %d, clusters %s\n',DATA.exptno,p,min(uid),max(uid),size(AllV,3),nc);
end
objs{3} = V;
details.Converged(3) = V.Converged;

ll(3) = V.NlogL;
quick = 1;
if newtemplate || ~isfield(DATA,'TemplateScores')
    [d, best] = max(distance);
    [xy, details.cid] = ProjectND(DATA, best, objs{best});
    %        [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
    %could use gmdistribution/cluster to assign to groups here.
    if quick & length(unique(details.cid)) > 1
        DATA.clid = find(details.cid  == 2);
        DATA.nid = find(details.cid  == 1);
        DATA.clst = details.cid;
    else
        [cluster.crit, a] = GMDip(xy(uid,:),DATA.energy(1,uid),'label',DATA.idstr);
        cluster.sign = a.sign;
        if a.sign >= 0
            DATA.clid = find(xy(:,1) > cluster.crit(1));
            DATA.nid = find(xy(:,1) <= cluster.crit(1));
        else
            DATA.clid = find(xy(:,1) < cluster.crit(1));
            DATA.nid = find(xy(:,1) >= cluster.crit(1));
        end
        DATA.clst(DATA.clid) = 2;
        DATA.clst(DATA.nid) = 1;
    end
    DATA.MeanSpike = PlotMeanSpike(DATA,'recalc');
    if DATA.currentcluster == 1
        DATA.cluster.MeanSpike = DATA.MeanSpike;
        DATA.cluster.spts = DATA.spts;
    else
        DATA.cluster.next{DATA.currentcluster-1}.MeanSpike = DATA.MeanSpike;
    end
    TemplatePlot(DATA,'nodip','usemean');
    DATA = get(DATA.toplevel,'UserData');
    %        imean = DATA.MeanSpike;  seems unused
    ntr = 1;
end
if isfield(DATA,'TemplateScores')
    
    [objs{4}, distance(4), a] = GMfit(DATA.TemplateScores(uid,DATA.tmplspace(1,:)),nc,ntr);
    %        objs{5} = gmdistribution.fit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),3,'Options',statset('MaxIter',1000));
    distance(4) = gmdistance(objs{4});
    ll(4) = objs{4}.NlogL;
    details.Converged(4) = objs{4}.Converged;
    
    tic;
    for j = 1:nloops
        T{j}= gmdistribution.fit(DATA.TemplateScores(uid,DATA.tmplspace(1,:)),nc,'Options',statset('MaxIter',1000));
        ll(j) = T{j}.NlogL;
    end
    if nloops
        toc
        tic;
        Tn= gmdistribution.fit(DATA.TemplateScores(uid,DATA.tmplspace(1,:)),nc,'Options',statset('MaxIter',1000),'Replicates',nloops);
        toc;
    end
    %        [a,b,c] = cluster(objs{4},DATA.TemplateScores(:,DATA.tmplspace(1,:)));
end

bestll = min(ll);
[d, best] = max(distance);
if ismember(best,[2 3])
    safeid = [1 4];
    [d, best] = max(distance(safeid));
    best = safeid(best);
    err = sprintf('Using Space %d (%.1f, %.1f) not Space 3 (%.1f. %.1f)\r\n',best,distance(best),ll(best),distance(3),ll(3));
    fprintf('%s',err);
    if DATA.logfid > 2
        fprintf(DATA.logfid,'%s',err);
    end
end
details.bestd = d;
details.besti = best;
if best == 4
    DATA.cluster.space = [6 4];
    DATA.cluster.shape = 2;
    if ~isfield(DATA,'clid') || isempty(DATA.clid)
        xy = ProjectND(DATA, best, objs{best});
        %        [cluster.crit, details] =
        %        FindDip(xy(:,1),DATA.energy(1,:));
        [cluster.crit, details] = GMDip(xy(uid,:),DATA.energy(1,uid),'label',DATA.idstr);
        cluster.sign = details.sign;
        if details.sign >= 0
            DATA.clid = find(xy(:,1) > cluster.crit(1));
            DATA.nid = find(xy(:,1) <= cluster.crit(1));
        else
            DATA.clid = find(xy(:,1) < cluster.crit(1));
            DATA.nid = find(xy(:,1) >= cluster.crit(1));
        end
    end
    olddistance = distance;
    
    objs{4} = IterateTemplateFit(DATA, objs{best});
    distance(4) = gmdistance(objs{4});
    DATA = get(DATA.toplevel,'UserData');
end



details.TemplateUsed = DATA.TemplateUsed;

obj = objs{best};
[xy, details.cid] = ProjectND(DATA, best, objs{best});
e(1) = mean(DATA.energy(1,details.cid == 1));
e(2) = mean(DATA.energy(1,details.cid == 2));
if e(1) > e(2)
    details.cid = 3 - details.cid;
end


%    [a,b,c] = BestAngle(xy(:,1),xy(:,2), 3); %Should not be necessary.
details.ll = ll;
%    details.bestangle = a;
%    xy = xyrotate(xy(:,1),xy(:,2),a);


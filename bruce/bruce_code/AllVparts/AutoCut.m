function [E, Scores, tbs, xy, details]  = AutoCut(DATA, varargin)
usev = 0;
refine = 0;
usegm = 0;
tbs = [];
Scores = [];
newtemplate = 0;
newDATA = 0;

if strcmp(DATA.autocutmode,'mahal')
    usegm = 1;
end
E.cutmode = DATA.autocutmode;
E.newscores = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'refine',4)
        refine  = 1;
    elseif strncmpi(varargin{j},'mahal',4)
        usegm = 1;
        E.cutmode = 'mahal';
    end
    j = j+1;
end
E.manual = 0;
AllV = GetAllV(DATA);
ctime = now;
if usegm
    DATA.cluster.spts = DATA.spts;
    DATA.cluster.probe = DATA.probe;
    if ~isfield(DATA.cluster,'next')
        DATA.cluster.next = {};
    end
    for j = 1:length(DATA.ncomponents)
        [bd,obj{j},xy,c] = BestSpace(DATA,'newtemplate','ncells',DATA.ncomponents(j));
    end
    E.newscores = 1;
    DATA.ndxy = xy;
    [a,b] = max(bd);
    E.bestspace = [c.bestd c.besti];
    E.bestd = c.besti;
    E.bestll = c.ll;
    E.bestcl = c.cid;
    E.bestfitres = c.Converged;
    if isfield(c,'TemplateUsed')
        E.TemplateUsed = c.TemplateUsed;
    end
    E.auto = 1;
    if c.besti == 4
        [gm,d]= max(bd(1:3));
        s = sprintf('GM %.2f for space %d (%.2f for %d)dt%d csd%d',a,c.besti,gm,d,DATA.dvdt,DATA.csd);
    else
        s = sprintf('GM %.2f for space %s dt%d csd%d',a,DATA.SpaceTypes{b},DATA.dvdt,DATA.csd);
    end
    fprintf('%s\n',s);
    PrintMsg(DATA.logfid,'P%d (S%d) %s %d events at %s\r\n',DATA.probe(1),DATA.savespikes,s,DATA.nevents,datestr(now));
    %    [dip, details] = FindDip(xy(:,1),DATA.energy(1,:),'gmix');
    [dip, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
    [E.gmfit2d, d] = GMfit(xy,2,1,'idlist',c.cid);
    E.mahal(4) = d;
    E.angle = 0;
    if d > 2  && d > details.gmdprime * 1.1 %% 2D significantly better - check rotation
        [a, gm,  dipres] = BestAngleGM(xy, E.gmfit2d, details);
        if gm > 2 && gm > details.gmdprime * 1.5
            %           xy = dipres.xy; %rotated values
            details.dipres = dipres.gmfit;
            details.gmdprime = gm;
            E.angle = a;
            [dip, details] = GMDip(dipres.xy,DATA.energy(1,:),'label',DATA.idstr);
        end
    end
    E.mahal(4) = details.gmdprime;
    details.pctook = c.pctook;
    
    E.space = [6 c.besti];
    E.bestd = bd;
    E.pcplot = [2 8 10 11 12];
    crit = dip(1);
    E.gmdip = dip;
    E.xyr = [dip(1) 0];
    E.shape = 2;
    E.sign = details.sign;
    %in case rotation was applied by BestAngle above
    x = xyrotate([crit crit],[min(xy(:,2)) max(xy(:,2))],-E.angle);
    E.pos = x([1 3 2 4]);
    %    E.pos = [crit min(xy(:,2)) crit max(xy(:,2))];
    E.plottype = DATA.plottype;
    E.gmfit = obj{end};
    if length(obj) > 1
        E.gmfits = obj(1:end-1);
    end
    E.gmfit1d = details.G{details.best};
    E.gmdipres = details.dipres;
    E.gmdprime = details.gmdprime;
    E.autodipsize = details.dipsize;
    E.dipsize = details.cdipsize;
    details.newDATA = 1;
    return;
end
p = DATA.pcplots;
for j =1:length(p)
    [as(j),bs(j),c] = BestAngle(DATA.pcs(:,p(j,1)),DATA.pcs(:,p(j,2)),1);
    gd(j) = c.mahal;
end
n = length(as);
[x,y] = GetClusterXYData(DATA,[]);
n = n+1;
vare = n;
[as(vare), bs(vare),c] = BestAngle(x,y,1);
bs(vare) = bs(vare).* 0.7;  %% only use varE if substantially better
gd(vare) = c.mahal;
n = n+1;

[bd,obj,xy] = BestSpace(DATA);
bs(n) = BimodalCoeff(xy(:,1));
[gd(n), besttype] = max(bd);
as(n) = 0;
E.bestspace(1) = bs(n);
E.bestspace(2) = gd(n);
E.bestd = bd;
n = n+1;


if usev
    p = DATA.vpts;
    for j =1:length(p)
        [as(j+n),bs(j+n),c] = BestAngle(AllV(p(j,1),p(j,2),:),AllV(p(j,3),p(j,4),:),2);
        bs(j+n) = c.bmc(c.besti);
    end
end
[a,j] = max(bs);
pcbii = a;
cluster.angle = as(j);
if j == vare %x,y already made
    p = DATA.vpts;
    cluster.space = [];
    E.plottype = DATA.plottype;
elseif j == vare+1  %BestSpace
    p = DATA.vpts;
    cluster.space = [6 besttype];
    cluster.angle = 0;
    E.plottype = DATA.plottype;
    E.shape = 2;
    E.space = cluster.space;
    x = xy(:,1);
    y = xy(:,2);
elseif j > size(DATA.pcplots,1)
    j = j-8;
    p = DATA.vpts;
    cluster.space = [2 DATA.vpts(j,:)];
    x = AllV(p(j,1),p(j,2),:);
    y = AllV(p(j,3),p(j,4),:);
    E.plottype = 2;
else
    p= DATA.pcplots;
    cluster.space = [1 DATA.pcplots(j,:)];
    x = DATA.pcs(:,p(j,1));
    y = DATA.pcs(:,p(j,2));
    E.plottype = 1;
end
xy = xyrotate(x,y,cluster.angle);
%    [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
[cluster.crit, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
cluster.sign = details.sign;
if refine
    if j == 9  %used var/e
        bettercrit = 1;
    else
        bettercrit = 1.2; %prefer PC cut if equal
    end
    if details.sign >= 0
        DATA.clid = find(xy(:,1) > cluster.crit(1));
        DATA.nid = find(xy(:,1) <= cluster.crit(1));
    else
        DATA.clid = find(xy(:,1) < cluster.crit(1));
        DATA.nid = find(xy(:,1) >= cluster.crit(1));
    end
    [Scores, details.TemplateUsed, details.DprimeUsed] = TemplatePlot(DATA);
    E.newscores = 1;
    cluster.firstspace = cluster.space;
    cluster.firstbmi = bs(j);
    if length(bd) < 4  %no template for the BestSpace calc; Do again
        DATA.TemplateScores = Scores;
        [bd, obj, bxy] = BestSpace(DATA);
        [gd(vare+1), besttype] = max(bd);
        bmi =  BimodalCoeff(bxy(:,1));
        E.bestspace(1) = bmi;
        E.bestspace(2) = gd(vare+1);
        E.bestd = bd;
        bs(vare+1) = bmi;
        if bmi > pcbii
            xy = bxy;
            cluster.space = [6 besttype];
            cluster.angle = 0;
            E.plottype = DATA.plottype;
            E.shape = 2;
            E.space = cluster.space;
            x = xy(:,1);
            y = xy(:,2);
            pcbii = bmi;
            %                [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
            [cluster.crit, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
            cluster.sign = details.sign;
        end
    end
    p = DATA.tmplots;
    for j =1:8
        [tas(j),tbs(j),c] = BestAngle(Scores(:,p(j,1)),Scores(:,p(j,2)),1);
        tgd(j) = c.mahal;
    end
    [a,j] = max(tbs);
    if a > pcbii * bettercrit
        cluster.angle = tas(j);
        cluster.firstspace = cluster.space;
        cluster.space = [3 p(j,:)];
        x = Scores(:,p(j,1));
        y = Scores(:,p(j,2));
        E.plottype = 3;
        E.shape = 1;
        xy = xyrotate(x,y,cluster.angle);
        %            [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
        [cluster.crit, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
        cluster.sign = details.sign;
    end
end
cluster.auto = 1;
E.autotook = mytoc(ctime);
if length(cluster.space) > 1 && cluster.space(1) == 6 && cluster.space(2) == 4
    E.plottype = 3;
end
details.newDATA = newDATA;
E = BoundaryFromCluster(E,cluster, DATA.currentcluster);

function [AllClusters, AllFullVData, details] = LoadCluster(dirname, eid, varargin)
%LoadCluster(dirname, expts, ...)
%Load Cluster info for an expt or list of expts
%if length(expts) = 1 returns a 1 x nprobes cell matrix, else 
%a 1 x length(expts) matrix each containing 1 x nprobes cells
%
%LoadCluster Combines info from Clusters and ClusterDetails files
%so that clst, t and Evec fields are in clusters
%LoadCluster(dirname, expts, 'gextxy') also inlcudes xy
%LoadCluster(dirname, expts, 'rawxy') also inlcudes xy with any roation
%removed (PlotClusters needs it this way
%LoadCluster(dirname, expts, 'alltimes') 
%replaces Clusters{}.times with ClusterDetails{}.t, saving memory
%(Also needed by PlotClusters)

AllClusters = {};
AllFullVData = {};
getauto = 0;
f = {'Evec' 'clst' 't'}; %fields to copy

if ~isdir(dirname) && exist(dirname) %given filename
    name = dirname;
    dirname = fileparts(name);
    if nargin > 1
        varargin = {eid varargin{:}};
    end
    eid = GetExptNumber(name);
end



rawxy = 0;
alltimes = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'auto',5)
        getauto = 1;
    elseif strncmpi(varargin{j},'plain',5)
        plain = 1;
    elseif strncmpi(varargin{j},'noauto',5)
        getauto = -1;
    elseif strncmpi(varargin{j},'getxy',5)
        f = {f{:} 'xy' 'triggerV'};
    elseif strncmpi(varargin{j},'rawxy',5)
        rawxy = 1;
        f = {f{:} 'xy' 'triggerV'};
    elseif strncmpi(varargin{j},'alltimes',5)
        alltimes = 1;
    end
    j = j+1;
end
details = {};

if ispc && dirname(1) == '/'
    dirname(1) = '\';
end
for j = 1:length(eid)
    ts = now;
    e = floor(eid(j));
    exptnos(j) = e;
    if round(rem(eid(j),e).*10) == 1
        xs = 'a';
    else
        xs = '';
    end
    name = [dirname '/Expt' num2str(e) xs 'ClusterTimes.mat'];
    dname = [dirname '/Expt' num2str(e) xs 'ClusterTimesDetails.mat'];
    daname = [dirname '/Expt' num2str(e) xs 'AutoClusterTimesDetails.mat'];
    aname = [dirname '/Expt' num2str(e) xs 'AutoClusterTimes.mat'];
    if getauto > 0
        dname = daname;
        name = aname;
    end
        
    if exist(aname,'file') && getauto >= 0
        tic;
        load(aname);
        details{j}.loaddur = toc;
        AutoClusters = Clusters;
        for p = 1:length(AutoClusters)
            AutoClusters{p}.auto = 1;
        end
        if exist('FullVData','var')
            AllFullVData{j} = FullVData;
        end
        details{j}.loadname = aname;
    else 
        AutoClusters = {};
    end
    if exist(name,'file')
        tic;
        load(name);
        details{j}.loaddur(1) = toc;
        details{j}.loadname = name;
        AllClusters{j} = Clusters;
        AllClusters{j}{1}.loadname = name;
        if exist('FullVData','var')
            AllFullVData{j} = FullVData;
        end
        for k = 1:length(AutoClusters)
            if k > length(Clusters) || isempty(Clusters{k})
                AllClusters{j}{k} = AutoClusters{k};
            end
        end
    elseif ~isempty(AutoClusters)
        AllClusters{j} = AutoClusters;
        fprintf('Can''t read %s\n',name);
    else
        fprintf('Can''t read %s or %s\n',name,aname);
        Clusters = {};
        return;
    end
    details{j}.loadtime = now;
    [ClusterDetails, CD] = LoadClusterDetails(name);
    details{j}.loaddur = cat(2,details{j}.loaddur, CD.loaddur);

    for k = 1:length(ClusterDetails)
        if isfield(ClusterDetails{k},'clst') && ...
                length(ClusterDetails{k}.clst) ~= AllClusters{j}{k}.nspks
            fprintf('Cluster %d E%d Details Clst (%d) does not match  Cluster nspks (%d)!!\n',k,e,length(ClusterDetails{k}.clst),AllClusters{j}{k}.nspks);
            AllClusters{j}{k}.needed = 1;
        elseif isfield(ClusterDetails{k},'ctime') && ...
                isfield(AllClusters{j}{k},'savetime') && ...
                ClusterDetails{k}.ctime > AllClusters{j}{k}.savetime(end) + 0.0001; %10 sec diff
            fprintf('Cluster %d Details newer than Cluster!!\n',k);
            AllClusters{j}{k}.needed = 1;
        end
        for n = 1:length(f)
            if isfield(ClusterDetails{k},f{n})
                AllClusters{j}{k}.(f{n}) = ClusterDetails{k}.(f{n});
            end
        end
        if isfield(ClusterDetails{k},'next')
        for c = 1:length(ClusterDetails{k}.next)
            %only load ClusterDetails if Cluster is still defined in next{c}
            if c <= length(AllClusters{j}{k}.next) && isfield(AllClusters{j}{k}.next{c},'space')
            for n = 1:length(f)
                if isfield(ClusterDetails{k}.next{c},f{n})
                    AllClusters{j}{k}.next{c}.(f{n}) = ClusterDetails{k}.next{c}.(f{n});
                end
            end
            end
        end
        end
 %some old files saved as row, not column
        if isfield(AllClusters{j},'clst') && size(AllClusters{j}(k).clst,1) == 1
            AllClusters{j}(k).clst = AllClusters{j}(k).clst';
        end
        if rawxy && isfield(ClusterDetails{k},'xy')
            xy = ClusterDetails{k}.xy;
            C = AllClusters{j}{k};
            if C.shape == 0
                AllClusters{j}{k}.xy = xy;
            else
                AllClusters{j}{k}.xy = xyrotate(xy(:,1),xy(:,2),-C.angle);
            end
        end
        if alltimes && isfield(AllClusters{j}{k},'t')
            AllClusters{j}{k}.times = AllClusters{j}{k}.t;
            AllClusters{j}{k} = rmfield(AllClusters{j}{k},'t');
        end
    end
    end
    for k = 1:length(AllClusters{j})
        C = FixCluster(AllClusters{j}{k});
        if isfield(C,'space')
            if ~isfield(C,'quick')
                C.quick = 0;
            end
            if ~isfield(C,'dropi')
                C.dropi = [0 0 0 0];
            end
            if ~isfield(C,'trigdt')
                C.trigdt = 0;
            end
            if ~isfield(C,'manual')
                C.manual = 0;
            end
            if ~isfield(C,'next')
                C.next = {};
            end
            if ~isfield(C,'clusterprog')
                C.clusterprog = '';
            end
            if ~isfield(C,'next')
                C.next = {};
            elseif isstruct(C.next)
                next = C.next;
                C = rmfield(C,'next');
                C.next{1} = next;
            end
            AllClusters{j}{k} = C;
            if isfield(C,'trigset')
                for c = 1:length(C.trigset{1}.next)
                    if c > length(C.next) || isempty(C.next{c}) || (isfield(C.next{c},'triggerset') &&C.next{c}.triggerset ==1)
                        if isfield(C.trigset{1}.next{c},'space')
                            C.next{c} = C.trigset{1}.next{c};
                            C.next{c}.trigset = 1;
                        end
                    end
                end
                AllClusters{j}{k} = rmfields(C, 'pcplot');
                if isfield(ClusterDetails{k},'trigset') && isfield(ClusterDetails{k}.trigset{1},'clst')
                    AllClusters{j}{k}.trigset{1}.clst = ClusterDetails{k}.trigset{1}.clst;
                    AllClusters{j}{k}.trigset{1}.t = ClusterDetails{k}.trigset{1}.t;
                    Cd = ClusterDetails{k}.trigset{1};
                    for c = 1:length(Cd.next)
                        if c > length(AllClusters{j}{k}.next) || ~isfield(AllClusters{j}{k}.next{c},'xy')
                            AllClusters{j}{k}.next{c}.xy = Cd.next{c}.xy;
                        end
                    end                    
                end
            end
                if isfield(AllClusters{j}{k},'times') && diff(size(AllClusters{j}{k}.times)) < 0
                    if isempty(AllClusters{j}{k}.times)
                        cprintf('red','%s P%d times is empty\n',name,k);
                    else
                        cprintf('red','%s P%d times is a row vector\n',name,k);
                        AllClusters{j}{k}.times = AllClusters{j}{k}.times';
                    end
                end
                if isfield(AllClusters{j}{k},'t') && diff(size(AllClusters{j}{k}.t)) < 0
                    cprintf('red','%s P%d t is a row vector\n',name,k);
                    AllClusters{j}{k}.t = AllClusters{j}{k}.t';
                end
            end
    details{j}.loaddur(end+1) = mytoc(ts);
    details{j}.exptno = eid;
end

if length(AllClusters) == 1
    AllClusters = FixCluster(AllClusters{1});
    if length(AllFullVData)
    AllFullVData = AllFullVData{1};
    end
    details = details{1};
end



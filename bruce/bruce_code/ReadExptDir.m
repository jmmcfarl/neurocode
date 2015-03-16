function [AllExpts, Ex] = ReadExptDir(name, varargin)
% [Expts, Idx] = ReadExptDir(name, varargin)
% Read all SPike2 .mat files in a directory and combines the Expts lists
% into one list
%
% Called by AplaySpkFile  if first argument is a directory

state.relist = 0;
state.resort = 0; %redo SortExpts, but not whole listing
state.online = 0;
state.quick = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'online',4)
        state.online=1;
    elseif strncmpi(varargin{j},'resort',4)
        state.resort=1;
    elseif strncmpi(varargin{j},'quick',4)
        state.quick=1;
    elseif strncmpi(varargin{j},'relist',4)
        state.relist = 1;
        state.resort = 1;
    end
    j=j+1;
end

expdates = [];
%first sort numerically by suffix number
%And find modification times of individual expts files
d = mydir([name '/*.mat']);
if isempty(d)
    AllExpts = {};
    Ex = {};
    cprintf('red','No Matlab files in %s\n',name);
    return;
end
% mnk =GetMonkeyName(name);
mnk = 'jbe';
suffixes = [];
expfiles = {};
for j = 1:length(d)
    if state.online
        if ~isempty(regexp(d(j).filename,'Expt[0-9]*.mat'))
            suffixes(j) = str2double(regexprep(d(j).filename,'Expt([0-9]*).mat','$1'));
        end
    elseif ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]*\.[0-9]+.mat']))
        suffixes(j) = str2double(regexprep(d(j).filename,'.*[\.,A-z,0-9]\.([0-9]+).mat','$1'));        
    elseif ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]+.mat']))
        suffixes(j) = str2double(regexprep(d(j).filename,'.*[\.,A-z]([0-9]+).mat','$1'));        
    elseif ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]*Expts.mat']))
        e = GetExptNumber(d(j).name);
        if e > 0
            expdates(e) = d(j).datenum;
            expfiles{e} = d(j).name;
        end
    end
end

if isempty(suffixes) %if no suffixes, must be in single files. 
    ne = 0;
    for j = 1:length(d)
        if ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]*idx.mat']))
            e = GetExptNumber(d(j).name);
            ne = ne+1;
            suffixes(j) = ne;
        end
    end
end

[a,b] =sort(suffixes);
sid = b(a> 0 & a < 2000); %get rid of Utah ns5->.mat files 

[a,b,c] = GetMonkeyName(name);
outname  = [name '/' a c 'Expts.mat'];
errfile = [name '/ClusterErrors.mat'];
    
if exist(outname) && state.resort == 0
        cerrexlist = [];
        xd = dir(outname);
    load(outname);
    if exist(errfile)
        load(errfile);
    else
        errorlist = [];
    end
    if isfield(errorlist,'ex')
        cerrexlist = [errorlist.ex];
    end
    if ~exist('Idx','var') %old file
        state.resort = 1;
    else
    reloaded = 0;
    for j = 1:length(Expts)
        if ~isfield(Expts{j}.Header,'suffix')
            Expts{j}.Header.suffix = GetExptNumber(Expts{j});
        end
        e = Expts{j}.Header.suffix;
        errid = find(ismember(cerrexlist,e));
        if ~isempty(errid)
            Expts{j}.clustererrs = errorlist(errid);
        end
        if e > 0 && length(expdates) >= e && expdates(e) > xd.datenum
            fprintf('Reloading %d\n',e);
            E = load(expfiles{e});
            if ~isempty(E.Expts)
                Expts{j} = E.Expts{1};
            end
            Idx{j} = E.Tidx;
            reloaded = reloaded+1;
        end
        if isfield(Expts{j}.Header,'loadname')
            a= fileparts(Expts{j}.Header.loadname);
           loadname = strrep(Expts{j}.Header.loadname,a,name);
           Expts{j}.Header.loadname = loadname;
        end
        if ~isfield(Expts{j},'errs') && length(Idx) == length(Expts) && isfield(Idx{j},'errs')
            Expts{j}.errs = Idx{j}.errs;
        end
    end
    if reloaded
        fprintf('Saving Expts with reloaded suffixes\n');
        save(outname,'Expts','Idx');
    end
    AllExpts = SetTrialOffsets(Expts);
    Ex = Idx;
    return;
    end
end

AllExpts = {};
AllErrs = {};
combineexpts = 0;

if isempty(sid) %No files with Expt data
    return;
end
nex = 1;
for j = 1:length(sid)
    fprintf('Reading %s\n',d(sid(j)).name);
    [Ex{nex}, Expts] = APlaySpkFile(d(sid(j)).name,'nospikes','noerrs', varargin{:});
    if length(Expts) > 1 && length(sid) > 2 %Usually this is an error, and the first expt is bad.  
%Could add a test later to try and combine these if possible        
        fprintf('%d Expts in %s\n',length(Expts),d(sid(j)).name);
        if combineexpts
            AllExpts = {AllExpts{:} CombineExpts(Expts)};
        else
            AllExpts = {AllExpts{:} Expts{end}};
        end
    elseif ~isempty(Expts)
        AllExpts = {AllExpts{:} Expts{:}};
    end
    if isfield(Ex{nex},'errs') && ~isempty(Ex{nex}.errs)
        cprintf('blue','Errors for %s\n',d(sid(j)).name);
        for k = 1:length(Ex{nex}.errs)
            cprintf('red','%s\n',Ex{nex}.errs{k});
        end
        if isempty(Expts)
            AllErrs = {AllErrs{:} Ex{nex}.errs};
        elseif isfield(AllExpts{end},'errs')
            if isfield(AllExpts{end}.errs,'msg')
                for k = 1:length(Ex{nex}.errs)
                    AllExpts{end}.errs(end+1).msg = Ex{nex}.errs{k};
                    AllExpts{end}.errs(end+1).t = 0;
                end
            else
                AllExpts{end}.errs = {AllExpts{end}.errs{:} Ex{nex}.errs{:}};
            end
        else
            AllExpts{end}.errs = Ex{nex}.errs;
        end
    end
    nex = length(AllExpts)+1; %so that it lines up with Expts,
end
AllExpts = SetTrialOffsets(AllExpts);


for j = 1:length(AllExpts)
    e = GetExptNumber(AllExpts{j});
    if isfield(AllExpts{j},'Header') && e > 1
        AllExpts{j}.Header.exptno = e;
    end
end

for j = 1:length(Ex)
    if isfield(Ex{j},'errs') && ~isempty(Ex{j}.errs)
        cprintf('blue','Errors for %s\n',d(sid(j)).name);
        for k = 1:length(Ex{j}.errs)
            cprintf('red','%s\n',Ex{j}.errs{k});
        end
    end
    e = GetExptNumber(AllExpts{j});
    if e > 0 && e <= length(expfiles) && exist(expfiles{e})
        E = load(expfiles{e});
        if isfield(E,'Expts') && ~isempty(E.Expts)
            exptno = GetExptNumber(E.Expts{1});
            if exptno > 1
                E.Expts{1}.Header.exptno = exptno;
            end
            E.Expts{1}.Header.trialoffset = AllExpts{j}.Header.trialoffset;
            save(expfiles{e},'-struct','E');
        end
        
    end
end
[a,b,c] = GetMonkeyName(name);
outname  = [name '/' a c 'Expts.mat'];
Expts = AllExpts;
Idx = Ex;

if ~exist(fileparts(outname))
    cprintf('error','Cannot save %s. Folder does not exist\n',outname);
elseif state.quick == 0 || ~exist(outname) %don't write out Expts if didn't load everything
    save(outname, 'Expts','Idx','AllErrs');
end

if nargout > 1 %? combine Ex{}
    
end
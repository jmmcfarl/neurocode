function n = GetExptNumber(name, varargin)
%n = GetExptNumber extracts the number of an experiment from a filename
%
n = 0;


if isfield(name,'Header')
    E = name;
    if isfield(E.Header,'suffix') && E.Header.suffix > 0
        n = E.Header.suffix;
    elseif ~isfield(E.Header,'bysuffix')
        if isfield(E.Header,'exptno')
            n = E.Header.exptno;
        end            
    elseif isfield(E.Header,'exptno') && E.Header.bysuffix == 0 
        n = E.Header.exptno;
    elseif isfield(E.Header,'Name')
        n = GetExptNumber(E.Header.Name);
    else
        n = NaN;
    end
    return;
end

if iscell(name)
    if isempty(name)
        n  =0;
        return;
    end
    for j = 1:length(name)
        n(j) = GetExptNumber(name{j}, varargin{:});
    end
%If its a set of clusters, just return number for the set    
    x = unique(n);
    if sum(x(x>0)) == 1
        n = x(x>0);
    end
    if sum(n) == 0 && isfield(name{1},'Header') %An old Expts without exptno
        n = 1:length(name);
    end
    return;
end

if isstruct(name)
    if isfield(name,'exptno')
        n = name.exptno;
        return;
    elseif isfield(name,'suffix')
        n = name.suffix;
        return;
    end
end
offset = 4;
if isempty(name)
    n = 0;
return;
end
id = regexp(name, 'Expt[0-9,a]*[ACF][lu].*.mat');
if isempty(id)
    id = regexp(name, 'Expt[0-9]*.p[0-9]*FullV.*.mat');
    
end    
if isempty(id)
    id = regexp(name, '\.[0-9]*Expts.*.mat');
    offset = 1;
end    
xid = regexp(name, 'Expt[0-9]*a[ACF][lu].*.mat');
if length(id) == 1
    n = sscanf(name(id+offset:end),'%d');
    if length(xid)
        n = n+0.1;
    end
elseif regexp(name,'\.p[0-9]+t[0-9]+')
    xid = regexp(name,'\.p[0-9]+t[0-9]+');
    a = sscanf(name(xid(1):end),'.p%dt%d');
    n = a(2);
else
    id = regexp(name, '\.[0-9]*.mat');
    if isempty(id)
        id = regexp(name, '\.[0-9]*.lfp.mat');
    end
    if ~isempty(id)
        n = sscanf(name(id+1:end),'%d');
    end
end

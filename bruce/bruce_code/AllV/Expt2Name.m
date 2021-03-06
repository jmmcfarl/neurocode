function [expname, exptypename, suff, stimname]  = Expt2Name(Expt, varargin)
% take an expt and return a name identifying the
% exptype, of the form rds.dxXce  
SpkDefs;
suff = [];
addsuff = 0;
crtrial = 0;
addprobe = 0;
nt = length(Expt.Trials);

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'addsuff',6)
        addsuff = 1;
    elseif strncmpi(varargin{j},'addprobe',6)
        addprobe = 1;
    end
    j = j+1;
end


if isfield(Expt.Stimvals,'st')
    if ischar(Expt.Stimvals.st)
    stimname = Expt.Stimvals.st;
    else
    stimname = stimnames{Expt.Stimvals.st+1};
    end
else
    stimname = '';
end

    if strmatch(Expt.Stimvals.e2, 'e0')
        exptypename = Expt.Stimvals.et;
        if strmatch(Expt.Stimvals.e3, 'e0')
            exptypename = Expt.Stimvals.et;
        else
            exptypename = [Expt.Stimvals.et 'X' Expt.Stimvals.e3];
        end
    else
        if strmatch('or',Expt.Stimvals.et) & Expt.Stimvals.ei > 100 & strcmp(Expt.Stimvals.e2,'me')
            exptypename = 'dirXme';
        else
        if strmatch(Expt.Stimvals.e3, 'e0')
            exptypename = [Expt.Stimvals.et 'X' Expt.Stimvals.e2];
        else
            exptypename = [Expt.Stimvals.et 'X' Expt.Stimvals.e2 'X' Expt.Stimvals.e3];
        end
        end
    end
    
    if isfield(Expt,'Header')
     
    if isfield(Expt.Header,'psych') && Expt.Header.psych
        exptypename = [exptypename 'P'];
    end
        if isfield(Expt.Header,'rc') && Expt.Header.rc
        suff = 'RC';
    end
    end
    if strcmp(Expt.Stimvals.et,'or') && Expt.Stimvals.ei > 100
        exptypename = strrep(exptypename, 'or','dir');
    end
    
    if strcmp(Expt.Stimvals.et,'tf') & crtrial > nt/2 %Counterphase
        exptypename = ['C' exptypename];
    end
    if isfield(Expt.Stimvals,'rb') & Expt.Stimvals.rb ~= 0
        exptypename = [exptypename 'RB'];
    end

    if strmatch(exptypename,{'dxXId' 'dxXIdP' 'orXId' 'dirXId'},'exact')
        jx = GetEval(Expt,'jx');
        if jx > 0
         exptypename = [exptypename 'D'];
        end
        if strmatch(Expt.Stimvals.Bs,'cylinder')
         exptypename = [exptypename 'B'];
        end
    end
    if strmatch(exptypename,{'OpXdx' 'PpXdx'},'exact')
        Bs = GetEval(Expt,'Bs');
        bh = GetEval(Expt,'bh');
        st = GetEval(Expt,'st');
        sz = GetEval(Expt,'sz');
        if Bs == st && bh >= sz
            exptypename = [exptypename];
        else
            exptypename = [exptypename 'noback'];
        end
            
        
    end
    if strncmp(exptypename,'dxXce',5)
        if isfield(Expt.Stimvals,'n2') && Expt.Stimvals.n2 > 2
            exptypename = strrep(exptypename,'dxXce','dxXces');
        elseif Expt.Stimvals.i2 < 2 %not +- 1
            exptypename = strrep(exptypename,'dxXce','dxXces');
        end
    end

    expname = [stimname '.' exptypename];
if addsuff
    expname = [expname suff];
end

if addprobe
    if isfield(Expt.Header,'cellnumber') && Expt.Header.cellnumber > 0
        expname = [expname '.cell' num2str(Expt.Header.cellnumber)];
    elseif isfield(Expt.Header,'probe')
        expname = [expname '.p' num2str(Expt.Header.probe)];
    end
end

    

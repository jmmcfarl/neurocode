function [fo] = setGNM_NLBFs(fo,Stim,p,n_bfs,space_type)
if nargin < 3
    p=.05;
end
if nargin < 4
    n_bfs = [];
end
if nargin < 5
    space_type = 1;
end

g = Stim*get_k_mat(fo);
nmods = length(fo.mods); %number of modules
for imod=1:nmods; 	%for each module
    
    if ~isempty(n_bfs)
        n_dxs = n_bfs;
    else
        %update the tent basis X-axis
        n_dxs = length(fo.mods(imod).nlx);
    end
    %         fo.mods(imod).nlx = prctile(g(:,imod),linspace(0.05,99.95,n_dxs)); %equipopulated
    if space_type==1
        fo.mods(imod).nlx = linspace(prctile(g(:,imod),p),prctile(g(:,imod),100-p),n_dxs); %equispacing
    elseif space_type==0
        fo.mods(imod).nlx = prctile(g(:,imod),linspace(p,100-p,n_dxs)); %equipopulated
    else
        error('Not supported')
    end
    %set nearest tent basis to 0
    [~,nearest] = min(abs(fo.mods(imod).nlx));
    fo.mods(imod).nlx(nearest) = 0;
    if strcmp(fo.mods(imod).nltype,'threshlin')
        fo.mods(imod).nly = fo.mods(imod).nlx;
        fo.mods(imod).nly(1:nearest) = 0;
    elseif strcmp(fo.mods(imod).nltype,'quad')
        fo.mods(imod).nly = fo.mods(imod).nlx.^2;
    elseif strcmp(fo.mods(imod).nltype,'lin')
        fo.mods(imod).nly = fo.mods(imod).nlx;
    end
     
end;


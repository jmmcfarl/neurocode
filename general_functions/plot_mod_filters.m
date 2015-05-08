function [h,fig_dims,stim_filters,mod_signs] = plot_mod_filters(model,sp_dx,dt)

if nargin < 2 || isempty(sp_dx);
    sp_dx = 1;
end
if nargin < 3 || isempty(dt)
    dt = 1;
end

stim_mods = find([model.mods(:).Xtarget] == 1);
flen = model.stim_params(1).stim_dims(1);
nPix = model.stim_params(1).stim_dims(2);
mod_signs = [model.mods(stim_mods).sign];
NLtypes = {model.mods(stim_mods).NLtype};
stim_filters = reshape([model.mods(stim_mods).filtK],flen,nPix,length(stim_mods));

n_rows = floor(sqrt(length(stim_mods)));
n_cols = ceil(length(stim_mods)/n_rows);
fig_dims = [n_rows n_cols];

xax = (1:nPix)*sp_dx;
xax = xax - nanmean(xax);
tax = (0:(flen-1))*dt + dt/2;
tax = tax*1e3;

h = figure();
for ii = 1:length(stim_mods)
    subplot(n_rows,n_cols,(ii-1)+1);
    imagesc(xax,tax,squeeze(stim_filters(:,:,ii)));
    set(gca,'ydir','normal');
    colormap(gray);
    cam = max(abs(caxis()));
    caxis([-cam cam]);
    if mod_signs(ii)==1
        stype = 'E';
        sc = 'r';
    else
        stype = 'I';
        sc = 'b';
    end
    title(sprintf('%s %s-Filt',NLtypes{ii},stype),'Color',sc);
    xlabel('Rel position (deg)');
    ylabel('Time lag (ms)');
end
figufy(h);
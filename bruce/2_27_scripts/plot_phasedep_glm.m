function [h] = plot_phasedep_glm(phase_filts,phase_offset,phase_kLin,Xstim,phase_sig,init_mod,is_vis)

n_phase_bins = length(phase_offset);
phase_bin_edges = linspace(-pi,pi,n_phase_bins + 1);
phase_bin_edges(end) = [];
PBx = phase_bin_edges;

stim_params = init_mod.stim_params;
n_bars = stim_params.stim_dims(2);
nLags = stim_params.stim_dims(1);
kLen = nLags*n_bars;
filtLen = nLags*n_bars;

%% CREATE TENT BASIS FUNCTION OUTPUTs
nPBs = length(PBx);
XPBx = zeros(length(phase_sig),nPBs);
dp = median(diff(PBx));
for j = 2:nPBs-1
    cur_set = find(phase_sig > PBx(j-1) & phase_sig < PBx(j)); %left side
    XPBx(cur_set,j) = XPBx(cur_set,j) + (1 - (PBx(j) - phase_sig(cur_set))/dp);
    cur_set = find(phase_sig >= PBx(j) & phase_sig < PBx(j+1)); %right side
    XPBx(cur_set,j) = XPBx(cur_set,j) + (1 - (phase_sig(cur_set)-PBx(j))/dp);
end
%first TB
cur_set = find(phase_sig >= PBx(1) & phase_sig < PBx(2)); %right side
XPBx(cur_set,1) = XPBx(cur_set,1) + (1 - (phase_sig(cur_set)-PBx(1))/dp);
cur_set = find(phase_sig >= PBx(end)); %left side
XPBx(cur_set,1) = XPBx(cur_set,1) + (1 - (pi-phase_sig(cur_set))/dp);
%last TB
cur_set = find(phase_sig > PBx(end-1) & phase_sig < PBx(end)); %left side
XPBx(cur_set,end) = XPBx(cur_set,end) + (1 - (PBx(end)-phase_sig(cur_set))/dp);
cur_set = find(phase_sig >= PBx(end)); %right side
XPBx(cur_set,end) = XPBx(cur_set,end) + (1 - (phase_sig(cur_set)-PBx(end))/dp);

%%
kmat = reshape(phase_filts,filtLen,nPBs);
filt_outs = Xstim*kmat;
all_filt_outs = filt_outs(:);
kax = linspace(prctile(all_filt_outs,0.2),prctile(all_filt_outs,99.8),200);
for pb = 1:n_phase_bins    
    f(pb,:) = ksdensity(filt_outs(:,pb),kax);   
    yy(pb,:) = log(1+exp(kax + mean(phase_kLin) + phase_offset(pb)));
    cumf = cumtrapz(kax,f(pb,:));
    hp = find(cumf > 0.998,1,'first');
    if ~isempty(hp)
        yy(pb,hp:end) = nan;
    end
end

%%
if is_vis == 0
    h = figure('visible','off');
else
h = figure;
end
cx = max(abs(phase_filts));
cax = [-cx cx];
fyr = [0 max(f(:))];
yyr = [0 nanmax(yy(:))];
for pb = 1:n_phase_bins
    subplot(3,n_phase_bins,pb)
    imagesc(flipud(reshape(kmat(:,pb),nLags,n_bars)));colormap(jet);
    caxis(cax);
    subplot(3,n_phase_bins,n_phase_bins + pb)
    [ax,h1,h2] = plotyy(kax,f(pb,:),kax,yy(pb,:));
    xlim(ax(1),kax([1 end]));
    xlim(ax(2),kax([1 end]));
    ylim(ax(1),fyr); ylim(ax(2),yyr);
set(ax(1),'ytick',[]); set(ax(2),'ytick',[]);
end
subplot(3,n_phase_bins,2*n_phase_bins + 1)
imagesc(kax,1:n_phase_bins,f);
subplot(3,n_phase_bins,2*n_phase_bins + 2)
imagesc(kax,1:n_phase_bins,yy);

init_k = init_mod.mods(1).filtK;
init_g = Xstim*init_k;
subplot(3,n_phase_bins,2*n_phase_bins + 4)
fi = ksdensity(init_g,kax);
yi = log(1+exp(kax + mean(init_mod.kLin) + init_mod.spk_NL_params(1)));
[ax,h1,h2] = plotyy(kax,fi,kax,yi);
xlim(ax(1),kax([1 end]));
xlim(ax(2),kax([1 end]));
ylim(ax(1),fyr); ylim(ax(2),yyr);
set(ax(1),'ytick',[]); set(ax(2),'ytick',[]);

subplot(3,n_phase_bins,2*n_phase_bins + 5);
imagesc(reshape(init_mod.mods(1).filtK,nLags,n_bars));
caxis(cax);

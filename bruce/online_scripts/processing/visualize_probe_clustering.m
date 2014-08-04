function [] = visualize_probe_clustering(probe_nums,expt_nums,stored_spkxy,autoclust,use_plots,dsize)

close all

if length(probe_nums) == 1
    probe_nums = probe_nums*ones(size(expt_nums));
end


n_expts = length(expt_nums);

n_cols = ceil(sqrt(n_expts));
n_rows = ceil(n_expts/n_cols);

plot_dens = use_plots(1);
plot_wvfrm = use_plots(2);
plot_scat = use_plots(3);

if plot_dens
    dens_fig = figure;
end
if plot_wvfrm
    wvfrm_fig = figure;
end
if plot_scat
    scat_fig = figure;
end

scrsz = get(0,'ScreenSize');

for ee = 1:n_expts
    if plot_dens
        figure(dens_fig);
        subplot(n_rows,n_cols,ee);
        if ~isempty(stored_spkxy{expt_nums(ee),probe_nums(ee)})
        [handles, details] = DensityPlot_jmm(stored_spkxy{expt_nums(ee),probe_nums(ee)}(:,1),stored_spkxy{expt_nums(ee),probe_nums(ee)}(:,2),'sqrtsc');
        hold on
        if autoclust(expt_nums(ee),probe_nums(ee)).man_code ~=4
        z = pdf(autoclust(expt_nums(ee),probe_nums(ee)).gmm,[details.x(:) details.y(:)]);
        z = sqrt(z);
        contour(details.x, details.y, reshape(z,size(details.x)),10,'w','linewidth',0.5);
        end
        set(gca,'xticklabel',[],'yticklabel',[],'ydir','normal');
        end
        title(sprintf('Expt %d',expt_nums(ee)));
    end
    if plot_wvfrm
        figure(wvfrm_fig);
        subplot(n_rows,n_cols,ee);
        if ~isempty(stored_spkxy{expt_nums(ee),probe_nums(ee)})
        hold on
        tax = (1:size(autoclust(expt_nums(ee),probe_nums(ee)).avg_wvfrm,1))/3e4;
      shadedErrorBar(tax,autoclust(expt_nums(ee),probe_nums(ee)).avg_wvfrm(:,2),autoclust(expt_nums(ee),probe_nums(ee)).std_wvfrm(:,2),{'color','k'});
        shadedErrorBar(tax,autoclust(expt_nums(ee),probe_nums(ee)).avg_wvfrm(:,1),autoclust(expt_nums(ee),probe_nums(ee)).std_wvfrm(:,1),{'color','r'});
        xlim([0 1.2e-3])
        set(gca,'xticklabel',[],'yticklabel',[]);
        end
        title(sprintf('Expt %d',expt_nums(ee)));
    end
    if plot_scat
        figure(scat_fig);
        subplot(n_rows,n_cols,ee);
        if ~isempty(stored_spkxy{expt_nums(ee),probe_nums(ee)})
        hold on
        plot(stored_spkxy{expt_nums(ee),probe_nums(ee)}(autoclust(expt_nums(ee),probe_nums(ee)).idx == 2,1),...
            stored_spkxy{expt_nums(ee),probe_nums(ee)}(autoclust(expt_nums(ee),probe_nums(ee)).idx == 2,2),'k.','markersize',dsize);
        plot(stored_spkxy{expt_nums(ee),probe_nums(ee)}(autoclust(expt_nums(ee),probe_nums(ee)).idx == 1,1),...
            stored_spkxy{expt_nums(ee),probe_nums(ee)}(autoclust(expt_nums(ee),probe_nums(ee)).idx == 1,2),'r.','markersize',dsize);
        set(gca,'xticklabel',[],'yticklabel',[]);
        title(sprintf('Expt %d',expt_nums(ee)));
        end
        axis tight
    end
end
if plot_dens && n_expts > 1
    set(dens_fig,'position',[1 1 3*scrsz(3)/4 scrsz(4)-100]);
end
if plot_scat && n_expts > 1
    set(scat_fig,'position',[100 1 3*scrsz(3)/4 scrsz(4)-100]);
end
if plot_wvfrm && n_expts > 1
    set(wvfrm_fig,'position',[200 1 3*scrsz(3)/4 scrsz(4)-100]);
end
if n_expts > 1
    dprimes = [autoclust(:,probe_nums(1)).dprime];
    mahal_ds = [autoclust(:,probe_nums(1)).mahal_d];
    dist_fig = figure;
    subplot(2,1,1)
    plot(dprimes,'o-');
    title('Dprime');
    subplot(2,1,2)
    plot(mahal_ds,'o-');
    title('Mahal dist');
    set(dist_fig,'position',[scrsz(3)-500 scrsz(4)/2 500 scrsz(4)/2-100]);
end
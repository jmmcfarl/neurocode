function [] = NIMdisplay_model(nim,Xstim,XLin,Robs)
%
% NIMdisplay_model(nim,<Xstim>,<XLin>,<Robs>)
%
% Creates a display of the elements of a given NIM
% INPUTS:
%     nim: model structure
%     <Xstim>: provide the stimulus matrix if you want to display the distributions of generating signals
%     <XLin>: If there are additional linear covariates in the model, provide them too (for displaying generating distributions)
%     <Robs>: If there are spike history terms (and you want to display generating distributions)

%%
if nargin < 2
    Xstim = [];
end
if nargin < 3
    XLin = [];
end
if nargin < 4
    Robs = [];
end

%%
nmods = length(nim.mods);
spkhstlen = nim.spk_hist.spkhstlen;
n_hist_bins = 500; %internal parameter determining histogram resolution

if ~isempty(Xstim)
    %compute generating signals (and internal generating signals)
    if ~isempty(Robs)
        [~, ~, ~, G, gint] = NIMmodel_eval(nim,Robs,Xstim,XLin);
    elseif spkhstlen == 0
        [~, ~, ~, G, gint] = NIMmodel_eval(nim,[],Xstim,XLin);
    else
        % error('Need to provide Robs to compute spike history filter output');
				disp('Need to provide Robs to compute spike history filter output. Spk-NL distribution will be off.');
				Robs = zeros(size(Xstim,1),1);
        [~, ~, ~, G, gint] = NIMmodel_eval(nim,Robs,Xstim,XLin);
		end
		% Remove spike history offset from G (since officially part of spiking nonlinearity
		theta = nim.spk_NL_params(1);
		G = G - theta;
else
    G = []; gint = [];
end

%%
% if nim.stim_params.stim_dims(3) > 1 %if more than 1 spatial dim
%     error('Plotting for 2 spatial dims not supported yet!');
% end

%% PLOT SPIKE HISTORY TERM
if spkhstlen > 0
    fig_spkhist = figure();
    subplot(2,1,1)
		stairs(nim.spk_hist.bin_edges(1:end-1)*nim.stim_params.dt,nim.spk_hist.coefs);
    xlim(nim.spk_hist.bin_edges([1 end])*nim.stim_params.dt)
    xl = xlim();
    line(xl,[0 0],'color','k','linestyle','--');
    xlabel('Time lag (bins)');
    ylabel('Spike history filter')
    title('Spike history term','fontsize',14)
    
    subplot(2,1,2)
    stairs(nim.spk_hist.bin_edges(1:end-1)*nim.stim_params.dt,nim.spk_hist.coefs);
    xlim(nim.spk_hist.bin_edges([1 end-1])*nim.stim_params.dt)
    set(gca,'xscale','log')
    xl = xlim();
    line(xl,[0 0],'color','k','linestyle','--');
    xlabel('Time lag');
    ylabel('Spike history filter')
    title('Spk Hist Log-axis','fontsize',14)
    
    %set horizontal position of figure
    cp = get(fig_spkhist,'Position');
    cp(1) = 600;
    set(fig_spkhist,'Position',cp);
end

%% PLOT LINEAR FILTER
if nim.stim_params.lin_dims > 0
    fig_linfilt = figure();
    plot(nim.kLin,'o-');
    xlabel('COVARIATE INDEX');
    ylabel('LINEAR FILTER')
    title('ADDITIONAL COVARIATE FILTER','fontsize',14)
    
    %set horizontal position of figure
    cp = get(fig_linfilt,'Position');
    cp(1) = 400;
    set(fig_linfilt,'Position',cp);
end

%% CREATE FIGURE SHOWING INDIVIDUAL SUBUNITS
fig_full = figure();
if nim.stim_params.stim_dims(3) > 1
    n_columns = nim.stim_params.stim_dims(1) + 1;
    n_rows = nmods;
else
    n_columns = max(round(sqrt(nmods/2)),1);
    n_rows = ceil(nmods/n_columns);
end
nLags = nim.stim_params.stim_dims(1);
dt = nim.stim_params.dt;
nPix = squeeze(nim.stim_params.stim_dims(2:end));

%create filter time lag axis
if isempty(nim.stim_params.tent_spacing)
    tax = (0:(nLags-1))*dt;
else
    tax = (0:nim.stim_params.tent_spacing:(nLags-1)*nim.stim_params.tent_spacing)*dt;
end
tax = tax * 1000; % put in units of ms

xr = tax([1 end]);
for imod=1:nmods
    thismod = nim.mods(imod);
    
    if nim.stim_params.stim_dims(3) == 1
        %PLOT FILTER
        subplot(n_rows,2*n_columns,(imod-1)*2+1);
        if nPix == 1 %if temporal-only stim
            plot(tax,thismod.filtK,'.-');
            line(xr,[0 0],'color','k','linestyle','--');
            xlim(xr);
            xlabel('Time lag')
            ylabel('Filter coef');
        elseif nPix(2) == 1
            imagesc(1:nPix(1),tax,reshape(thismod.filtK,nLags,nPix(1)));
            cl = max(abs(thismod.filtK));
            caxis([-cl cl]);
            %colormap(jet);
						colormap(gray);
            set(gca,'ydir','normal');
            xlabel('Pixels')
            ylabel('Time lags');
        end
        if strcmp(thismod.NLtype,'lin')
            title('Linear stimulus filter','fontsize',14)
        elseif thismod.sign == 1
            title('Excitatory stimulus filter','fontsize',14);
        elseif thismod.sign == -1
            title('Suppressive stimulus filter','fontsize',14);
        end
    else
        for tt = 1:nim.stim_params.stim_dims(1)
            subplot(n_rows,n_columns,(imod-1)*n_columns + tt);
            cur_fdims = tt - 1 + (1:nim.stim_params.stim_dims(1):prod(nim.stim_params.stim_dims));
            imagesc(1:nPix(1),1:nPix(2),reshape(thismod.filtK(cur_fdims),nim.stim_params.stim_dims(2:end)));
						colormap(gray)
        if strcmp(thismod.NLtype,'lin')
            title(sprintf('Lin-input Lag %d',tt-1),'fontsize',10);
        elseif thismod.sign == 1
            title(sprintf('E-Input Lag %d',tt-1),'fontsize',10);
        elseif thismod.sign == -1
            title(sprintf('S-Input Lag %d',tt-1),'fontsize',10);
        end
       end
    end
    
    %PLOT UPSTREAM NL
    if nim.stim_params.stim_dims(3) == 1
        subplot(n_rows,2*n_columns,(imod-1)*2+2);
    else
        subplot(n_rows,n_columns,(imod)*n_columns);
    end
    if ~isempty(gint) %if computing distribution of filtered stim
        [gendist_y,gendist_x] = hist(gint(:,imod),n_hist_bins);
				
				% Sometimes the gendistribution has a lot of zeros (dont want to screw up plot)
				[a b] = sort(gendist_y);
				if a(end) > a(end-1)*1.5
					gendist_y(b(end)) = gendist_y(b(end-1))*1.5;
				end
    else
        gendist_x = linspace(-3,3,n_hist_bins); %otherwise, just pick an arbitrary x-axis to plot the NL
    end
    
    if strcmp(thismod.NLtype,'nonpar')
        cur_modx = thismod.NLx; cur_mody = thismod.NLy;
    elseif strcmp(thismod.NLtype,'lin')
        cur_modx = gendist_x; cur_mody = cur_modx;
    elseif strcmp(thismod.NLtype,'quad')
        cur_modx = gendist_x;
        cur_mody = cur_modx.^2;
    elseif strcmp(thismod.NLtype,'threshlin')
        cur_modx = gendist_x;
        cur_mody = cur_modx;
        cur_mody(cur_mody < 0) = 0;
    end
    cur_xrange = cur_modx([1 end]);
    
    if ~isempty(gint)
        [ax,h1,h2] = plotyy(cur_modx,cur_mody,gendist_x,gendist_y);
        if strcmp(thismod.NLtype,'nonpar')
            set(h1,'Marker','o');
        end
        %         set(h2,'Color','k')
        set(h1,'linewidth',1)
        xlim(ax(1),cur_xrange)
        xlim(ax(2),cur_xrange);
        ylim(ax(1),[min(cur_mody) max(cur_mody)]);
        set(ax(2),'ytick',[])
        yl = ylim();
        line([0 0],yl,'color','k','linestyle','--');
        ylabel(ax(1),'Subunit output','fontsize',12);
        ylabel(ax(2),'Probability','fontsize',12)
    else
        h = plot(cur_modx,cur_mody,'linewidth',1);
        if strcmp(thismod.NLtype,'nonpar')
            set(h,'Marker','o');
        end
        xlim(cur_xrange)
        ylim([min(cur_mody) max(cur_mody)]);
        ylabel('Subunit output','fontsize',12);
    end
    box off
    xlabel('Internal generating function')
    title('Upstream NL','fontsize',14)
end

%set horizontal position of fig
cp = get(fig_full,'Position');
cp(1) = 0;
set(fig_full,'Position',cp);

%% if temporal only stim, also create a plot directly comparing temporal filters
if nPix(1) == 1 && nmods > 1
    filt_compare = figure();hold on;
    if nmods == 2 %if 2 filters, use red-blue
        cmap = [0 0 1; 1 0 0];
    else
        cmap = colormap(jet(nmods));
    end
    for imod = 1:nmods
        plot(tax,nim.mods(imod).filtK,'.-','color',cmap(imod,:));
    end
    xl = xlim();
    line(xl,[0 0],'color','k','linestyle','--');
end

%% PLOT SPIKING NL FUNCTION
if ~isempty(G)
    
    fig_spk_nl = figure();
    n_bins = 1000;
    [Gdist_y,Gdist_x] = hist(G,n_hist_bins);
    
    %this is a hack to deal with cases where the threshold linear terms
    %create a min value of G
    if Gdist_y(1) > 2*Gdist_y(2)
        Gdist_y(1) = 1.5*Gdist_y(2);
    end
    
    cur_xrange = Gdist_x([1 end]);
    
    if strcmp(nim.spk_NL_type,'logexp')
        cur_y = nim.spk_NL_params(3)*log(1 + ...
            exp(nim.spk_NL_params(2)*Gdist_x + theta));
    elseif strcmp(nim.spk_NL_params,'exp')
        cur_y = exp(Gdist_x + theta);
    else
        error('Unsupported spk NL type');
    end
    cur_y = cur_y/nim.stim_params.dt; %convert to correct firing rate units
    [ax,h1,h2] = plotyy(Gdist_x,cur_y,Gdist_x,Gdist_y);
    %     set(h2,'Color','k')
    set(h1,'linewidth',1)
    yr = [min(cur_y) max(cur_y)];
    xlim(ax(1),cur_xrange)
    xlim(ax(2),cur_xrange);
    ylim(ax(1),yr);
    
    xlabel('Generating function')
    ylabel(ax(1),'Predicted firing rate','fontsize',14);
    ylabel(ax(2),'Probability','fontsize',14)
    set(ax(2),'ytick',[]);
    title('Spiking NL','fontsize',14)
    
    %set horizontal position
    cp = get(fig_spk_nl,'Position');
    cp(1) = 1100;
    set(fig_spk_nl,'Position',cp);
end


function [h] = plot_NMM_filters_1d(nim,pix_ax,lag_ax,stim_Xtarg)


if nargin < 4 || isempty(stim_Xtarg)
    stim_Xtarg = 1;
end
stim_params = nim.stim_params(stim_Xtarg);
nLags = stim_params.stim_dims(1);
dt = stim_params.dt;
nPix = stim_params.stim_dims(2);

if nargin < 2 || isempty(pix_ax)
    pix_ax = 1:nPix;
end
if nargin < 3 || isempty(lag_ax)
    lag_ax = ((1:nLags)*dt - dt/2)*1e3;
end
%%
Xtargs = [nim.mods(:).Xtargets];
stim_mods = find(Xtargs == stim_Xtarg);

%% CREATE FIGURE SHOWING INDIVIDUAL SUBUNITS

h = figure();
n_columns = max(round(sqrt(length(stim_mods)/2)),1);
n_rows = ceil(length(stim_mods)/n_columns);

for imod = 1:length(stim_mods)
    thismod = nim.mods(stim_mods(imod));
    
    %PLOT FILTER
    subplot(n_rows,2*n_columns,(imod-1)*2+1);
    if nPix == 1 %if temporal-only stim
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
end



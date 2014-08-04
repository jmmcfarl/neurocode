function [] = NMMCdisplay_model( nim, Xstim, Robs, disp_Xtargs )
%
% Usage: NMMCdisplay_model( nim, <Xstim>, <Robs>, <disp_Xtargs> )
%
% Creates a display of the elements of a given NIM
% INPUTS:
%     nim: model structure
%     <Xstim>: provide the stimulus matrix if you want to display the distributions of generating signals
%     <Robs>: If there are spike history terms (and you want to display generating distributions)
%     <disp_Xtargs>: vector of Xtargets for which filters will be displayed
%       (default all)
%%
if nargin < 2
	Xstim = [];
end
if nargin < 3
	Robs = [];
end

Xtargets = [nim.mods(:).Xtarget];
nTargets = length(unique(Xtargets));

if nargin < 4 
	disp_Xtargs = 1:nTargets;
end
    
%%
Nmods = length(nim.mods);
spkhstlen = nim.spk_hist.spkhstlen;
n_hist_bins = 500; %internal parameter determining histogram resolution

if ~isempty(Xstim)
	% Compute generating signals (and internal generating signals)
	if ~isempty(Robs)
    [~, ~, ~, G, gint] = NMMmodel_eval(nim,Robs,Xstim);
	elseif spkhstlen == 0
    [~, ~, ~, G, gint] = NMMmodel_eval(nim,[],Xstim);
  else
    % error('Need to provide Robs to compute spike history filter output');
    disp('Need to provide Robs to compute spike history filter output. Spk-NL distribution will be off.');
    Robs = zeros(size(Xstim,1),1);
    [~, ~, ~, G, gint] = NMMmodel_eval(nim,Robs,Xstim);
  end
	% Remove spike history offset from G (since officially part of spiking nonlinearity
	theta = nim.spk_NL_params(1);
	G = G - theta;
else
	G = []; gint = [];
end

%%
% Check if spike-history and/or spiking nonlinearity to plot
if (spkhstlen > 0) || ~isempty(G)
	Xtra_col = 1;
else
	Xtra_col = 0;
end

%% CREATE FIGURE SHOWING INDIVIDUAL SUBUNITS
plotconv = zeros(2,1);
for tt = Nmods
	if isfield(nim.mods(tt),'Xconv') && (length(nim.mods(tt).Xconv) > 1)
		plotconv(1) = 1;
	end
	if isfield(nim.mods(tt),'psc') && ~isempty(nim.mods(tt).psc)
		plotconv(2) = 1;
	end
end
NPL = 2 + sum(plotconv);

for tt = disp_Xtargs
	cur_mods = find(Xtargets == tt);
    
	cur_fig = figure();
	if nim.stim_params(tt).stim_dims(3) > 1
		n_columns = nim.stim_params(tt).stim_dims(1) + 1;
		n_rows = length(cur_mods);
	else
		n_columns = max(round(sqrt(length(cur_mods)/2)),1);
		n_rows = ceil(length(cur_mods)/n_columns);
	end
	nLags = nim.stim_params(tt).stim_dims(1);
	dt = nim.stim_params(tt).dt;
	nPix = squeeze(nim.stim_params(tt).stim_dims(2:end));
    
	% Create filter time lag axis
	if isempty(nim.stim_params(tt).tent_spacing)
		tax = (0:(nLags-1))*dt;
	else
		tax = (0:nim.stim_params(tt).tent_spacing:(nLags-1)*nim.stim_params(tt).tent_spacing)*dt;
	end
	tax = tax * 1000; % put in units of ms
    
  if (n_rows == 1) && (spkhstlen > 0) && ~isempty(G)
		Xtra_col = 2;
	end
	
	for imod = 1:length(cur_mods)
		thismod = nim.mods(cur_mods(imod));
        
		if nim.stim_params(tt).stim_dims(3) == 1
			% Plot filter
			subplot( n_rows, NPL*n_columns + Xtra_col, (imod-1)*(NPL+Xtra_col)+1 );
			if nPix < 4 % if dimensions easily visualized as temporal kernels
				if isfield(thismod, 'keat_basis')
					kblen = size(thismod.keat_basis,2);                    
					tax = (0:kblen-1)*dt*1000;
					plot(tax,thismod.filtK(:)'*thismod.keat_basis,'.-');                    
        else
					plot(tax,reshape(thismod.filtK,nim.stim_params(tt).stim_dims(1:2)),'.-');
				end
				xr = tax([1 end]);
				line(xr,[0 0],'color','k','linestyle','--');
        if diff(xr) > 0
          xlim(xr);
          xlabel('Time lag')
        end
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
			for jj = 1:nim.stim_params(tt).stim_dims(1)
				subplot( n_rows, n_columns+Xtra_col, (imod-1)*(n_columns+Xtra_col) + jj );
				cur_fdims = jj - 1 + (1:nim.stim_params(tt).stim_dims(1):prod(nim.stim_params(tt).stim_dims));
				imagesc(1:nPix(1),1:nPix(2),reshape(thismod.filtK(cur_fdims),nim.stim_params(tt).stim_dims(2:end)));
				colormap(gray)
                
				if strcmp(thismod.NLtype,'lin')
					title(sprintf('Lin-input Lag %d',jj-1),'fontsize',10);
				elseif thismod.sign == 1
					title(sprintf('E-Input Lag %d',jj-1),'fontsize',10);
				elseif thismod.sign == -1
					title(sprintf('S-Input Lag %d',jj-1),'fontsize',10);
				end
			end	
		end
        
		% Plot upstream nonlinearity
		if nim.stim_params(tt).stim_dims(3) == 1
			subplot( n_rows, NPL*n_columns+Xtra_col, (imod-1)*(NPL+Xtra_col)+2 );
		else
			subplot( n_rows, n_columns+Xtra_col, (imod)*(n_columns+Xtra_col) );
		end
		if ~isempty(gint) %if computing distribution of filtered stim
			[gendist_y,gendist_x] = hist(gint(:,imod),n_hist_bins);
            
			% Sometimes the gendistribution has a lot of zeros (dont want to screw up plot)
			[a b] = sort(gendist_y);
			if a(end) > a(end-1)*1.5
				gendist_y(b(end)) = gendist_y(b(end-1))*1.5;
			end
			
			% Make plot not reflect heavy tails
			cur_modx = gendist_x((gendist_x >= my_prctile(gint(:,imod),2)) & (gendist_x <= my_prctile(gint(:,imod),99)));
		else
			gendist_x = linspace(-3,3,n_hist_bins); %otherwise, just pick an arbitrary x-axis to plot the NL
			cur_modx = gendist_x;
		end
        
		if strcmp(thismod.NLtype,'nonpar')
			cur_modx = thismod.NLx; cur_mody = thismod.NLy;
		elseif strcmp(thismod.NLtype,'lin')
			%cur_modx = gendist_x; 
			cur_mody = cur_modx;
		elseif strcmp(thismod.NLtype,'quad')
			%cur_modx = gendist_x;
			cur_mody = cur_modx.^2;
		elseif strcmp(thismod.NLtype,'threshlin')
			%cur_modx = gendist_x;
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
			%ylabel(ax(2),'Probability','fontsize',12)
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
				
		% Plot Xconv
		if plotconv(1) 
			subplot( n_rows, NPL*n_columns+Xtra_col, (imod-1)*NPL+Xtra_col+3 );
			hold on
			plot( nim.mods(imod).Xshifts, nim.mods(imod).Xconv, 'b' )
			plot( nim.mods(imod).Xshifts, nim.mods(imod).Xconv, 'b.' )
			axis([nim.mods(imod).Xshifts([1 end]) 0 1])
			xlabel('Position')
			title('Spatial Conv.','fontsize',14)
		end
		if plotconv(2) && ~isempty(nim.mods(imod).psc)
			subplot( n_rows, NPL*n_columns+Xtra_col, (imod-1)*(NPL+Xtra_col)+NPL );
			hold on
			plot( (0:(length(nim.mods(imod).psc)-1))*nim.mods(imod).psc_spacing, nim.mods(imod).psc, 'b' )
			plot( (0:(length(nim.mods(imod).psc)-1))*nim.mods(imod).psc_spacing, nim.mods(imod).psc, 'b.' )
			axis tight
			xlabel('Time lag')
			title('PSC term','fontsize',14)		
		end
	end
    
end

%% PLOT SPIKE HISTORY TERM
if spkhstlen > 0
	%fig_spkhist = figure();
	if n_rows > 1
		subplot( n_rows, NPL*n_columns+Xtra_col, NPL*n_columns+Xtra_col )
	else
		if (spkhstlen > 0) && ~isempty(G)
			subplot( n_rows, NPL*n_columns+Xtra_col, NPL*n_columns+Xtra_col-1 )
		else
			subplot( n_rows, NPL*n_columns+Xtra_col, NPL*n_columns+Xtra_col )
		end
	end
	stairs(nim.spk_hist.bin_edges(1:end-1)*nim.stim_params(1).dt,nim.spk_hist.coefs);
	xlim(nim.spk_hist.bin_edges([1 end])*nim.stim_params(1).dt)
	xl = xlim();
	line(xl,[0 0],'color','k','linestyle','--');
	xlabel('Time lag (bins)');
	ylabel('Spike history filter')
	title('Spike history term','fontsize',12)
    
  %  subplot(2,1,2)
  %  stairs(nim.spk_hist.bin_edges(1:end-1)*nim.stim_params(1).dt,nim.spk_hist.coefs);
  %  xlim(nim.spk_hist.bin_edges([1 end-1])*nim.stim_params(1).dt)
  %  set(gca,'xscale','log')
  %  xl = xlim();
  %  line(xl,[0 0],'color','k','linestyle','--');
  %  xlabel('Time lag');
  %  ylabel('Spike history filter')
  %  title('Spk Hist Log-axis','fontsize',14)
    
	%set horizontal position of figure
	%cp = get(fig_spkhist,'Position');
	%cp(1) = 600;
	%set(fig_spkhist,'Position',cp);
end

%% PLOT SPIKING NL FUNCTION
if ~isempty(G)
    
	%fig_spk_nl = figure();
	if n_rows > 1
		subplot( n_rows, NPL*n_columns+Xtra_col, 2*(NPL*n_columns+Xtra_col) )
	else
		subplot( n_rows, NPL*n_columns+Xtra_col, NPL*n_columns+Xtra_col )
	end
	
	n_bins = 1000;
	[Gdist_y,Gdist_x] = hist(G,n_hist_bins);
    
	%this is a hack to deal with cases where the threshold linear terms
	%create a min value of G
	if Gdist_y(1) > 2*Gdist_y(2)
		Gdist_y(1) = 1.5*Gdist_y(2);
	end
    
	%cur_xrange = Gdist_x([1 end]);
	cur_xrange = [my_prctile(Gdist_x,2) my_prctile(Gdist_x,99)];
  xplot = Gdist_x((Gdist_x >= cur_xrange(1)) & (Gdist_x <= cur_xrange(2)));
	
	if strcmp(nim.spk_NL_type,'logexp')
		cur_y = nim.spk_NL_params(3)*log(1 + exp(nim.spk_NL_params(2)*xplot + theta));
	elseif strcmp(nim.spk_NL_type,'exp')
		cur_y = exp(xplot + theta);
  elseif strcmp(nim.spk_NL_type,'linear')
    cur_y = xplot + theta;
	else
		error('Unsupported spk NL type');
	end
	cur_y = cur_y/nim.stim_params(1).dt; %convert to correct firing rate units
	[ax,h1,h2] = plotyy(xplot,cur_y,Gdist_x,Gdist_y);
	%     set(h2,'Color','k')
	set(h1,'linewidth',1)
	yr = [min(cur_y) max(cur_y)];
	xlim(ax(1),cur_xrange)
	xlim(ax(2),cur_xrange);
	ylim(ax(1),yr);
    
	xlabel('Generating function')
	ylabel(ax(1),'Predicted firing rate','fontsize',12);
	%ylabel(ax(2),'Probability','fontsize',14)
	set(ax(2),'ytick',[]);
	title('Spiking NL','fontsize',14)
    
	%set horizontal position
  %  cp = get(fig_spk_nl,'Position');
  %  cp(1) = 1100;
  %  set(fig_spk_nl,'Position',cp);
end


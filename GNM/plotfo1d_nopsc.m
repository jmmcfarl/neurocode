function f1 = plotfo1d_nopsc(fo,max_rows,cmap_type,to_transpose)
% USAGE: test = plotfo1d(fo)
% plots fitted model for 1d stimulus (e.g. bar stimuli used in Touryan2001)

if nargin < 4
    to_transpose = 0;
end

% f1 = figure('Name',fo.mname);
f1 = figure();
colormap(gray);

nmods = length(fo.mods);

n_columns = ceil(nmods/max_rows);
n_rows = min(max_rows,nmods);

k_mat = get_k_mat(fo);
w = arrayfun(@(x) x.w,fo.mods);
k_mat = bsxfun(@times,k_mat,w);
mval = max(abs(k_mat(:)));
if strcmp(cmap_type,'equal')
    mval = 0.9*mval;
    mrange = [-mval mval]; %this is the overall centered caxis range
end
sdim = fo.stim_params.sdim;
for imod=1:nmods
    thismod = fo.mods(imod);
    
    cur_mod = imod;
    if cur_mod > max_rows
        cur_col = floor((cur_mod-1)/max_rows);
        cur_mod = mod((cur_mod-1),max_rows)+1;
    else
        cur_col = 0;
    end
    maxval = max(thismod.k);
    minval = min(thismod.k);
    mval = max(abs(maxval),abs(minval));
    if mval == 0
        mval = 1;
    end
    %     thismod.k = thismod.k*thismod.w;
    
    if strcmp(cmap_type,'centered')
        mrange = [-mval mval]; %this is a centered cmap range for this neuron
    elseif strcmp(cmap_type,'uncon')
        mrange = [minval maxval];
    end
    
    if unique(mrange)==0
        mrange = [-1 1];
    end
    
    subplot(n_rows,2*n_columns,(cur_mod-1)*2*n_columns+1+2*cur_col);
    plot2drf(thismod.k,sdim,mrange,to_transpose);
    %     caxis([min(thismod.k) max(thismod.k)])
    
    subplot(n_rows,2*n_columns,(cur_mod-1)*2*n_columns+2+2*cur_col);
    if strcmp(thismod.nltype,'lin')
        xx = linspace(thismod.nlx(1),thismod.nlx(end),1000);
        yy = xx;
        plot(xx,yy,'b')
    elseif strcmp(thismod.nltype,'quad')
         xx = linspace(thismod.nlx(1),thismod.nlx(end),1000);
       yy = xx.^2;
        plot(xx,yy,'b')
    elseif strcmp(thismod.nltype,'threshlin')
         xx = linspace(thismod.nlx(1),thismod.nlx(end),1000);
        yy = xx;
        yy( xx < 0) = 0;
        plot(xx,yy,'b')
    elseif strcmp(thismod.nltype,'uncon')
        xx = thismod.nlx;
        yy = thismod.nly - min(thismod.nly);
        plot(xx,yy,'b'); %grid();
    end
    axis tight;
    %     set(gca,'XTick',[])
    % set(gca,'YTick',[])
    
    hold on;
    if length(unique(thismod.k)) ~= 1
        if(isfield(thismod,'fox'))
            %         cur_foys = thismod.foys;
            %         cur_foys = cur_foys*scale;
            %         plot(thismod.foxs,cur_foys,'-r');
            
            cur_foy = thismod.foy;
            cumdist = cumsum(cur_foy)/sum(cur_foy);
            bound(1) = thismod.fox(find(cumdist < 0.001,1,'last'));
            bound(2) = thismod.fox(find(cumdist > 0.999,1,'first'));
            uset = find(xx >= bound(1) & xx <= bound(2));
            
            cla
            yy = yy-yy(uset(1));
            plot(xx,yy,'b')
            scale = max(yy(uset));
            cur_foy = thismod.foy;
            scale = scale/range(cur_foy);
            cur_foy = cur_foy*scale*0.9;
            plot(thismod.fox,cur_foy,'k');
            axis tight
            ylim([min(yy(uset)) max(yy(uset))])
            xlim(bound)
            box off
            %                 xlim([-4 4])
        end;
    end
    hold off;
    set(gca,'fontname','arial')
    % 	hline(0,'--k'); vline(0,'--k');
%     title(sprintf('weight: %0.4g',thismod.w));
cur_weight = range(ylim());
    title(sprintf('weight: %0.4g',cur_weight));

end


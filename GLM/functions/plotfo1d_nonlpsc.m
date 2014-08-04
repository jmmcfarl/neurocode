function f1 = plotfo1d_nonlpsc(fo,max_rows)
% USAGE: test = plotfo1d(fo)
% plots fitted model for 1d stimulus (e.g. bar stimuli used in Touryan2001)


f1 = figure('Name',fo.mname);
colormap(gray);

nmods = length(fo.mods);

k_mat = get_k_mat(fo);
max_k = max(abs(k_mat(:)));
z_range = [-max_k max_k];

n_columns = ceil(nmods/max_rows);
n_rows = min(max_rows,nmods);
for imod=1:nmods
    thismod = fo.mods(imod);
    sdim = thismod.SDIM;
    
    cur_mod = imod;
    if cur_mod > max_rows
        cur_col = floor((cur_mod-1)/max_rows);
        cur_mod = mod((cur_mod-1),max_rows)+1;
    else
        cur_col = 0;
    end
    
    subplot(n_rows,n_columns,(cur_mod-1)*n_columns+1+cur_col); plot2drf(thismod.k,sdim,z_range); 
       title(sprintf('weight: %0.4g',thismod.w));
end


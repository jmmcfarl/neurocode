function boxplot_capped(data,groups,bounds)

capped_data = data;
un_groups = unique(groups(~isnan(groups)));
n_groups = length(un_groups);

for ii = 1:n_groups
    curset = find(groups == un_groups(ii));
    
    cur_data = data(curset);
    cur_bounds = prctile(cur_data,bounds);
    cur_data(cur_data > cur_bounds(2)) = cur_bounds(2);
    cur_data(cur_data < cur_bounds(1)) = cur_bounds(1);
    
    capped_data(curset) = cur_data;
    
end

boxplot(capped_data,groups,'whisker',Inf);
function [] = mean_rate_order_view(data_matrix, x_vals, mean_rates, target_cells)

used_data = data_matrix(target_cells,:);
used_rates = mean_rates(target_cells);

[dummy,rate_order] = sort(used_rates);

cmap = colormap(jet(length(target_cells)));

for i = 1:length(target_cells)
    
   plot(x_vals,used_data(rate_order(i),:),'Color',cmap(i,:))
   hold on
    
end


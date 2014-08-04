function [] = anatomy_sorted_view(data_mat, sess_data, target_cells, property, xvals)

used_data = data_mat(target_cells,:);

for i = 1:length(target_cells)
    eval(['cell_prop(i) = sess_data(target_cells(i)).' property ';']);
end

[dummy, sorted_vals] = sort(cell_prop);

imagesc(xvals,1:length(target_cells),used_data(sorted_vals,:));
shading flat
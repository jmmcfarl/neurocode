function [mean_vec, se_vec] = get_data_mean_se(data_mat, target_cells)


used_data = data_mat(target_cells,:);

mean_vec = nanmean(used_data);
se_vec = nanstd(used_data)/sqrt(length(target_cells));

clear all
close all

load F:\WC_Germany\overall_EC\overall_EC_dir.mat
drive_letter = 'F';

%%
load overall_EC_spike_data
load overall_EC_heka_UDS_data
load corresponding_lfp_state_data
load overall_EC_UDS_sep

mec = find_struct_field_vals(sess_data,'region','MEC');
lec = find_struct_field_vals(sess_data,'region','LEC');
layer3 = find_struct_field_vals(sess_data,'layer','3');
layer2 = find_struct_field_vals(sess_data,'layer','2');
layer5 = find_struct_field_vals(sess_data,'layer','56');
l5mec = intersect(layer5,mec);



%%
uds_amp = upstate_mean - downstate_mean;
uds_amp(upstate_mean > -10) = nan;

for i = 1:length(sess_data)
    mean_updur(i) = nanmean(mp_updur_lfpc{i});
end

features(:,1) = uds_amp(:);
features(:,2) = up_rate(:);
% features(:,3) = median_up_lag(:);
features(:,3) = median_down_lag(:);
features(:,4) = mean_updur(:);
n_features = size(features,2);
bad_points = find(any(isnan(features),2));
features(bad_points,:) = [];
features = zscore(features);

%%
n_clusts = 2;
D = pdist(features);
Y = linkage(D,'average');
% dendrogram(Y,0,'colorthreshold','default');
T = cluster(Y,'MaxClust',n_clusts);

%%
cmap = colormap(jet(n_clusts));
cnt = 1;
for i = 1:n_features-1
    for j = i+1:n_features
        subplot(3,4,cnt)
        for c = 1:n_clusts
            plot(features(T==c,i),features(T==c,j),'.','color',cmap(c,:))
            hold on
            xlabel(sprintf('feature %d',i))
            ylabel(sprintf('feature %d',j))
            cnt = cnt + 1;
        end
    end
end
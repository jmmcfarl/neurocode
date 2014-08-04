function fin_autoclust = flip_cluster_assignment(probe_num,autoclust,target_expts,stored_spkxy)

new_autoclust = autoclust;

for ee = target_expts
    flip = [2 1];
    prev_clust_idx = new_autoclust(ee,probe_num).idx;
    new_autoclust(ee,probe_num).idx = flip(prev_clust_idx);
    new_autoclust(ee,probe_num).avg_wvfrm = new_autoclust(ee,probe_num).avg_wvfrm(:,flip);
    new_autoclust(ee,probe_num).std_wvfrm = new_autoclust(ee,probe_num).std_wvfrm(:,flip);
end

n_expts = length(target_expts);
n_cols = ceil(sqrt(n_expts));
n_rows = ceil(n_expts/n_cols);
f1 = figure('Name','previous');
f2 = figure('Name','new');
% f3 = figure('Name','previous');
% f4 = figure('Name','new');

for ee = 1:length(target_expts)
    figure(f1);
    subplot(n_rows,n_cols,ee);
    plot(stored_spkxy{target_expts(ee),probe_num}(autoclust(target_expts(ee),probe_num).idx == 1,1),stored_spkxy{target_expts(ee),probe_num}(autoclust(target_expts(ee),probe_num).idx == 1,2),'r.','markersize',0.5);
    hold on
    plot(stored_spkxy{target_expts(ee),probe_num}(autoclust(target_expts(ee),probe_num).idx == 2,1),stored_spkxy{target_expts(ee),probe_num}(autoclust(target_expts(ee),probe_num).idx == 2,2),'k.','markersize',0.5);
    
    figure(f2);
    subplot(n_rows,n_cols,ee);
    plot(stored_spkxy{target_expts(ee),probe_num}(new_autoclust(target_expts(ee),probe_num).idx == 1,1),stored_spkxy{target_expts(ee),probe_num}(new_autoclust(target_expts(ee),probe_num).idx == 1,2),'r.','markersize',0.5);
    hold on
    plot(stored_spkxy{target_expts(ee),probe_num}(new_autoclust(target_expts(ee),probe_num).idx == 2,1),stored_spkxy{target_expts(ee),probe_num}(new_autoclust(target_expts(ee),probe_num).idx == 2,2),'k.','markersize',0.5);

%     figure(f1);
%     subplot(n_rows,n_cols,ee);
%     plot(stored_spkxy{target_expts(ee),probe_num}(autoclust(target_expts(ee),probe_num).idx == 1,1),stored_spkxy{target_expts(ee),probe_num}(autoclust(target_expts(ee),probe_num).idx == 1,2),'r.','markersize',0.5);
%     hold on
%     plot(stored_spkxy{target_expts(ee),probe_num}(autoclust(target_expts(ee),probe_num).idx == 2,1),stored_spkxy{target_expts(ee),probe_num}(autoclust(target_expts(ee),probe_num).idx == 2,2),'k.','markersize',0.5);
%     
%     figure(f2);
%     subplot(n_rows,n_cols,ee);
%     plot(stored_spkxy{target_expts(ee),probe_num}(new_autoclust(target_expts(ee),probe_num).idx == 1,1),stored_spkxy{target_expts(ee),probe_num}(new_autoclust(target_expts(ee),probe_num).idx == 1,2),'r.','markersize',0.5);
%     hold on
%     plot(stored_spkxy{target_expts(ee),probe_num}(new_autoclust(target_expts(ee),probe_num).idx == 2,1),stored_spkxy{target_expts(ee),probe_num}(new_autoclust(target_expts(ee),probe_num).idx == 2,2),'k.','markersize',0.5);
end

ui = input('Use flip (y)?\n','s');
if strcmpi(ui,'y')
    fin_autoclust = new_autoclust;
else
    fin_autoclust = autoclust;
end
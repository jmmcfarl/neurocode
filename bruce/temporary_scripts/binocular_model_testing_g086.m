%%
clear all
cd ~/Analysis/bruce/G086/ET_final/
load binoc_datadump_dual

% load cur_best_binoc_mods

flen = 13; %number of lags
dt = 0.01; %in sec
dx = 0.0283; %in deg
spatial_usfac = 2; %spatial up-sampling
use_nPix_us = size(left_stim,2); %number of (upsampled) pixels
se_stim_params = NMMcreate_stim_params([flen use_nPix_us]);
Xmat = [create_time_embedding(left_stim,se_stim_params) create_time_embedding(right_stim,se_stim_params)];

%%
%make custom regularization matrix
%create L2 mat for left eye
L2_params = create_L2_params([],[1 flen*use_nPix_us],[flen use_nPix_us],2,3,[0 0]);
L2_mat_us = generate_L2_mat(L2_params,2*flen*use_nPix_us);
%add L2 mat for right eye
L2_params = create_L2_params([],[flen*use_nPix_us+1 2*flen*use_nPix_us],[flen use_nPix_us],2,3,[0 0]);
L2_mat_us = L2_mat_us + generate_L2_mat(L2_params,2*flen*use_nPix_us);

cc = 98; %SUs are cc > 96 (99 is best)

%get relevant spike and stim variables
cur_Robs = Robs_mat(:,cc);
cc_uinds = find(~isnan(cur_Robs));
cur_Robs = cur_Robs(cc_uinds);
cur_Xmat = Xmat(used_inds(cc_uinds),:);

nim_stim_params = NMMcreate_stim_params([flen use_nPix_us*2]);
% init_lambda_custom = 500;
init_lambda_custom = 20;
nim_reg_params = NMMcreate_reg_params('lambda_custom',init_lambda_custom);

optim_params.optTol = 1e-4;
optim_params.progTol = 1e-6;
optim_params.maxIter = 200;

mod_signs = [1 1 1 -1 -1 -1];
NL_types = {'threshlin','quad','quad','threshlin','quad','quad','quad'};
% NL_types = {'lin','quad','quad','quad','quad','quad'};
GQM = NMMinitialize_model(nim_stim_params,mod_signs,NL_types,nim_reg_params);
GQM = NMMfit_filters(GQM,cur_Robs,cur_Xmat,[],[],0,optim_params,L2_mat_us);

%adjust regularization and refit
% base_lambda_custom = 700; %300
% base_lambda_L1 = 25; %15
base_lambda_custom = 50; %300
base_lambda_L1 = 5; %15
[LL, penLL, pred_rate, G, gint] = NMMmodel_eval(GQM,cur_Robs,cur_Xmat);
GQM = NMMadjust_regularization(GQM,[],'lambda_custom',base_lambda_custom./var(gint)');
GQM = NMMadjust_regularization(GQM,[],'lambda_L1',base_lambda_L1./std(gint)');
GQM = NMMfit_filters(GQM,cur_Robs,cur_Xmat,[],[],0,optim_params,L2_mat_us);

% GQM = NMMadjust_regularization(GQM,4,'lambda_custom',4e3);
% % % lambda_custom = 10;
% % % lambda_L1 = 2;
% % GQM2 = NMMadjust_regularization(GQM,[1 4],'lambda_custom',[2e3 2e3]);
% GQM2 = NMMadjust_regularization(GQM,[1 4],'lambda_L1',[75 75]);
% % % GQM2.mods(end).filtK = 0.01*randn(size(GQM2.mods(end).filtK));
% GQM2 = NMMfit_filters(GQM2,cur_Robs,cur_Xmat,[],[],0,optim_params,L2_mat_us);

NMMdisplay_model(GQM);

all_GQM(cc) = GQM;
%%
close all
fig_dir = '/home/james/Desktop/grantfigs/';

cc = 98;
cur_GQM = all_GQM(cc);
cur_Robs = Robs_mat(:,cc);
cc_uinds = find(~isnan(cur_Robs));
cur_Robs = cur_Robs(cc_uinds);
cur_Xmat = Xmat(used_inds(cc_uinds),:);

[LL, penLL, pred_rate, G, gint] = NMMmodel_eval(cur_GQM,cur_Robs,cur_Xmat);


stim_filters = [cur_GQM.mods(:).filtK];
left_filters = stim_filters(1:flen*32,:);
right_filters = stim_filters((flen*32+1):end,:);
stim_mod_signs = [cur_GQM.mods(:).sign];
mod_stim_params = cur_GQM.stim_params;
mod_stim_params.stim_dims(2) = mod_stim_params.stim_dims(2)/2;
left_filt_data = get_filter_properties_v2(left_filters,mod_stim_params,dx);
right_filt_data = get_filter_properties_v2(right_filters,mod_stim_params,dx);

n_filters = size(left_filters,2);
left_filters = reshape(left_filters,[flen use_nPix_us n_filters]);
right_filters = reshape(right_filters,[flen use_nPix_us n_filters]);

f1 = figure;
t_ax = ((0:flen-1)*dt + dt/2)*1e3;
p_ax = (1:use_nPix_us)*dx - use_nPix_us/2*dx;
p_ax = -0.43 - p_ax;
% p_ax = 0.34 - p_ax;
n_hist_bins = 500; %internal parameter determining histogram resolution
for ii = 1:n_filters
    subplot(n_filters,3,3*(ii-1)+1);
    imagesc(p_ax,t_ax,squeeze(left_filters(:,:,ii))); set(gca,'ydir','normal');
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    subplot(n_filters,3,3*(ii-1)+2);
    imagesc(p_ax,t_ax,squeeze(right_filters(:,:,ii)));set(gca,'ydir','normal');
    caxis([-cam cam]);
    colormap(gray);
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    subplot(n_filters,3,3*(ii-1)+3);
    [gendist_y,gendist_x] = hist(gint(:,ii),n_hist_bins);
    thismod = cur_GQM.mods(ii);
    if strcmp(thismod.NLtype,'nonpar')
        cur_modx = thismod.NLx; cur_mody = thismod.NLy;
    elseif strcmp(thismod.NLtype,'lin')
        cur_modx = gendist_x; cur_mody = cur_modx;
    elseif strcmp(thismod.NLtype,'quad')
        cur_modx = gendist_x;
        cur_mody = cur_modx.^2;
    elseif strcmp(thismod.NLtype,'threshlin')
        cur_modx = gendist_x;
        cur_mody = cur_modx;
        cur_mody(cur_mody < 0) = 0;
    end
    cur_xrange = cur_modx([1 end]);
    ma = max(abs(cur_xrange))*0.95;
    cur_xrange = [-ma ma];
    [ax,h1,h2] = plotyy(cur_modx,cur_mody,gendist_x,gendist_y);
    if strcmp(thismod.NLtype,'nonpar')
        %                 set(h1,'Marker','o');
    end
    %         set(h2,'Color','k')
    set(h1,'linewidth',1)
    xlim(ax(1),cur_xrange)
    xlim(ax(2),cur_xrange);
    ylim(ax(1),[min(cur_mody) max(cur_mody)]);
    set(ax(2),'ytick',[])
    yl = ylim();
    line([0 0],yl,'color','k','linestyle','--');
    if thismod.sign == 1
        title('Excitatory');
    else
        title('Suppressive');
    end
%     ylabel(ax(1),'Subunit output','fontsize',12);
%     ylabel(ax(2),'Probability','fontsize',12)
end


ufilts = 2:3;
best_tfreq = mean(left_filt_data.FFt(ufilts));
best_tfreq_ind = find(left_filt_data.ax_Ft >= best_tfreq,1);
best_xfreq = mean(left_filt_data.FFx(ufilts));
best_xfreq_ind = find(left_filt_data.ax_Fx >= best_xfreq,1);
left_FFt = squeeze(mean(left_filt_data.filt_FFts(:,:,ufilts),3));
right_FFt = squeeze(mean(right_filt_data.filt_FFts(:,:,ufilts),3));

f2 = figure; hold on
plot(left_filt_data.ax_Fx,left_FFt(best_tfreq_ind,:));
plot(left_filt_data.ax_Fx,right_FFt(best_tfreq_ind,:),'r');
xlim([-12 12]);
xlabel('Spatial frequency (cyc/deg)');
legend('Left eye','Right eye');

% f3 = figure; hold on
% plot(left_filt_data.ax_Ft/dt,left_FFt(:,best_xfreq_ind));
% plot(left_filt_data.ax_Ft/dt,right_FFt(:,best_xfreq_ind),'r');
% xlim([-40 40]);
% xlabel('Temporal frequency (cyc/deg)');
% legend('Left eye','Right eye');

% f4 = figure;
% imagesc(left_filt_data.ax_Fx,left_filt_data.ax_Ft/dt,left_FFt);shading flat
% set(gca,'ydir','normal');
% xlim([-12 12]);
% ylim([-35 35]);
% yl = ylim(); xl = xlim();
% line(xl,[0 0],'color','w');
% line([0 0],yl,'color','w');
% line(xl,[best_tfreq best_tfreq]/dt,'color','r')
% xlabel('Spatial frequency (cyc/deg)');
% ylabel('Temporal frequency (cyc/s)');

% fname = [fig_dir sprintf('C%d_modfilts3.pdf',cc)];
% fig_width = 7;
% rel_height = n_filters/3;
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

% fname = [fig_dir sprintf('C%d_spatialfreqdual.pdf',cc)];
% fig_width = 5;
% rel_height = 0.8;
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f2);
% 
% fname = [fig_dir sprintf('C%d_FFTdual.pdf',cc)];
% fig_width = 5;
% rel_height = 0.8;
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f2);

%%
cur_GQM = all_GQM(cc);
modFilts = reshape([cur_GQM.mods(:).filtK],[flen use_nPix_us 2 length(cur_GQM.mods)]);
temp_kerns = squeeze(std(modFilts,[],2));
% temp_kerns = mean(temp_kerns,2);
efilts = find([cur_GQM.mods(:).sign] == 1);
ifilts = find([cur_GQM.mods(:).sign] == -1);
f5 = figure;
subplot(2,1,1)
% plot(t_ax,temp_kerns(:,efilts),'b');
% hold on
% plot(t_ax,temp_kerns(:,ifilts),'r');
plot(t_ax,squeeze(mean(temp_kerns(:,1,efilts),3)),'b','linewidth',2);
hold on
plot(t_ax,squeeze(mean(temp_kerns(:,1,ifilts),3)),'r','linewidth',2);
xlim([0 max(t_ax)]);
xlabel('Time lag (ms)');
title('Left eye')
subplot(2,1,2)
plot(t_ax,squeeze(mean(temp_kerns(:,2,efilts),3)),'b','linewidth',2);
hold on
plot(t_ax,squeeze(mean(temp_kerns(:,2,ifilts),3)),'r','linewidth',2);
    xlim([0 max(t_ax)]);
    xlabel('Time lag (ms)');
title('Right eye')
    fname = [fig_dir sprintf('C%d_tempkerns.pdf',cc)];
    fig_width = 5;
    rel_height = 1.6;
    figufy(f5);
    exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    % close(f2);

%% COMPUTE PREDICTED DISPARITY TUNING
cc = 98;
cur_GQM = all_GQM(cc);

poss_dxs = -15:15; %range of disparities (in pxiels)
Nframe_test = 1e4; %about 5000 RLS frames seems to be sufficient
dds = 12/100; %sparsity of test stim (real data is either 12/100 or 67/100)

full_nPix = 36;
use_nPix = 16;
full_nPix_us = full_nPix*spatial_usfac;
%repeat for up-sampled versions of the Xmatrix
buffer_pix = floor((full_nPix - use_nPix)/2);
[Xinds_up,~] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
cur_use_pix = (1/spatial_usfac:1/spatial_usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

for xx = 1:length(poss_dxs)
    fprintf('Computing avg response to disparity %d/%d\n',xx,length(poss_dxs));
    %make an RLS stimulus
    rand_stim = zeros(Nframe_test,full_nPix_us/2);
    rand_mat = rand(Nframe_test,full_nPix_us/2);
    rand_stim(rand_mat <= dds/2) = -1;
    rand_stim(rand_mat >= 1-(dds/2)) = 1;
    
    %up-sample bars by factor of spatial_usfac
    rand_stimmat_up = zeros(size(rand_stim,1),full_nPix_us);
    for ii = 1:size(rand_stim,2)
        for jj = 1:spatial_usfac
            rand_stimmat_up(:,spatial_usfac*(ii-1)+jj) = rand_stim(:,ii);
        end
    end
    
    rand_stimmat_up = rand_stimmat_up/std(rand_stimmat_up(:));
    
    %make left and right stims perfectly correlated but spatially shifted
    rand_left_stim = rand_stimmat_up;
    rand_right_stim = shift_matrix_Nd(rand_left_stim,poss_dxs(xx),2);
    left_X = create_time_embedding(rand_left_stim,NMMcreate_stim_params([flen full_nPix_us]));
    right_X = create_time_embedding(rand_right_stim,NMMcreate_stim_params([flen full_nPix_us]));
    sim_Xmat = [left_X(:,use_kInds_up) right_X(:,use_kInds_up)];
    
    %compute average firing rate output predicted by the model
    [~,~,mod_out] = NMMmodel_eval(cur_GQM,[],sim_Xmat);
    sim_disp_tune(xx) = mean(mod_out);
    
end

f1 = figure();
plot(poss_dxs*dx,sim_disp_tune/dt,'k','linewidth',2)
xlabel('Vertical disparity (deg)');
ylabel('Predicted rate (Hz)');
axis tight
xlim([-0.4 0.4]);
yl = ylim();
line([0 0],yl,'color','k');
fname = [fig_dir sprintf('C%d_disparity2.pdf',cc)];
fig_width = 5;
rel_height = 0.8;
figufy(f1);
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

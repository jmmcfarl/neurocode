clear all
close all
cd /home/james/Data/bruce/7_15_12/G029
% load ./Expt1_newcompiled_data_flen8_d1_new 
load ./corrected_Xmat_flen6_fullres
fullX = fullX_cor/std(fullX_cor(:));
% fullX = fullX/std(fullX(:));

load ./expt1_eyecor_d1p25_nosac_v2 gabor*
gabor_params = gabor_params_f{end};
load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
load ./oned_fixation_fits_v3.mat

 %%
NT = size(fullX,1);
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);

idier = [full_expt_vec full_trial_vec];
unique_ids = unique(idier,'rows');
n_trials = size(unique_ids,1);
xv_frac =0.2;
xv_num = round(xv_frac*n_trials);
xv_set = randperm(n_trials);
xv_set(xv_num+1:end) = [];
tr_set = setdiff(1:n_trials,xv_set);

xv_inds = [];
for ii = 1:xv_num
   cur_inds = find(full_expt_vec == unique_ids(xv_set(ii),1) & full_trial_vec == unique_ids(xv_set(ii),2));
   xv_inds = [xv_inds; cur_inds];
end
tr_inds = setdiff(1:NT,xv_inds);

%%
close all
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
sdim = length(xpatch_inds);
flen = 6;
dt = 118*2/1e4;

for t = 1:96
% for t = single_units
% t = 2;
gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);

gabor_emp1_vec = [gabor_emp1(:)'; zeros(flen-1,sdim^2)];
gabor_emp1_mat = makeStimRows(gabor_emp1_vec,flen);

gabor_emp2_vec = [gabor_emp2(:)'; zeros(flen-1,sdim^2)];
gabor_emp2_mat = makeStimRows(gabor_emp2_vec,flen);

gabor_outs1 = fullX*gabor_emp1_mat';
gabor_outs2 = fullX*gabor_emp2_mat';

%initialize model
defmod.lambda_L2x = 0;
defmod.lambda_dT = 0;
defmod.lambda_d2X = 0;
defmod.lambda_L1x = 0;

%initialize model
defmod2.lambda_L2x = 5000;
defmod2.lambda_dT = 0;
defmod2.lambda_d2X = 0;
defmod2.lambda_L1x = 0;


stim_params.flen = flen;
stim_params.spatial_dims = 0;
stim_params.sdim = 1;


tr_spikebins = convert_to_spikebins(full_binned_spks(tr_inds,t));
xv_spikebins = convert_to_spikebins(full_binned_spks(xv_inds,t));


% X2 = sqrt(gabor_outs1.^2 + gabor_outs2.^2);
X2 = gabor_outs1.^2 + gabor_outs2.^2;
[NT,klen] = size(X2);
k0 = zeros(klen,1);
cur_ndims = 1;
n_inits = 1;
clear temp_gnm LL xvLL
for in = 1:n_inits
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.05*randn(klen,cur_ndims);
    for i = 1:cur_ndims
        kern_types{i} = 'lin';
    end
    temp_gnm(in) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    temp_gnm(in) = fitGNM_filters(temp_gnm(in),X2(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6);
    LL(in) = temp_gnm(in).LL;
    xvLL(in) = getLL_GNM(temp_gnm(in),X2(xv_inds,:),xv_spikebins,'none');
end
[~,best] = min(xvLL);
gnmold(t) = temp_gnm(best);
best_gnmold_xvLL(t) = xvLL(best);
best_gnmold_LL(t) = LL(best);

k = get_k_mat(gnmold(t));
% k0 = [k; k];

X = [gabor_outs1 gabor_outs2];
[NT,klen] = size(X);
k0 = zeros(klen,1);
stim_params.spatial_dims = 1;
stim_params.sdim = 2;

time_ax = [1:flen 1:flen];
clear i
phase_ax = [ones(1,flen) i*ones(1,flen)];

cur_ndims = 2;
n_inits = 1;
clear temp_gnm LL xvLL
for in = 1:n_inits
    init_signs = ones(cur_ndims,1);
%     init_signs(end) = -1;
    init_kerns = 0.05*randn(klen,cur_ndims);
    for i = 1:cur_ndims
        kern_types{i} = 'quad';
    end
%     kern_types{1} = 'lin';
    temp_gnm(in) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    temp_gnm(in) = fitGNM_filters(temp_gnm(in),X(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6);
    LL(in) = temp_gnm(in).LL;
    getLL_GNM(temp_gnm(in),X(tr_inds,:),tr_spikebins,'none')
    xvLL(in) = getLL_GNM(temp_gnm(in),X(xv_inds,:),xv_spikebins,'none');

    temp_gnm2(in) = createGNM(init_kerns,init_signs,kern_types,defmod2,stim_params);
    temp_gnm2(in) = fitGNM_filters(temp_gnm2(in),X(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6);
    LL2(in) = temp_gnm2(in).LL;
    getLL_GNM(temp_gnm2(in),X(tr_inds,:),tr_spikebins,'none')
    xvLL2(in) = getLL_GNM(temp_gnm2(in),X(xv_inds,:),xv_spikebins,'none');
end
[~,best] = min(xvLL);
best_gnm_LL(t) = LL(best);
gnm(t) = temp_gnm(best);
best_gnm_xvLL(t) = xvLL(best);

[~,best] = min(xvLL2);
best_gnm_LL2(t) = LL2(best);
gnm2(t) = temp_gnm2(best);
best_gnm_xvLL2(t) = xvLL2(best);

kmat = get_k_mat(gnm(t));

%%
% conv_mat = [gabor_emp1_mat' gabor_emp2_mat'];
% k_pix = conv_mat*kmat;
% % plotfilterbank(k_pix,sdim,1:sdim^2)
% 
% clear i
% form_k = kmat(1:flen,:) + kmat((flen+1):end,:)*i;
% k_amps = abs(form_k);
% k_phase = angle(form_k);
% 
% %
% 
% avg_best_theta = gabor_params(t,3)+pi/2;
% avg_best_x = gabor_params(t,1);
% avg_best_y = gabor_params(t,2);
% 
% fsdim = sdim^2;
% % f1 = figure;
% 
% nxax = xax(xpatch_inds);
% nyax = yax(ypatch_inds);
% n_filts = size(k_pix,2);
% k_norm = zeros(n_filts,1);
% for ii = 1:n_filts
%     k_norm(ii) = norm(k_pix(:,ii));
% end
% [~,ord] = sort(k_norm,'descend');
% k_pix = k_pix(:,ord);
% k_norm = k_norm(ord);
% [~,best_slice_ids] = get2dMaxSlices(k_pix,flen,sdim,1:(sdim^2));
% figure
% for ii = 1:n_filts
%     cur_filtmat = reshape(k_pix(:,ii),flen,fsdim);
%     
%     filt_projprof = project_2drf(cur_filtmat(:),avg_best_theta,sdim);
%     
%     subplot(n_filts,2,(ii-1)*2+1)
%     
%     imagesc(nxax,nyax,reshape(cur_filtmat(best_slice_ids(ii),:),sdim,sdim));colormap(jet);
%     zmax = max(abs(cur_filtmat(best_slice_ids(ii),:))); caxis([-zmax zmax]); %colorbar
%     title(sprintf('Norm: %.3f',k_norm(ii)))
%     set(gca,'ydir','normal')
%     hold on
%     plot(avg_best_x,avg_best_y,'wo','linewidth',2)
%     cur_linex = linspace(nxax(1),nxax(end),50);
%     cur_liney = tan(avg_best_theta)*(cur_linex - avg_best_x)+avg_best_y;
%     plot(cur_linex,cur_liney,'color','w')
%     xlim(nxax([1 end]));ylim(nyax([1 end]));
%     
%     subplot(n_filts,2,(ii-1)*2+2)
%     imagesc(1:sdim,-(0:(flen-1))*dt,filt_projprof);
%     title(sprintf('Norm: %.3f',k_norm(ii)))
%     xm = max(abs(filt_projprof(:)));
%     caxis([-xm xm]/1.5);
%     title('Space-time projection')   
%     
% end
% t
% set(gcf,'Position',[600 200 800 600])
% 
% figure
% plot(norm_mu_ori_tuning(t,:))
% set(gcf,'Position',[1400 200 500 400])
% 
% figure
% plot(k)
% set(gcf,'Position',[200 200 500 400])
% 
% disp('XVLL')
% fprintf('Temporal: %.3f\n',best_gnm_xvLL(t));
% fprintf('Spatial: %.3f\n',best_gnmold_xvLL(t));
% disp('LL')
% fprintf('Temporal: %.3f\n',best_gnm_LL(t));
% fprintf('Spatial: %.3f\n',best_gnmold_LL(t));
% 
% pause
% close all
% 

end

% save spatiotemporal_mods gnm
%%
% % [~,~,~,~,g] = getLL_GNM(gnm(t),X(tr_inds,:),tr_spikebins,'none');
% % quadr_mod_ov = fitGNM_spkNL(gnm(t),g,tr_spikebins,1);
% for t = single_units
% % t = 69;
% gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
% gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
% 
% mv = max(abs(gabor_emp1(:)));
% t
% imagesc(gabor_emp1)
% caxis([-mv mv]);
% colorbar
% pause
% clf
% end
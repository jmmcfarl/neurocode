clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
%%
stim_fs = 100; %in Hz
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));
good_cellnum = CellList(1,good_sus,1);
use_sus = 1:96;

%%
for i = 1:length(Expts)
    if strcmp(Expts{i}.Header.expname,'grating.OpXseRC') | strcmp(Expts{i}.Header.expname,'grating.OpRC')
        is_bar_expt(i) = 1;
    else
        is_bar_expt(i) = 0;
    end
    
    if strcmp(Expts{i}.Stimvals.Bs,'image')
        expt_image_back(i) = 1;
    else
        expt_image_back(i) = 0;
    end
    
    expt_sim_sacs(i) = Expts{i}.Stimvals.ijump;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
    
end
expt_bar_ori(expt_bar_ori == -45) = 135;

load ./all_un_bar_pos
n_bar_pos = size(all_un_bar_pos,1);
%%
flen = 15;
beg_buffer = round(stim_fs*0.15);
bar_oris = [0 45 90 135];
for EE = [1 2 3 4];
    cd ~/Data/bruce/G081/
    fprintf('Analyzing %d ori expts\n',bar_oris(EE));
    
    %     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE));
    
    %expts with X deg bars and gray back (sim sacs)
    cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris(EE));
%     cur_expt_set(cur_expt_set==2) = []; %this expt had continuous bar luminance
    cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
    cur_expt_set(cur_expt_set<= 8) = []; %this expt had continuous bar luminance
    
    un_bar_pos = all_un_bar_pos(:,EE);
    %%
    
    all_stim_times = [];
    all_rel_stimes = [];
    all_rel_etimes = [];
    all_phase = [];
    all_Op = [];
    all_bar_mat = [];
    all_used_inds = [];
    all_binned_spks = [];
    all_exptvec = [];
    all_trialvec = [];
    for ee = 1:length(cur_expt_set);
        fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
        cur_expt = cur_expt_set(ee);
        fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
        load(fname);
        
        n_trials = length(Expts{cur_expt}.Trials);
        trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
        trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
        trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
        
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_expt}.Trials(tt).Start/1e4;
            cur_Op = Expts{cur_expt}.Trials(tt).Op;
            cur_phase = Expts{cur_expt}.Trials(tt).ph;
            
            cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(tt).End(end)/1e4];
            cur_binned_spks = nan(length(cur_stim_times),length(use_sus));
            for cc = 1:length(use_sus)
                cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges);
                cur_binned_spks(:,cc) = cur_hist(1:end-1);
            end
            
            cur_bar_mat = zeros(length(cur_stim_times),n_bar_pos);
            for bb = 1:n_bar_pos
                cur_set = find(cur_Op==un_bar_pos(bb));
                pset = cur_phase(cur_set) == 0;
                nset = cur_phase(cur_set) == pi;
                cur_bar_mat(cur_set(pset),bb) = 1;
                cur_bar_mat(cur_set(nset),bb) = -1;
            end
            bar_Xmat = makeStimRows(cur_bar_mat,flen);
            cur_used_inds = ones(length(cur_stim_times),1);
            cur_used_inds(1:flen) = 0;
            cur_used_inds(1:beg_buffer) = 0;
            
            all_stim_times = [all_stim_times; cur_stim_times];
            all_rel_stimes = [all_rel_stimes; cur_stim_times- trial_start_times(tt)];
            all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_stim_times];
            all_Op = [all_Op; cur_Op];
            all_phase = [all_phase; cur_phase'];
            all_used_inds = [all_used_inds; cur_used_inds];
            all_bar_mat = [all_bar_mat; bar_Xmat];
            all_binned_spks = [all_binned_spks; cur_binned_spks];
            all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
        end
        
    end
    
    %%
    [c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
    n_trials = length(ia);
    
    xv_frac = 0.2;
    n_xv_trials = round(n_trials*xv_frac);
    xv_set = randperm(n_trials);
    xv_set(n_xv_trials+1:end) = [];
    xv_inds = find(ismember(ic,xv_set));
    tr_inds = find(~ismember(ic,xv_inds))';
    
    tr_inds(all_used_inds(tr_inds) == 0) = [];
    xv_inds(all_used_inds(xv_inds) == 0) = [];
    
    Xmat = all_bar_mat;
    
    lin_X = zeros(length(all_stim_times),length(cur_expt_set)-1);
    for i = 1:length(cur_expt_set)-1
        cur_set = find(all_exptvec==i);
        lin_X(cur_set,i) = 1;
    end
    null_rate = zeros(96,length(all_stim_times));
    for i = 1:length(cur_expt_set)
        cur_set = find(all_exptvec==i);
        cur_set2 = cur_set(ismember(cur_set,tr_inds));
        expt_avg_rate(i,:) = mean(all_binned_spks(cur_set2,:));
        null_rate(:,cur_set) = repmat(expt_avg_rate(i,:),[1 1 length(cur_set)]);
    end
    true_avg_rate(EE,:) = mean(all_binned_spks(tr_inds,:));
    
    %%
%     cd ~/James_scripts/bruce/G081/model_fits2/
%     
%     stim_params.spatial_dims = 1;
%     stim_params.sdim = n_bar_pos;
%     SDIM = n_bar_pos;
%     stim_params.flen = flen;
%     klen = length(un_bar_pos)*flen;
%     clear defmod
%     
%     NLtype = 0;
%     for cc = 1:length(use_sus)
%         fprintf('Cell %d of %d\n',cc,length(use_sus));
%         
%         Xmat_tr = Xmat(tr_inds,:);
%         tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,cc));
%         
%         sta(EE,cc,:) = mean(Xmat_tr(tr_spkbns,:)) - mean(Xmat_tr);
%         
%         nneg = 6;
%         npos = 6;
%         spike_cond_stim = Xmat_tr(tr_spkbns,:);
%         cur_sta = squeeze(sta(EE,cc,:));
%         cur_sta = cur_sta/norm(cur_sta);
%         proj_mat = cur_sta'/(cur_sta*cur_sta')*cur_sta;
% %         stim_proj = Xmat_tr - Xmat_tr*proj_mat;
%         stim_proj =  Xmat_tr;
%         stvcv = cov(stim_proj(tr_spkbns,:));  utvcv = cov(stim_proj);
%         [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%         evals = diag(evals);
%         all_evals(EE,cc,:) = evals;
%         
%         stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%         mf = max(npos,nneg);
%         figure('visible','off');
%         subplot(3,mf,1)
%         imagesc(reshape(cur_sta,flen,SDIM));
%         for i = 1:npos
%             subplot(3,mf,mf+i)
%             imagesc(reshape(stcs(:,i),flen,SDIM));
%         end
%         for i = 1:nneg
%             subplot(3,mf,2*mf+i)
%             imagesc(reshape(stcs(:,npos+i),flen,SDIM));
%         end
%         colormap(gray)
%         fname = sprintf('Unit%d_Ori%d_STC',use_sus(cc),bar_oris(EE));
%         fillPage(gcf,'papersize',[15 9]);
%         print(fname,'-dpng');
%         close 
%         
%         defmod.lambda_L1x = 10;
%         defmod.lambda_d2XT = 200;
%         kern_types{1} = 'lin';
%         init_kerns = 0.01*randn(klen,1);
%         init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
%         init_linK = zeros(length(cur_expt_set)-1,1);
%         glm = createGNM_v2(init_kerns,1,kern_types,init_linK,defmod,stim_params);
%         glm_fit(EE,cc) = fitGNM_filters_v2(glm,Xmat(tr_inds,:),lin_X(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6);
%         
%         defmod.lambda_L1x = 10;
%         defmod.lambda_d2XT = 200;
%         npq = 2;
%         nnq = 1;
%         kern_types{1} = 'lin';
%         for i = 2:(1+npq+nnq)
%             kern_types{i} = 'quad';
%         end
%         init_kerns = 0.01*randn(klen,1+npq+nnq);
%         init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
%         init_linK = zeros(length(cur_expt_set)-1,1);
%         quad = createGNM_v2(init_kerns,[1 ones(1,npq) -1*ones(1,nnq)],kern_types,init_linK,defmod,stim_params);
%         for i = 2:(1+npq+nnq)
%             quad.mods(i).lambda_L1x = 1;
%             quad.mods(i).lambda_d2XT = 5;
%         end
%         quad_fit(EE,cc) = fitGNM_filters_v2(quad,Xmat(tr_inds,:),lin_X(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6);
%  
%         plotfo1d_nopsc(quad_fit(EE,cc),2,'centered',0,0);colormap(jet);
%         fname = sprintf('Unit%d_Ori%d_Quad',use_sus(cc),bar_oris(EE));
%         fillPage(gcf,'papersize',[15 9]);
%         print(fname,'-dpng');
%         close 
%         
%         
%         %         defmod.lambda_L1x = 2;
%         %         defmod.lambda_d2XT = 50;
%         %         npq = 8;
%         %         nnq = 2;
%         %         for i = 1:(npq+nnq)
%         %             kern_types{i} = 'threshlin';
%         %         end
%         %         init_kerns = 0.01*randn(klen,npq+nnq);
%         %         init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
%         % %         init_linK = zeros(length(cur_expt_set)-1,1);
%         %         gnm = createGNM_v2(init_kerns,[ones(1,npq) -1*ones(1,nnq)],kern_types,init_linK,defmod,stim_params);
%         % %         gnm_fit(EE,cc) = fitGNM_filters_v2(gnm,Xmat(tr_inds,:),lin_X(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6);
%         %         gnm_fit(EE,cc) = fitGNM_filters(gnm,Xmat(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6);
%         %
%         %         gnm_fit(EE,cc) = adjust_all_reg(gnm_fit(EE,cc),'lnl2',500); %1000
%         %         gnm_fit(EE,cc) = setGNM_NLBFs(gnm_fit(EE,cc),Xmat(tr_inds,:));
%         %         gnm_fit(EE,cc) = adjust_all_reg(gnm_fit(EE,cc),'nltype','uncon');
%         %         gnm_fit(EE,cc) = adjust_all_reg(gnm_fit(EE,cc),'nlmon',1);
%         %         gnm_fit(EE,cc) = fitGNM_internal_NLs(gnm_fit(EE,cc),Xmat(tr_inds,:),tr_spkbns,0,2);
%         %         gnm_fit(EE,cc) = fitGNM_filters(gnm_fit(EE,cc),Xmat(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6);
%         %         [~, ~, ~, ~, g] = getLL_GNM(gnm_fit(EE,cc),Xmat(tr_inds,:),tr_spkbns,'none');
%         %         gnm_fit(EE,cc) = fitGNM_spkNL(gnm_fit(EE,cc),g,tr_spkbns,0);
%         
%         
%         Robs = all_binned_spks(xv_inds,cc);
%         xv_spkbns = convert_to_spikebins(all_binned_spks(xv_inds,cc));
%         [glm_xvLL(EE,cc), pnll, lpen, prate, g, int_g] = getLL_GNM_v2(glm_fit(EE,cc),Xmat(xv_inds,:),lin_X(xv_inds,:),xv_spkbns,'none');
%         [quad_xvLL(EE,cc), pnll, lpen, prate, g, int_g] = getLL_GNM_v2(quad_fit(EE,cc),Xmat(xv_inds,:),lin_X(xv_inds,:),xv_spkbns,'none');
%         
%         xv_null(EE,cc) = -sum(Robs.*log(null_rate(cc,xv_inds)') - null_rate(cc,xv_inds)')/sum(Robs);
%         
%     end
end

%%
% save bar_avg_rates true_avg_rate
%%
save all_unit_stim_fits glm_fit flen quad_fit xv_null *xvLL tr_* 
%%
for cc = 1:96
    if ismember(cc,good_sus)
        fprintf('SU: %d\n',cc);
    else
        fprintf('MU: %d\n',cc);
    end
    
    subplot(2,2,1)
    imagesc(un_bar_pos,1:flen,reshape(get_k_mat(glm_fit(1,cc)),flen,n_bar_pos))
    ca = max(abs(caxis())); caxis([-ca ca]*0.9);
    title(sprintf('Ori: 0, Mag: %.2d',ca));
    
    subplot(2,2,2)
    imagesc(un_bar_pos,1:flen,reshape(get_k_mat(glm_fit(2,cc)),flen,n_bar_pos))
    ca = max(abs(caxis())); caxis([-ca ca]*0.9);
    title(sprintf('Ori: 45, Mag: %.2d',ca));
    
    subplot(2,2,3)
    imagesc(un_bar_pos,1:flen,reshape(get_k_mat(glm_fit(3,cc)),flen,n_bar_pos))
    ca = max(abs(caxis())); caxis([-ca ca]*0.9);
    title(sprintf('Ori: 90, Mag: %.2d',ca));
    
    subplot(2,2,4)
    imagesc(un_bar_pos,1:flen,reshape(get_k_mat(glm_fit(4,cc)),flen,n_bar_pos))
    ca = max(abs(caxis())); caxis([-ca ca]*0.9);
    title(sprintf('Ori: 135, Mag: %.2d',ca));
    
    pause
    clf
end
clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
%%
stim_fs = 100; %in Hz
Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[2 40]/niqf);
use_lfps = [1:4:96];

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
flen = 30;
beg_buffer = round(stim_fs*0.15);
bar_oris = [0 45 90 135];
% for EE = 1:length(bar_oris);
for EE = 1
    fprintf('Analyzing %d ori expts\n',bar_oris(EE));
    
    %     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE));
    
    %expts with X deg bars and gray back (sim sacs)
    cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris(EE));
    cur_expt_set(cur_expt_set==2) = []; %this expt had continuous bar luminance
    cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
    
    un_bar_pos = all_un_bar_pos(:,EE);
    
    %%
    all_stim_times = [];
    all_rel_stimes = [];
    all_rel_etimes = [];
    all_phase = [];
    all_Op = [];
    all_bar_mat = [];
    all_used_inds = [];
    %     all_binned_spks = [];
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
            
            %                 cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(tt).End(end)/1e4];
            %                 cur_binned_spks = nan(length(cur_stim_times),length(use_sus));
            %                 for cc = 1:length(use_sus)
            %                     cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges);
            %                     cur_binned_spks(:,cc) = cur_hist(1:end-1);
            %                 end
            
            cur_bar_mat = zeros(length(cur_stim_times),n_bar_pos);
            for bb = 1:n_bar_pos
                cur_set = find(cur_Op==un_bar_pos(bb));
                pset = cur_phase(cur_set) == 0;
                nset = cur_phase(cur_set) == pi;
                cur_bar_mat(cur_set(pset),bb) = 1;
                cur_bar_mat(cur_set(nset),bb) = -1;
                %                 cur_bar_mat(cur_set,bb) = 1;
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
            %                 all_binned_spks = [all_binned_spks; cur_binned_spks];
            all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
        end
        
    end
    
    %%
    all_Vmat_interp = [];
    for ee = 1:length(cur_expt_set);
        cur_expt = cur_expt_set(ee);
        
        Vmat = [];
        for ll = 1:length(use_lfps)
            fprintf('Electrode %d of %d\n',ll,length(use_lfps));
            filename = sprintf('Expt%d.p%dFullV.mat',cur_expt,use_lfps(ll));
            load(filename);
            V = double(FullV.V);
            %     V = V + FullV.sumscale*sumv;
            V = V*FullV.intscale(1)/FullV.intscale(2);
            nparts = length(FullV.blklen);
            dV = [];
            %splice together multiple blocks
            cur_pt = 1;
            for pp = 1:nparts
                cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
                cur_range(cur_range > length(V)) = [];
                curV = decimate(V(cur_range),dsf);
                curV = filtfilt(filt_b,filt_a,curV);
                dV = [dV curV];
                cur_pt = cur_pt + FullV.blklen(pp);
            end
            Vmat(:,ll) = dV;
        end
        
        t_ax = [];
        for pp = 1:nparts
            cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
            t_ax = [t_ax downsample(cur_t_ax,dsf)];
        end
        t_ax(size(Vmat,1)+1:end) = [];
        
        fprintf('LFP len: %d\n',range(t_ax));
        
        use_set = find(all_exptvec==ee);
        Vmat_interp = interp1(t_ax,Vmat,all_stim_times(use_set));
        all_Vmat_interp = [all_Vmat_interp; Vmat_interp];
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
    Yobs = zscore(all_Vmat_interp);
    
    %%
    stim_params.spatial_dims = 1;
    stim_params.sdim = n_bar_pos;
    SDIM = n_bar_pos;
    stim_params.flen = flen;
    klen = length(un_bar_pos)*flen;
    clear defmod
    
    for cc = 1:length(use_lfps)
        fprintf('Cell %d of %d\n',cc,length(use_lfps));
                
        defmod.lambda_L1x = 0;
        defmod.lambda_d2T = 20;
        defmod.lambda_d2X = 20;
        kern_types{1} = 'lin';
        init_kerns = 0.01*randn(klen,1);
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        glm = createGNM(init_kerns,1,kern_types,defmod,stim_params,'gauss');
        glm_fit(EE,cc) = fitGNM_filters(glm,abs(Xmat(tr_inds,:)),Yobs(tr_inds,cc),'none',[],1e-4,1e-6,0);
        
%         defmod.lambda_L1x = 0;
%         defmod.lambda_d2T = 20;
%         npq = 2;
%         nnq = 2;
%         kern_types{1} = 'lin';
%         for i = 2:(1+npq+nnq)
%             kern_types{i} = 'quad';
%         end
%         init_kerns = 0.01*randn(klen,1+npq+nnq);
%         init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
%         quad = createGNM(init_kerns,[1 ones(1,npq) -1*ones(1,nnq)],kern_types,defmod,stim_params,'gauss');
%         for i = 2:(1+npq+nnq)
%             quad.mods(i).lambda_L1x = 0;
%             quad.mods(i).lambda_d2T = 10;
%         end
%         quad_fit(EE,cc) = fitGNM_filters(quad,Xmat(tr_inds,:),Yobs(tr_inds,cc),'none',[],1e-4,1e-6,0);
        
        [glm_sse, ~, ~, prate] = getLL_GNM(glm_fit(EE,cc),abs(Xmat(xv_inds,:)),Yobs(xv_inds,cc),'none');
        
%         [quad_sse, ~, ~, prate] = getLL_GNM(quad_fit(EE,cc),Xmat(xv_inds,:),Yobs(xv_inds,cc),'none');
        
        yobs_var = var(Yobs(xv_inds,cc));
        glm_r2(EE,cc) = 1-glm_sse/yobs_var;
%         quad_r2(EE,cc) = 1-quad_sse/yobs_var;
        
    end
    
end

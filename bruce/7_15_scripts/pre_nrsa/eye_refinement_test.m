
diff_used_inds = diff(used_inds);
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);
x_drifts_l = nan(size(avg_x_pos));
y_drifts_l = nan(size(avg_y_pos));
x_drifts_r = nan(size(avg_x_pos));
y_drifts_r = nan(size(avg_y_pos));
for i = 1:n_fixs
    cur_set = rel_fix_start_inds(i):rel_fix_stop_inds(i);
    medx_pos = median(full_eyepos(cur_set,1));
    medy_pos = median(full_eyepos(cur_set,2));
    x_drifts_l(cur_set) = full_eyepos(cur_set,1) - medx_pos;
    y_drifts_l(cur_set) = full_eyepos(cur_set,2) - medy_pos;
    medx_pos = median(full_eyepos(cur_set,3));
    medy_pos = median(full_eyepos(cur_set,4));
    x_drifts_r(cur_set) = full_eyepos(cur_set,3) - medx_pos;
    y_drifts_r(cur_set) = full_eyepos(cur_set,4) - medy_pos;
end

%%
max_shift = 16;
dshift = 1;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

fix_durs = rel_fix_stop_inds - rel_fix_start_inds;
usable_fixs = find(fix_durs >= 80);

ov_eps = squeeze(eps(2,:,:));

%%
chunk_dur = 5;
n_chunks = 80/chunk_dur;

eps_prior_sigma = 0.25; %0.2
leps_prior = zeros(n_shifts,1);
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
% leps_prior = ones(n_shifts,1);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior));
% cur_cent = [ov_eps(cur_fix,1) ov_eps(cur_fix,2)];
% cdist = sum((bsxfun(@minus,SH,cur_cent)/Fsd).^2,2);
% leps_prior = -cdist/(2*eps_prior_sigma^2);
% leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior));

deps_sigma = 0.075; %0.06
cdist = squareform(pdist(SH/Fsd));
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@plus,lA,leps_prior');
% lA = repmat(leps_prior',n_shifts,1);
lA = bsxfun(@minus,lA,logsumexp(lA,2));

for cf = 1:length(usable_fixs)
    fprintf('Fixation %d of %d\n',cf,length(usable_fixs));
    cur_fix = usable_fixs(cf);
    cur_im_nums = rel_fix_start_inds(cur_fix):rel_fix_stop_inds(cur_fix);
    LL_fun = nan(length(cur_im_nums),n_shifts);
    
    shifted_ims = nan(length(x_shifts)*length(y_shifts),sdim^2);
    for ii = 1:length(cur_im_nums)
        Robs = full_binned_spks(cur_im_nums(ii),:)';
        cur_im = reshape(fullX(cur_im_nums(ii),:),sdim,sdim);
        
        cnt = 1;
        for xx = 1:length(x_shifts)
            for yy = 1:length(y_shifts)
                d2 = dist_shift2d(cur_im, x_shifts(xx),2,0);
                d2 = dist_shift2d(d2,y_shifts(yy),1,0);
                shifted_ims(cnt,:) = d2(:);
                cnt = cnt + 1;
            end
        end
        
        gabor_outs1 = gabor_emp1_filt_r*shifted_ims';
        gabor_outs2 = gabor_emp2_filt_r*shifted_ims';
        gabor_outs = gabor_outs1.^2 + gabor_outs2.^2;
        gabor_outs = sqrt(gabor_outs);
        pred_Rs = bsxfun(@times,gabor_outs,gabor_params_r(:,7));
        pred_Rs = bsxfun(@plus,pred_Rs,gabor_params_r(:,8));
        pred_Rs = log(1+exp(pred_Rs));
        LLs = bsxfun(@times,log(pred_Rs),Robs);
        LLs = LLs - pred_Rs;
        LL_fun(ii,:) = sum(LLs);
    end
    %
    %     n_chunks = floor(length(cur_im_nums)/chunk_dur);
    chunk_starts = (0:n_chunks-1)*chunk_dur + 1;
    chunk_stops = chunk_starts + chunk_dur;
    lB = nan(n_chunks,n_shifts);
    for i = 1:n_chunks
        lB(i,:) = sum(LL_fun(chunk_starts(i):chunk_stops(i),:));
    end
    
    lalpha=zeros(n_chunks,n_shifts);
    lbeta = zeros(n_chunks,n_shifts);
    lscale=zeros(n_chunks,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior' + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:n_chunks
        lalpha(t,:) = logmulexp(lalpha(t-1,:),lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(n_chunks,:)=log(ones(1,n_shifts)) - lscale(n_chunks);
    for t=n_chunks-1:-1:1
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,lA') - lscale(t);
    end
    
    %compute posteriors over hidden states
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    max_loc = nan(n_chunks,1);
    for i = 1:n_chunks
        cur_set = (i-1)*chunk_dur + (1:chunk_dur);
        [~,max_loc(i)] = max(lgamma(i,:));
    end
      
    est_corrections(cf,:,:) = SH(max_loc,:);
    
end
orig_corrections = [ov_eps(usable_fixs,1) ov_eps(usable_fixs,2)];
%%
xxx = round(chunk_dur*(1:n_chunks)-chunk_dur/2);
clear meas_*
for cf = 1:length(usable_fixs)
    fprintf('Fixation %d of %d\n',cf,length(usable_fixs));
    cur_fix = usable_fixs(cf);
    cur_im_nums = rel_fix_start_inds(cur_fix):rel_fix_stop_inds(cur_fix);
    cur_im_nums(81:end) = [];
    meas_x_l(cf,:) = x_drifts_l(cur_im_nums(xxx));
    meas_x_r(cf,:) = x_drifts_r(cur_im_nums(xxx));
    meas_y_l(cf,:) = y_drifts_l(cur_im_nums(xxx));
    meas_y_r(cf,:) = y_drifts_r(cur_im_nums(xxx));
    
    xlcorr(cf) = corr(meas_x_l(cf,:)',squeeze(est_corrections(cf,:,1))');
    xrcorr(cf) = corr(meas_x_r(cf,:)',squeeze(est_corrections(cf,:,1))');
    ylcorr(cf) = corr(meas_y_l(cf,:)',squeeze(est_corrections(cf,:,2))');
    yrcorr(cf) = corr(meas_y_r(cf,:)',squeeze(est_corrections(cf,:,2))');
end



%%
% for i = 1:8
%     subplot(4,2,i)
%     imagesc(x_shifts/Fsd,y_shifts/Fsd,reshape(lgamma(i,:),33,33)); colorbar
% end
xxx = chunk_dur*(1:n_chunks)-chunk_dur/2;
figure
plot(x_drifts_l(cur_im_nums));hold on
plot(y_drifts_l(cur_im_nums),'r')
plot(x_drifts_r(cur_im_nums),'--');hold on
plot(y_drifts_r(cur_im_nums),'r--')
plot(xxx,SH(max_loc,1)/Fsd,'o-')
plot(xxx,SH(max_loc,2)/Fsd,'ro-')
xl = xlim();
line(xl,[ov_eps(cur_fix,1) ov_eps(cur_fix,1)]/Fsd)
line(xl,[ov_eps(cur_fix,2) ov_eps(cur_fix,2)]/Fsd,'color','r')



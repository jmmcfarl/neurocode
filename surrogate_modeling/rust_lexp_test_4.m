cd ~/James_scripts/surrogate_modeling/
addpath(genpath('~/James_scripts'));
clear all
close all

dt = 0.01; %in s
max_rate = Inf; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

NT = 100000; SDIM = 16; flen = 14;
stim = round(2*rand(NT,SDIM))-1;
xvNT = 30000; 
xvstim = round(2*rand(xvNT,SDIM))-1;

% GENERATE A FILTER
nfilts = 8;
x = repmat(1:SDIM,flen,1);
sigma1 = repmat(1.5,flen,SDIM);
sigma2 = repmat(1.5,flen,SDIM);
lambda = repmat(4.5,flen,SDIM);
% amp_vec = gampdf((1:flen)-1,9,0.5);
% amp_vec = gampdf((1:flen)-1,8,0.90);
amp_vec = ncx2pdf((1:flen-2)-1,4,1);
amp_vec = [0 0 amp_vec];
amp_vec = fliplr(amp_vec);
b = repmat(amp_vec',1,SDIM);
a = repmat(0,flen,SDIM);
psi1 = repmat(linspace(0,5*pi,flen)',1,SDIM);
psi2 = flipud(psi1);
psi1(1:5,:) = repmat(psi1(6,:),5,1);
psi2(1:5,:) = repmat(psi2(6,:),5,1);

% xovals = linspace(4,13,nfilts);
xovals = linspace(3,15,nfilts);
% xovals = rand(1,nfilts)*11+3;
for i = 1:nfilts
    x0 = repmat(xovals(i),flen,SDIM);
%     if mod(i,2)==1
        cur_psi = mod(psi1 + rand*2*pi,2*pi);
        temp = b.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
%     else
%         cur_psi = mod(psi2 + rand*2*pi,2*pi);
%         temp = b.*exp(-((x-x0).^2./2./sigma2.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi2))+a;
%     end
    filt(i,:,:) = temp/norm(temp(:));
    filt_mat(i,:) = temp(:)/norm(temp(:));
end
filt_mat = filt_mat';
cvals = xovals;
cvals = (cvals-2);
% cvals(7:end) = fliplr(cvals(1:6));
cvals(ceil(nfilts/2)+1:end) = fliplr(cvals(1:floor(nfilts/2)));

% CREATE TIME-EMBEDDED STIMULUS
stim_emb = makeStimRows(stim,flen);
xvstim_emb = makeStimRows(xvstim,flen);

% FILTER STIMULUS
for i = 1:nfilts
    temp = filt(i,:,:);
    filt_stim(i,:) = zscore(stim_emb*temp(:));
    xvfilt_stim(i,:) = zscore(xvstim_emb*temp(:));
end
%% create spike function
%pass through internal NLs
beta = 4;
% coefs = [1 -1 1 -1 1 -1 1 -1];
% coefs = [1 -1 1 -1 1 -1 1 -1].*cvals;
% coefs = [1 -0.75 1 -0.75 1 -0.75 1 -0.75 1 -0.75].*cvals;
% coefs = [1 -0.75 1 -0.75 1 -0.75 1 -0.75].*cvals;
% coefs = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1].*cvals;
coefs = ones(1,nfilts).*cvals;
% coefs = ones(1,nfilts);
% coefs(mod(1:nfilts,2)==0) = -coefs(mod(1:nfilts,2)==0);
% coefs(coefs < 0) = coefs(coefs < 0)*0.75;
% coefs = [1 -1.5 1 -1.5]/2;
% coefs = [2 -1]/2;
g = zeros(1,NT);
xvg = zeros(1,xvNT);
for i = 1:nfilts
%     lfilt_stim(i,:) = 1/beta*log(1+exp(beta*filt_stim(i,:)));
%     xvlfilt_stim(i,:) = 1/beta*log(1+exp(beta*xvfilt_stim(i,:)));
    lfilt_stim(i,:) = filt_stim(i,:);
    lfilt_stim(i,filt_stim(i,:) < 0) = 0;   
    xvlfilt_stim(i,:) = xvfilt_stim(i,:);
    xvlfilt_stim(i,xvfilt_stim(i,:) < 0) = 0;
%     lfilt_stim(i,:) = filt_stim(i,:).^2;
%     xvlfilt_stim(i,:) = xvfilt_stim(i,:).^2;
    g = g + coefs(i)*lfilt_stim(i,:);
    xvg = xvg + coefs(i)*xvlfilt_stim(i,:);
end
g = g/std(g);

target_rate = 50; %in Hz
target_pspike = target_rate*dt;
cur_theta = 3;
cur_beta = 4;
p_spike = 1/cur_beta*log(1+exp(cur_beta*(g-cur_theta)));
xvp_spike = log(1+exp((xvg-cur_theta)));
scale_f = target_rate/(mean(p_spike)/dt);
p_spike = p_spike*scale_f;
xvp_spike = xvp_spike*scale_f;
% 
% figure
% [y,x] = ksdensity(g);
% plot(x,y); hold on
% yl = ylim();
% tx = -10:.02:10;
% plot(tx,scale_f/cur_beta*log(1+exp(cur_beta*(tx-cur_theta))),'k')
% ylim(yl)

%%
spikes = poissrnd(p_spike);
rbins = (find(spikes>0.5));
nsp = spikes(rbins);
spk_vals = unique(spikes); spk_vals(spk_vals==0) = [];
spikebins = [];
for i = 1:length(spk_vals)
    cur_set = find(spikes == spk_vals(i));
    spikebins = [spikebins; repmat(cur_set(:),spk_vals(i),1)];
end
% spikebins = unique(spikebins); %max 1 spike per bin
fprintf('Nspks: %d\n',length(spikebins));

xvspikes = poissrnd(xvp_spike);
xvrbins = (find(xvspikes>0.5));
xvnsp = spikes(xvrbins);
xvspk_vals = unique(xvspikes); xvspk_vals(xvspk_vals==0) = [];
xvspikebins = [];
for i = 1:length(xvspk_vals)
    cur_set = find(xvspikes == xvspk_vals(i));
    xvspikebins = [xvspikebins; repmat(cur_set(:),xvspk_vals(i),1)];
end
% spikebins = unique(spikebins); %max 1 spike per bin
% fprintf('Nspks: %d\n',length(xvspikebins));

%
nneg =5;
npos = 5;
spike_cond_stim = stim_emb(spikebins,:);
sta      = mean(spike_cond_stim) - mean(stim_emb);
sta = sta/norm(sta);
proj_mat = sta'/(sta*sta')*sta;
stim_proj = stim_emb - stim_emb*proj_mat;
stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
figure
subplot(3,5,1)
imagesc(reshape(sta,flen,SDIM));
for i = 1:5
    subplot(3,5,5+i)
    imagesc(reshape(stcs(:,i),flen,SDIM));
end
for i = 1:5
    subplot(3,5,10+i)
    imagesc(reshape(stcs(:,5+i),flen,SDIM));
    colormap(gray)
end

% stc_dims = [sta' stcs(:,1) stcs(:,5)];
stc_dims = [sta' stcs];

%% FIND SUBSPACE USING LEXP MODEL
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
% sdim = 12; flen = 10;
sdim = 16; flen = 14;
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 0;
defmod.lambda_L1x = 0; %2
defmod.lambda_dT = 0;
defmod.lambda_L2x = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

nmods = 9;
cur_basis = stc_dims;
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);

%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
init_signs = [1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];
glm_stcb = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
% glm_stcb.mods(1) = [];
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_quad2 = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots');
xvLL_quad = getLLGLM_lexp(glm_quad,xvstim_emb,xvspikebins,'none');

%%
nmods = 8;
cur_basis = stc_dims;
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);
%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
% init_signs = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];
init_signs = ones(1,nmods);
init_betas = 2*ones(nmods,1);
init_thetas = zeros(nmods,1);
% glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,[],[]);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
test = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots',[],1e-4,1e-6);
[test,coms,peak_locs] = get_filter_coms_1d(test);
[~,ord] = sort(peak_locs);
test.mods = test.mods(ord);
%%
temp = test;
for i = 1:length(temp.mods)
    temp.mods(i).nltype = 'uncon';
    temp.mods(i).lnl = 0;
    temp.mods(i).lnl2 = 20;
    temp.mods(i).nlcon = 0;
    temp.mods(i).nlmon = 1;
end
k_mat = get_k_mat(temp);
input = stim_emb*k_mat;
temp = normalizeRFs_full(temp,stim_emb); %This normalizes filter coefs so that the filtered stimulus distribution is standard normal (not really necessary)
temp = fitWeights_full(temp,input,spikebins,1);
fo = fitNL_full(temp,input,spikebins,1);
fo = fitWeights_full(fo,input,spikebins,1);
fo = fitNL_full(fo,input,spikebins,1);
fo = fitWeights_full(fo,input,spikebins,1);

for i = 1:6
subplot(3,2,i)
plot(test.mods(i).nlx,test.mods(i).nly/max(test.mods(i).nly))
hold on
plot(fo.mods(i).nlx,fo.mods(i).nly/max(fo.mods(i).nly),'r')
xlim([-3.5 3.5])
ylim([0 1])
end
%%
figure
for i = 1:8
    subplot(4,2,i)
    imagesc(squeeze(filt(i,:,:)*abs(coefs(i))))
    set(gca,'ydir','normal')
    colormap(gray)
%     caxis([-0.4 0.2])
    caxis([-0.5 0.3]*8)
end


test = get_filter_coms_1d(test);
coms = arrayfun(@(x) x.filt_com,test.mods);
w = arrayfun(@(x) x.w,test.mods);
efilts = find(w > 0);
ifilts = find(w < 0);
[~,eord] = sort(coms(efilts));
eord = efilts(eord);
[~,iord] = sort(coms(ifilts));
iord = ifilts(iord);
k_mat = get_k_mat(test);
figure
for i = 1:4
    subplot(4,2,(i-1)*2+1)
    imagesc(reshape(k_mat(:,eord(i)),flen,sdim))
    set(gca,'ydir','normal')
    colormap(gray)
    caxis([-2 1.25])
end
for i = 1:4
    subplot(4,2,(i-1)*2+2)
    imagesc(reshape(k_mat(:,iord(i)),flen,sdim))
    set(gca,'ydir','normal')
    colormap(gray)
    caxis([-2 1.25])
end

test = get_filter_coms_1d(glm_quad2);
coms = arrayfun(@(x) x.filt_com,test.mods);
w = arrayfun(@(x) x.w,test.mods);
efilts = find(w > 0);
ifilts = find(w < 0);
k_mat = get_k_mat(test);
k_norms = sqrt(sum(k_mat.^2));
[~,eord] = sort(k_norms(efilts),'descend');
eord = efilts(eord);
[~,iord] = sort(k_norms(ifilts),'descend');
iord = ifilts(iord);
figure
for i = 1:4
    subplot(4,2,(i-1)*2+1)
    imagesc(reshape(k_mat(:,eord(i)),flen,sdim))
    set(gca,'ydir','normal')
    colormap(gray)
    caxis([-0.6 0.6])
end
for i = 1:4
    subplot(4,2,(i-1)*2+2)
    imagesc(reshape(k_mat(:,iord(i)),flen,sdim))
    set(gca,'ydir','normal')
    colormap(gray)
    caxis([-0.6 0.6])
end

%%


nmods = 1; 
cur_basis = stc_dims;
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);
init_signs = 1;
init_kerns = cur_basis*STCcf_0;
glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_tlin(1) = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots');
xvLL_tlin(1) = getLLGLM_lexp(glm_tlin(1),xvstim_emb,xvspikebins,'none');

for i = 2:10
    fprintf('Mod %d\n',i);
    if mod(i,2)==0
        cur_weight = -1;
    else
        cur_weight = 1;
    end
   glm0 = add_null_filter(glm_tlin(i-1),'zero',cur_weight,'threshlin'); 
   glm_tlin(i) = fitGLM_lexp(glm0,stim_emb,spikebins,'tots');
   xvLL_tlin(i) = getLLGLM_lexp(glm_tlin(i),xvstim_emb,xvspikebins,'none');
end

% %normalize
% for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
% init_kerns = cur_basis*STCcf_0;
% init_signs = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];
% init_betas = 4*ones(nmods,1);
% init_thetas = zeros(nmods,1);
% 
% % glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,[],[]);
% glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
% glm_stcb.image_type = '1d';
% glm_stcb.spk_nl = 'logexp';
% glm_tlin = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots');
% xvLL_tlin = getLLGLM_lexp(glm_tlin,xvstim_emb,xvspikebins,'none');

%%
n_it = 3;
for i = 1:n_it
    glm_lexp2 = update_ARD_priors(glm_lexp2);
    l2pens(i,:) = arrayfun(@(x) x.lambda_L2x,glm_lexp2.mods);
    glm_lexp2 = fitGLM_lexp(glm_lexp2,stim_emb,spikebins,'tots');
end

%%
used_mods = [6:10];
for i = 1:length(used_mods)
   temp_glm = glm_tlin;
   temp_glm.mods(used_mods(i)+1:end) = [];
   mod_fit(i) = fitGLM_lexp(temp_glm,stim_emb,spikebins,'tots');
   mod_xvLL(i) = getLLGLM_lexp(mod_fit(i),stim_emb,spikebins,'none');
end


%%
for i = 1:nmods
    glm_lexp.mods(i).nltype = 'uncon';
    glm_lexp.mods(i).nlmon = 1;
end
full_fit = fitNLHI2d_fullbf(glm_lexp,stim_emb,spikebins,'tots');

%% SUBSPACE OVERLAP ANALYSIS
lexp_k_mat = get_k_mat(glm_tlin);
% lexp2_k_mat = get_k_mat(glm_lexp2);
quad_k_mat = get_k_mat(glm_quad);
rquad_k_mat = get_k_mat(glm_rquad);


overlap_rquad = subspace_overlap(rquad_k_mat,filt_mat);
overlap_lexp = subspace_overlap(lexp_k_mat,filt_mat);
% overlap_lexp2 = subspace_overlap(lexp2_k_mat,filt_mat);
overlap_quad = subspace_overlap(quad_k_mat,filt_mat);
overlap_stc = subspace_overlap(stc_dims,filt_mat);

soverlap_lexp = subspace(filt_mat,lexp_k_mat);
% soverlap_lexp2 = subspace(filt_mat,lexp2_k_mat);
soverlap_quad = subspace(filt_mat,quad_k_mat);
soverlap_stc = subspace(filt_mat,stc_dims);

%%
nmods = 16;
cur_basis = get_k_mat(glm_rquad);
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);
%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
init_signs = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];
% init_signs = ones(1,nmods);
init_betas = 2*ones(1,nmods);

defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 0;
defmod.lambda_L1x = 0; %2
defmod.lambda_dT = 0;
defmod.SDIM = SDIM;
defmod.fsdim = SDIM;
defmod.pids = 1:SDIM;
basis = 'pix';

glm_stcb = createGLM_lexp_rot(cur_basis,STCcf_0,init_signs,init_betas,defmod,basis);
glm_stcb.image_type = '1d';
for i = 1:length(glm_stcb.mods)
    glm_stcb.mods(i).nltype = 'rquad';
end
glm_stcb.spk_nl = 'logexp';
% glm_stcb.spk_nl = 'exp';pl
[glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,stim_emb*cur_basis);
kern_output = stim_emb*cur_basis;
glm_lexp_rot = fitGLM_lexp_rot(glm_stcb,kern_output,spikebins,'tots',100);

%%
for i = 1:length(glm_lexp_rot.mods)
    glm_lexp_rot.mods(i).locLambda = 20;
end
glm_lexp_rot2 = fitGLM_lexp_rot(glm_lexp_rot,kern_output,spikebins,'tots');


%%
% 
% %%
% nmods = 4;
% % cur_basis = get_k_mat(glm_lexp2);
% cur_basis = stc_dims;
% n_bvs = size(cur_basis,2);
% STCcf_0 = randn(n_bvs,nmods);
% %normalize
% for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
% init_kerns = cur_basis*STCcf_0;
% init_signs = [1 -1 1 -1];
% % init_signs = ones(1,nmods);
% init_betas = 4*ones(1,nmods);
% 
% defmod.lnl = 0;
% defmod.lh = 0;
% defmod.lnl2 = 0;
% defmod.lh2 = 0;
% defmod.nlcon = 0;
% defmod.nlmon = 0;
% defmod.locLambda = 0;
% defmod.lambda_dX = 0;
% defmod.lambda_L1x = 0; %2
% defmod.lambda_dT = 0;
% defmod.locLambda = 0;
% defmod.SDIM = sdim;
% defmod.fsdim = sdim;
% defmod.pids = 1:sdim;
% basis = 'pix';
% 
% init_nls{1} = 'pquad';
% init_nls{2} = 'nquad';
% init_nls{3} = 'pquad';
% init_nls{4} = 'nquad';
% for i = 1:nmods
%     nltypes{i} = 'quad';
% end
% Robs = zeros(1,NT);
% ftable = tabulate(spikebins);
% Robs(ftable(:,1)) = ftable(:,2);
% 
% glm_stcb = createGLM0_stcb_connl(cur_basis,STCcf_0,defmod,init_nls,nltypes,'test')
% glm_stcb.image_type = '1d';
% glm_stcb.spk_nl = 'logexp';
% for i = 1:nmods
%     glm_stcb.mods(i).w = init_signs(i);
% end
% kern_output = stim_emb*cur_basis;
% glm_lexp = fitSTCBF_connl(glm_stcb,kern_output,Robs,0);

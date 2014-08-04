function [mod_seq,LL_seq,kscale_seq] = iterateGNM_filts_NLs(mod0,X,spikebins,max_iter)

if nargin < 4
    max_iter = 10;
end

mod_seq = [mod0];
nfilts = length(mod0.mods);

[LL_0, ~, ~, ~, ~,ig] = getLL_GNM(mod0,X,spikebins,'none');
out_scales = zeros(nfilts,1);
for i = 1:nfilts
    out_scales(i) = std(ig{i});
end
cur_out_scales = out_scales;

k_scales = sqrt(sum(get_k_mat(mod0).^2));

LL_seq = [LL_0];
kscale_seq = [k_scales];
max_dLL = 1e-3;
dLL = Inf;
curLL = LL_0;
curmod = mod0;
it_cnt = 1;
while dLL > max_dLL
    fprintf('Iteration %d\n',it_cnt);
    curmod = fitGNM_internal_NLs(curmod,X,spikebins,1,0);
    [~, ~, ~, ~, ~,ig] = getLL_GNM(curmod,X,spikebins,'none');
    for i = 1:nfilts
        cur_out_scales(i) = std(ig{i});
        curmod.mods(i).nly = curmod.mods(i).nly/(cur_out_scales(i)/out_scales(i));
    end
        
    curmod = fitGNM_filters(curmod,X,spikebins,'none',[],1e-4,1e-6);
    [newLL] = getLL_GNM(curmod,X,spikebins,'none');
    k_scales = sqrt(sum(get_k_mat(curmod).^2));
    kscale_seq = [kscale_seq; k_scales];
    LL_seq = [LL_seq newLL];
    mod_seq = [mod_seq curmod];
    dLL = curLL-newLL;
    curLL = newLL;
    fprintf('LL improvement %.3f\n',dLL);
    
    it_cnt = it_cnt + 1;
end

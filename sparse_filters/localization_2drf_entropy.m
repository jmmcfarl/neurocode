function L = localization_2drf_entropy(T,loading_mat)

[N_k,N_f] = size(loading_mat);
flen = 9;
npix = 32;
spixids = 1:32^2;

L = loading_mat*inv(T');

norm_L = L;
for i = 1:N_f
    norm_L(:,i) = l2norm(L(:,i));
end

% pow_profiles = get2dPowerProfs(norm_L,flen,npix,spixids); %get spatial profiles of temporal power
pow_profiles = zeros(length(spixids),N_f);
for itk =1:N_f; 
	tkmat     = reshape(norm_L(:,itk),flen,length(spixids));
	pow_profiles(:,itk)   = sum(tkmat.^2)';
end

pow_profiles = pow_profiles./repmat(sum(pow_profiles),length(spixids),1);

entropies = nan(N_f,1);
for i = 1:N_f
    entropies(i) = sum(pow_profiles(:,i).*log(pow_profiles(:,i)));
end
L = -sum(entropies)-log(det(T*T'));
% L = -sum(entropies)/sqrt(det(T*T'));
function stim = generate_RLS_stim(NT,npix,dds,usfac,seed)
%
% makes a random bar stimulus with NT frames, npix pixels, and with a
% sparsity of dds. can incorporate spatial up-sampling by factor usfac
%

if nargin < 4 || isempty(usfac)
    usfac = 1;
end
if nargin < 5 || isempty(seed)
    seed = nan;
end
if ~isnan(seed)
    rng(seed);
end

stim = randi(2,NT,npix);
stim(stim == 2) = -1;

to_zero = rand(NT,npix);
stim(to_zero > dds/100) = 0;

if usfac > 1
    stim_up = zeros(NT, npix*usfac);
    for ii = 1:npix
        for jj = 1:usfac
            stim_up(:,usfac*(ii-1)+jj) = stim(:,ii);
        end
    end
    stim = stim_up;
end

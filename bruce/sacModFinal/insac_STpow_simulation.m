NT = 1e3;
wi = 2;
bar_width = 0.0565;
usfac = 10;
nbars = round(wi/bar_width);

dd = 12;

is_zero = rand(NT,nbars) > dd/100;
stim = randi(2,NT,nbars);
stim(stim==2) = -1;
stim(is_zero) = 0;


%%
npix = nbars*usfac;
pix_dx = bar_width/usfac;
if usfac > 1
    stim_up = zeros(NT,npix);
    for ii = 1:nbars
        for jj = 1:usfac
            stim_up(:,usfac*(ii-1)+jj) = stim(:,ii);
        end
    end
elseif usfac == 1
    stim_up = all_stim_mat;
end

%%
poss_smooth_widths = 1:1.5*usfac;
for pp = 1:length(poss_smooth_widths)
    smooth_width = poss_smooth_widths(pp);
stim_up_conv = stim_up;
for ii = 1:smooth_width
    uset = 1:(npix-ii);
    uset_map = (ii+1):npix;
    stim_up_conv(:,uset) = stim_up_conv(:,uset) + stim_up(:,uset_map);
end
stim_up_conv = stim_up_conv/(smooth_width + 1);

niqf_x = 1/(2*pix_dx);
niqf_t = 1/(2*dt);
fx = linspace(-niqf_x,niqf_x,npix);
ft = linspace(-niqf_t,niqf_t,NT);
[FX,FT] = meshgrid(fx,ft);

dft = 2*niqf_t/NT;
dfx = 2*niqf_x/npix;
PP = abs(fftshift(fft2(stim_up)));
PP_conv = abs(fftshift(fft2(stim_up_conv)));

H = fspecial('gaussian',[50 50],3);
PP_filt = filter2(H,PP);
PP_conv_filt = filter2(H,PP_conv);

avg_PP = mean(PP);
avg_PP_conv = mean(PP_conv);
dom(pp,:) = -(avg_PP_conv - avg_PP)./avg_PP;

end
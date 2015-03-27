% NT = 1e3;
% wi = 2;
% bar_width = 0.0565;
% usfac = 10;
% nbars = round(wi/bar_width);
% 
% dd = 12;
% 
% 
% is_zero = rand(NT,nbars) > dd/100;
% stim = randi(2,NT,nbars);
% stim(stim==2) = -1;
% stim(is_zero) = 0;
% 
% 
% %%
% npix = nbars*usfac;
% pix_dx = bar_width/usfac;
% if usfac > 1
%     stim_up = zeros(NT,npix);
%     for ii = 1:nbars
%         for jj = 1:usfac
%             stim_up(:,usfac*(ii-1)+jj) = stim(:,ii);
%         end
%     end
% elseif usfac == 1
%     stim_up = all_stim_mat;
% end
% 
% %%
% poss_smooth_widths = 1:1.5*usfac;
% for pp = 1:length(poss_smooth_widths)
%     smooth_width = poss_smooth_widths(pp);
% stim_up_conv = stim_up;
% for ii = 1:smooth_width
%     uset = 1:(npix-ii);
%     uset_map = (ii+1):npix;
%     stim_up_conv(:,uset) = stim_up_conv(:,uset) + stim_up(:,uset_map);
% end
% stim_up_conv = stim_up_conv/(smooth_width + 1);
% 
% dt = 0.01;
% niqf_x = 1/(2*pix_dx);
% niqf_t = 1/(2*dt);
% fx = linspace(-niqf_x,niqf_x,npix);
% ft = linspace(-niqf_t,niqf_t,NT);
% [FX,FT] = meshgrid(fx,ft);
% 
% dft = 2*niqf_t/NT;
% dfx = 2*niqf_x/npix;
% PP = abs(fftshift(fft2(stim_up)));
% PP_conv = abs(fftshift(fft2(stim_up_conv)));
% 
% H = fspecial('gaussian',[50 50],3);
% PP_filt = filter2(H,PP);
% PP_conv_filt = filter2(H,PP_conv);
% 
% avg_PP = mean(PP);
% avg_PP_conv = mean(PP_conv);
% dom(pp,:) = -(avg_PP_conv - avg_PP)./avg_PP;
% 
% end
% 
% 






%%
clear all
NT = 1e3;
dt = 0.01;
temp_usfac = 10;
up_dt = dt/temp_usfac;
phos_pframes = round(0.004/up_dt);



wi = 2;
bar_width = 0.0565;
% space_usfac = round(bar_width/drift_frame_dist);
space_usfac = 30;
nbars = round(wi/bar_width);
npix = nbars*space_usfac;
pix_dx = bar_width/space_usfac;

poss_drift_rates = (1:10)*pix_dx/up_dt;
for pp = 1:length(poss_drift_rates)
    pp
    drift_rate = poss_drift_rates(pp);
% drift_rate = 0.06/.004;
drift_frame_dist = drift_rate*up_dt;
trans_vel = round(drift_frame_dist/pix_dx);
act_drift_rate(pp) = trans_vel;


dd = 12;


is_zero = rand(NT,nbars) > dd/100;
stim = randi(2,NT,nbars);
stim(stim==2) = -1;
stim(is_zero) = 0;


%%
if space_usfac > 1
    stim_up = zeros(NT,npix);
    for ii = 1:nbars
        for jj = 1:space_usfac
            stim_up(:,space_usfac*(ii-1)+jj) = stim(:,ii);
        end
    end
elseif space_usfac == 1
    stim_up = stim;
end

%%
stim_up = repmat(stim_up,[1 1 temp_usfac]);
stim_up_trans = stim_up;
for ii = 2:phos_pframes
    stim_up_trans(:,:,ii) = shift_matrix_Nd(stim_up(:,:,1),trans_vel*(ii-1),2);
end
stim_up(:,:,(phos_pframes+1):end) = 0;
stim_up_trans(:,:,(phos_pframes+1):end) = 0;

stim_up = permute(stim_up,[3 1 2]);
stim_up = reshape(stim_up,[],npix);

stim_up_trans = permute(stim_up_trans,[3 1 2]);
stim_up_trans = reshape(stim_up_trans,[],npix);

%%
up_NT = NT*temp_usfac;
niqf_x = 1/(2*pix_dx);
niqf_t = 1/(2*up_dt);
fx = linspace(-niqf_x,niqf_x,npix);
ft = linspace(-niqf_t,niqf_t,up_NT);
[FX,FT] = meshgrid(fx,ft);

dft = 2*niqf_t/up_NT;
dfx = 2*niqf_x/npix;
PP = abs(fftshift(fft2(stim_up)));
PP_conv = abs(fftshift(fft2(stim_up_trans)));

[Hx,Ht] = meshgrid(-30:30,-30:30);
Hx = Hx*dfx; Ht = Ht*dft;
sig_x = 0.5; sig_t = 0.5;
H = exp(-(Hx.^2/(2*sig_x^2) + Ht.^2/(2*sig_t^2)));
H = H/sum(H(:));

% H = fspecial('gaussian',[50 50],10);
PP_filt = filter2(H,PP);
PP_conv_filt = filter2(H,PP_conv);

PP_orig{pp} = PP_filt;
PP_change{pp} = (PP_conv_filt - PP_filt)./PP_filt;
all_ft{pp} = ft;
all_fx{pp} = fx;

end

%%
f1 = figure();
for pp = 1:length(poss_drift_rates)
subplot(2,1,1)
imagesc(all_fx{pp},all_ft{pp},PP_change{pp});
caxis([-0.5 0.5])
ylim([-25 25]);
xlim([-15 15]);colorbar
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Temporal freq (Hz)');

subplot(2,1,2)
imagesc(all_fx{pp},all_ft{pp},PP_orig{pp}/max(PP_orig{pp}(:)));
caxis([0 1]);
ylim([-25 25]);
xlim([-15 15]); colorbar
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Temporal freq (Hz)');

pause
clf
end
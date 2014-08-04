function fstim = kfilterInput(stim,k)
%%USAGE: fstim = kfilterInput(stim,k)
% filter input with filter k of ONE specific modul

[NT SDIM]  = size(stim);
rfsiz      = length(k); %dimensionality of RF kernels
rflen      = rfsiz/SDIM; %number of lags in kernel
mk         = flipud(reshape(k,rflen,SDIM)); %matrix kernel (N_tauxD) (time dimension reflected)

%NOTE: need to reflect RF kernel in space before convolving with stim
%matrix!
fstim      = conv2(stim,fliplr(mk),'valid'); %convolve stimulus with appropriately reflected matrix kernel
end




%% and conv2 is even faster! 

% fstim = conv(stim(:,1),mk(:,1),'valid'); 
% for idim = 2:SDIM;  
% 	fstim = fstim + conv(stim(:,idim),mk(:,idim),'valid'); 
% end; 

%% since 2009, conv seems to be faster... 

% x = filter(mk(:,1),1,stim(:,1));
% for i = 2:SDIM;  
% 	x = x + filter(mk(:,i),1,stim(:,i)); 
% end; 
% fstim = x(rflen:end); 

%% older, slower version

%fstim = zeros(NT-rflen,1);
% for i = rflen:NT;
%     fstim(i-rflen+1,1) = reshape(stim((i-rflen+1):i,:),1,rfsiz)*k; 
% end


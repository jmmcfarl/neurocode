function [myglm] = fitNLHI_jmm(myglm,ustim,spkbs,lltol)
% USAGE: [myglm ll0] = fitNLHI(myglm,ustim,spkbs,lltol)
%   iteratively fit nonlinearity and history term
%disp(sprintf('*** Fitting NL and HI of %s ***',myglm.mname'))
%disp('*** Fitting NL and HI ***')



%compute output of stimulus convolved with the internel kernels of each of
%the NL modules
inputs = applykfilters(myglm,ustim);


SDIM   = size(ustim,2); %stimulus dimension

nmax = 10; %max iterations
[ll0, ll0p] = getLLGLM(myglm,ustim,spkbs);

% i    = 0;  
% lldiff=ll0;

% while((lldiff > lltol) && (i < nmax))
%     i = i+1;

nextNLfit = fitNL_jmm(myglm,inputs,spkbs);
nextHIfit = fitHI_jmm(nextNLfit,inputs,spkbs);

%%
lambda_range = linspace(0,1000,40);
[weight_mat,ppls] = scan_weights_lambda(nextHIfit,ustim,spkbs,lambda_range);

subplot(2,1,1)
plot(lambda_range,ppls,'.-')
xlabel('Lambda','fontsize',14)
ylabel('LL','fontsize',14)
subplot(2,1,2)
plot(lambda_range,weight_mat)
xlabel('Lambda','fontsize',14)
ylabel('Weighting','fontsize',14)

%%
nextHIfit.lambdaW = 200;
nextWfit  = fitGLMweights_jmm(nextHIfit,ustim,spkbs);
[llnow, pllnow] = getLLGLM_jmm(nextWfit,ustim,spkbs);

disp(sprintf(' -- NLHI iteration %d -- ll: %0.5f',i,llnow));
drawnow;

lldiff = ll0-llnow;
if(lldiff > 0);
    ll0   = llnow;
    myglm = nextWfit;
end;
% end;
% disp(sprintf(' -- Fitted NL & HI in %d iterations -- ll: %0.5f',i,ll0));
myglm.LL = ll0;
hilen = length(myglm.mods(1).h);
myglm = addFOdens(myglm,ustim,spkbs+hilen-1);
end


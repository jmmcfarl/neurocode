function [theg,x] = fitLFPpowerSpectrum(y,Fl,Fh,Fs)
    %[g,params] = fitLFPpowerSpectrum(y,Fl,Fh,Fs)
    %
    %Fits the PSD of an LFP y with a function:
    %psd = a*log(1+exp(b*(logfreq-c)) + d
    %
    %The function is fit only in the range Fl to Fh. It goes to d at low
    %frequencies and to d+a*b*(logfreq-c)) (log-linear) at high frequencies
    %Fs is the sampling frequency
    %
    % g can then be fed to despikeLFP
    %
    % History:
    %
    % 07/07/2010: Updated first line of documentation and fixed second
    % argument
    % 28/06/2010: v1.0
    %
    % Created by Patrick Mineault
    %
    % See also despikeLFP


    ffty = fft(y);
    allfreqs = (1:length(y)/2)'/length(y)*Fs;
    freqs = allfreqs(allfreqs > Fl & allfreqs < Fh);
    rg = 1+round(freqs*length(y)/Fs);
    logpsdy = log(ffty(rg).*conj(ffty(rg)));
    
    g = @(fs,x) x(1)*log( 1 + exp(x(2)*(log(fs)-x(3)))) + x(4);
    
    %Guestimate the optimal parameters
    smoothpsdy = ksr(log(freqs),logpsdy);
    [gm,maxidx] =  max(diff(smoothpsdy.f));
    c0 = smoothpsdy.x(maxidx+1);
    d0 = max(smoothpsdy.f);
    b0 = 1;
    a0 = (min(smoothpsdy.f)-d0)/(log(Fh)-c0);
    
    %Initialize guesses for parameters
    x0 = [a0;b0;c0;d0];
      
    %The inverse density of points in a spot is ~ diff(log(freqs)), meaning that
    %higher frequencies disproportionally affect the fit on a log scale. 
    %Correct for this.
    idensity = sqrt(diff(log(freqs)));
    idensity = [idensity(1);idensity];
    
    objectivefun = @(x) idensity.*(g(freqs,x) - logpsdy);
    %Feed that into lsqnonlin
    x = lsqnonlin(objectivefun,x0);
    
    allrg = 1+round(allfreqs*length(y)/Fs);
    logpsdy = log(ffty(allrg).*conj(ffty(allrg)));
    smoothpsdy = ksr(log(allfreqs),logpsdy);
    dbfactor = 10/log(10);
    smallfreqs = exp(log(min(allfreqs)):.05:log(max(allfreqs)));
    
    themax = max(dbfactor*smoothpsdy.f);
    
%     figure;
%     semilogx(allfreqs,dbfactor*logpsdy-themax,exp(smoothpsdy.x),dbfactor*smoothpsdy.f-themax,smallfreqs,dbfactor*g(smallfreqs,x)-themax);
%     box off;
%     xlabel('Frequency (Hz)');
%     ylabel('PSD (dB relative to max)');
%     legend('Empirical PSD of wideband signal','Smoothed PSD of wideband signal','Estimated PSD of LFP');
%     legend(gca,'Location','SouthWest');
%     legend boxoff;
%     xlim([min(allfreqs)/3,max(allfreqs)*3]);

    gx = exp(g(allfreqs,x));
    theg = [gx(1);gx;gx(end-1:-1:1)];
end

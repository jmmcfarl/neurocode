function [y] = delineSignal(dat,Fs,freqs,freqrange)
    %Deline a signal by fitting the PSD around the time of peaks in
    %frequency to a function, then filtering the signal with the square root 
    %inverse of that function so that the spectrum will be flat around that
    %frequency
    %
    % dat: the data, a vector
    % Fs: the sampling frequency in Hz
    % freqs: a vector of peak frequencies to equalize, in Hz (default:
    % [60,180]
    % freqrange: the length over which we expect the PSD of dat to be
    % stable, default 4Hz (+- 2Hz)
    if nargin < 3
        freqs = [60,180];
    end
    if nargin < 4
        freqrange = 4;
    end
    
    a = fft(double(dat(:)));

    %Now remove line noise from datlo to obtain y
    %At this point, ran curve fitting tool
    %
    %Eliminate line noise at 60 Hz
    thefilt = ones(size(a));
    
    figure;
    
    idx = 1;
    freqs = freqs(:)';
    for tgtr = freqs
    
        peak = round(tgtr*length(dat)/Fs);

        rg = round(((peak-freqrange/2*length(dat)/Fs):(peak+freqrange/2*length(dat)/Fs)))';
        a60 = abs(a(rg));
        x0 = [max(a60)-median(a60),.05,median(a60),peak,0.4]';
        thefun = @(p,x) p(1)*exp(-p(2)*abs(x-p(4)).^(p(5)))+p(3);
        [ps resnorm] = lsqcurvefit(thefun,x0,rg,a60);
        %plot(x,a60,x,thefun(ps,x))
        thefilt(rg) = ps(3)./thefun(ps,rg);
        thefilt(end-rg-2) = thefilt(rg);
        
        subplot(length(freqs),1,idx);
        plot(rg/length(dat)*Fs,a60,rg/length(dat)*Fs,thefun(ps,rg),'r');
        idx = idx+1;
        title(sprintf('PSD around %d and function fit',round(tgtr)));
    end

    % y is datlo with line noise removed
    y = real(ifft(a.*thefilt));
end
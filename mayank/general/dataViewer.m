%plot large data set in segments
function [] = dataViewer(synct,dt,sig1,sig2,sig3,sig4)

segLength = 10 %segment length in seconds

pointsPerSeg = round(segLength/dt);

numSegs = floor(length(sig1)/pointsPerSeg);
% Fs = 2016;
%         lof = .1; hif = 2; nyqf = Fs/2; %set bandpass range
%         lof = lof/nyqf; hif = hif/nyqf;
%         [b,a] = butter(1, [lof hif]);
sig1 = sig1 - mean(sig1);
sig1 = sig1/std(sig1);
if nargin>3
sig2 = sig2 - mean(sig2);
sig2 = sig2/std(sig2);
end


sig1_cent = (quantile(sig1,0.9) + quantile(sig1,0.1))/2;

if nargin>3
    sig2_cent = (quantile(sig2,0.9) + quantile(sig2,0.1))/2;
end

for i = 1:numSegs

    if nargin == 3
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),sig1(pointsPerSeg*(i-1)+1:pointsPerSeg*i))
%         ylim([-1500 2000])
    elseif nargin == 4
        seg1 = sig1(pointsPerSeg*(i-1)+1:pointsPerSeg*i);
        seg2 = sig2(pointsPerSeg*(i-1)+1:pointsPerSeg*i);
%         seg1 = filtfilt(b,a,seg1);
%         seg2 = filtfilt(b,a,seg2);
%         seg1 = detrend(seg1);
%         seg2 = detrend(seg2);
%         seg1=  seg1 - mean(seg1);
%         seg1 = seg1/std(seg1);
%         seg2 = seg2 - mean(seg2);
%         seg2 = seg2/std(seg2);
%         max1 = max(seg1);
%         min1 = min(seg1);
%         max2 = max(seg2);
%         min2 = min(seg2);
%         thr1 = (max1+min1)/2;
%         thr2 = (max2+min2)/2;
        subplot(2,1,1)
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),seg1','linewidth',2);hold on
                        line([synct(pointsPerSeg*(i-1)+1) synct(pointsPerSeg*i)],[sig1_cent sig1_cent],'Color','b')
% line([synct(pointsPerSeg*(i-1)+1) synct(pointsPerSeg*i)],[0.5 0.5],'Color','b')
        subplot(2,1,2)
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),seg2,'r','linewidth',2)
        line([synct(pointsPerSeg*(i-1)+1) synct(pointsPerSeg*i)],[sig2_cent sig2_cent],'Color','r')
%         line([synct(pointsPerSeg*(i-1)+1) synct(pointsPerSeg*i)],[0.5 0.5],'Color','r')

%         line([synct(pointsPerSeg*(i-1)+1) synct(pointsPerSeg*i)],[0 0],'Color','k')
%                 line([synct(pointsPerSeg*(i-1)+1) synct(pointsPerSeg*i)],[-0.5 -0.5],'Color','k')

    elseif nargin == 5
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),sig1(pointsPerSeg*(i-1)+1:pointsPerSeg*i));hold on
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),sig2(pointsPerSeg*(i-1)+1:pointsPerSeg*i),'r')
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),sig3(pointsPerSeg*(i-1)+1:pointsPerSeg*i),'g')
        ylim([-1500 2000])
    elseif nargin == 6
        sig1R = max(sig1(pointsPerSeg*(i-1)+1:pointsPerSeg*i)) - min(sig1(pointsPerSeg*(i-1)+1:pointsPerSeg*i));hold on
        sig2R = max(sig2(pointsPerSeg*(i-1)+1:pointsPerSeg*i)) - min(sig2(pointsPerSeg*(i-1)+1:pointsPerSeg*i));
        sig3R = max(sig3(pointsPerSeg*(i-1)+1:pointsPerSeg*i)) - min(sig3(pointsPerSeg*(i-1)+1:pointsPerSeg*i));
        sig4R = max(sig1(pointsPerSeg*(i-1)+1:pointsPerSeg*i)) - min(sig4(pointsPerSeg*(i-1)+1:pointsPerSeg*i));
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),1/2*sig2R+sig1(pointsPerSeg*(i-1)+1:pointsPerSeg*i),'linewidth',2)
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),sig1R+1/2*sig2R+sig2(pointsPerSeg*(i-1)+1:pointsPerSeg*i),'r','linewidth',2)
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),-1/2*sig3R+sig3(pointsPerSeg*(i-1)+1:pointsPerSeg*i),'g','linewidth',2)
        plot(synct(pointsPerSeg*(i-1)+1:pointsPerSeg*i),-1/2*sig4R-sig3R+sig4(pointsPerSeg*(i-1)+1:pointsPerSeg*i),'r','linewidth',2)
        ylim([-sig3R-1.5*sig4R sig1R+1.5*sig2R])
    end

    
    pause
    clf

end



%plot large data set in segments
function [] = dataViewer_t(synct,dt,sig1,times1,times2,sig2,times1b,times2b)

segLength = 20 %segment length in seconds
lnwdth = 2;

firstTime = synct(1);
synct = synct-firstTime;
times1 = times1 - firstTime;
if nargin > 4
times2 = times2 - firstTime;
end
% times1ID = zeros(size(sig1));
% times2ID = times1ID;
% 
% times1ID(round(times1/dt)) = 1;
% times2ID(round(times2/dt)) = 1;

if nargin > 6
    times1b = times1b - firstTime;
    times2b = times2b - firstTime;
end
pointsPerSeg = round(segLength/dt);

numSegs = floor(length(sig1)/pointsPerSeg);

for i = 1:numSegs

    begID = pointsPerSeg*(i-1)+1;
    endID = pointsPerSeg*i;

    %zscore signals within segment window
%     sig1(begID:endID) = sig1(begID:endID) - mean(sig1(begID:endID));
%     sig1(begID:endID) = sig1(begID:endID)/std(sig1(begID:endID));
%     if nargin > 5
%         sig2(begID:endID) = sig2(begID:endID) - mean(sig2(begID:endID));
%         sig2(begID:endID) = sig2(begID:endID)/std(sig2(begID:endID));
%     end

    begTime = (i-1)*segLength;
    endTime = i*segLength;
    plot(synct(begID:endID),sig1(begID:endID),'linewidth',lnwdth)
    hold on
    plot(times1(times1>begTime&times1<endTime),sig1(round(times1(times1>begTime&times1<endTime)/dt)),'go','linewidth',4,'MarkerSize',4)
    if nargin > 4
    plot(times2(times2>begTime&times2<endTime),sig1(round(times2(times2>begTime&times2<endTime)/dt)),'ro','linewidth',4,'MarkerSize',4)
    end

    if nargin > 5
        plot(synct(begID:endID),sig2(begID:endID),'r','linewidth',lnwdth)
        if nargin > 6
            plot(times1b(times1b>begTime&times1b<endTime),sig2(round(times1b(times1b>begTime&times1b<endTime)/dt)),'g*','linewidth',4,'MarkerSize',4)
            plot(times2b(times2b>begTime&times2b<endTime),sig2(round(times2b(times2b>begTime&times2b<endTime)/dt)),'r*','linewidth',4,'MarkerSize',4)
        end
    end


    pause
    clf


end






%first downsample to 504Hz
dsf = 4;
Fsd = 2016/dsf;
wcv_d = downsample(wcv_minus_spike,dsf);

min_effect_size = 2;

%zscore
wcv_z = zscore(wcv_d);

%filter low and medium
niqf = Fsd/2;
hcf1 = 1/niqf;
hcf2 = 5/niqf;

[b1,a1] = butter(2,hcf1,'low');
[b2,a2] = butter(2,hcf2,'low');

wcv_f1 = filtfilt(b1,a1,wcv_z);
wcv_f2 = filtfilt(b2,a2,wcv_z);

segTime = 60;
segLength = segTime*Fsd; 
numSegs = floor(length(wcv_z)/segLength);
seg_Taxis = [1:numSegs]*segTime-round(segTime)/2;
up_cross = [];
down_cross = [];
for i = 1:numSegs
    
   segBeg = (i-1)*segLength+1;
   segEnd = segBeg+segLength;
   
   wseg_f1 = wcv_f1(segBeg:segEnd);
   wseg_f2 = wcv_f2(segBeg:segEnd);
   
   %get GMM for wcv in segment
   
    wseg_d2 = downsample(wseg_f2,2);
    [u,sig,t,iter] = fit_mix_gaussian(wseg_d2,2);
    effect_size(i) = get_effect_size(u,sig)
    if effect_size(i) > 1.5
        [threshold(i)] = find_GMM_xp(u,sig,t)
        if abs(threshold(i)) > 2
            threshold(i) = (max(wseg_f1)+min(wseg_f1))/2;
        end
    else
    threshold(i) = (max(wseg_f1) + min(wseg_f1))/2;
    end
    
    %find threshold crossings
    thresh_det = wseg_f1;
    thresh_det(thresh_det>threshold(i)) = 1;
    thresh_det(thresh_det<threshold(i)) = 0;
    
    dtd = [0;diff(thresh_det)];
    up_cross = [up_cross; (find(dtd==1) + segBeg)];
    down_cross = [down_cross; (find(dtd == -1) + segBeg)];
    
%     plot(wseg_f1)
%     hold on
%     plot(down_cross,wseg_f1(down_cross),'ro')
%     plot(up_cross,wseg_f1(up_cross),'go')
%     pause;clf
    
        
end

up_cross_times = up_cross/Fsd;
down_cross_times = down_cross/Fsd;
up_thresh = spline(seg_Taxis,threshold,up_cross_times);
down_thresh = spline(seg_Taxis,threshold,down_cross_times);
up_thresh(up_cross_times<seg_Taxis(1)) = threshold(1);
up_thresh(up_cross_times>seg_Taxis(end)) = threshold(end);
down_thresh(down_cross_times<seg_Taxis(1)) = threshold(1);
down_thresh(down_cross_times>seg_Taxis(end)) = threshold(1);

for i = 1:length(up_cross)

    if wcv_f2(up_cross(i)) > up_thresh(i)
        temp = find(wcv_f2(up_cross(i)-2*Fsd:up_cross(i))<up_thresh(i),1,'last');
        if isempty(temp)
            up_trans(i) = NaN;
        else
            up_trans(i) = up_cross(i)-2*Fsd+temp;
        end
    else
        temp = find(wcv_f2(up_cross(i):up_cross(i)+2*Fsd)>up_thresh(i),1,'first');
        if isempty(temp)
            up_trans(i) = NaN;
        else
            up_trans(i) = up_cross(i)+temp;
        end
    end
        
end


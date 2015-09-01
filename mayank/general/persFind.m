%% FIND PERSISTENT ACTIVITY POINTS
windowSlide = 5;
windowSize = 50;
minCrossDur = 10;
Fs = 2016;

for i = 1:43

    upperCI = lf8_DCI{i}+2*std(lf8_DCI{i});

    persVec = wcv_DCI{i} - upperCI;


    %if all persVec vals are above 0 then all are pers points
    if min(persVec) > 0
        leftCrossings = 1;
        rightCrossings = length(persVec);
    elseif max(persVec) < 0
        leftCrossings = [];
        rightCrossings = [];
    else

        persCheck = zeros(size(persVec));
        persCheck(persVec>0) = 1;

        diffPersCheck = [0 diff(persCheck)];

        leftCrossings = find(diffPersCheck == 1);
        rightCrossings = find(diffPersCheck == -1);

        if isempty(leftCrossings)
            leftCrossings = 1;
        end

        %make sure you start with an up crossing
        if min(rightCrossings) < min(leftCrossings)

            leftCrossings = [1 leftCrossings];

        end

        %make sure you end with a right crossing
        if max(leftCrossings) > max(rightCrossings)

            rightCrossings = [rightCrossings length(persVec)];

        end


        persDur = rightCrossings - leftCrossings;

        %get rid of noise crossings
        leftCrossings(persDur < minCrossDur) = [];
        rightCrossings(persDur < minCrossDur) = [];

        downDur = leftCrossings(2:end) - rightCrossings(1:end-1);
        badDowns = find(downDur < minCrossDur);
        leftCrossings(badDowns+1) = [];
        rightCrossings(badDowns) = [];

    end

    %get avg DCI during pers
    avgDCI(i) = mean(wcv_DCI{i});
    if length(leftCrossings) > 0
        persDCI = [];
        for t = 1:length(leftCrossings)
            persDCI = [persDCI wcv_DCI{i}(leftCrossings(t):rightCrossings(t))];
        end
        avg_pers_DCI(i) = mean(persDCI);
        clear persDCI
    
    end
    
    if ~isempty(leftCrossings)
        %cross times in sec
        leftCrossTimes = (leftCrossings-1)*windowSlide+round(windowSize/2);
        rightCrossTimes = (rightCrossings-1)*windowSlide+round(windowSize/2);

        %get points where wcv is in persistent activity
        persPoints{i} = [];
        for t = 1:length(leftCrossTimes)

            tempLeftCrossPoint = round(Fs*leftCrossTimes(t));
            tempRightCrossPoint = round(Fs*rightCrossTimes(t));

            persPoints{i} = [persPoints{i} tempLeftCrossPoint:tempRightCrossPoint];

        end

    else
        persPoints{i} = [];
    end

%     plot(persVec)
%     hold on
%     plot(leftCrossings,zeros(size(leftCrossings)),'ro')
%     plot(rightCrossings,zeros(size(rightCrossings)),'go')
%     pause
%     clf

end


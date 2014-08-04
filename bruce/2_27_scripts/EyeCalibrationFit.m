function [EyePara DataPts] = EyeCalibrationFit(eyecal, FitMode, DISPLAY)


if nargin<2
    FitMode='linear';
end
if nargin<3
    DISPLAY = 0;
end
FixXu = unique(eyecal(5,:));
FixYu = unique(eyecal(6,:));

FixXu = FixXu(~isnan(FixXu));
FixYu = FixYu(~isnan(FixYu));
FixXu = FixXu(abs(FixXu)<7);
FixYu = FixYu(abs(FixYu)<7);
EyeCluster = [];
EyeClusterStd = [];
FixCenter = [];
EyeLdatas = [];
EyeRdatas = [];

DataPts = [];
for i=1:length(FixXu)
    for j=1:length(FixYu)
        inds = find(eyecal(5,:)==FixXu(i) & eyecal(6,:)==FixYu(j));
        if isempty(inds)
            continue;
        end
        
        eyecalij = eyecal(1:4, inds);
        centerij = median(eyecalij');  % estimated center of the eye data
        
        Ndata = size(eyecalij,2);
        dists = abs(eyecalij-centerij'*ones(1, Ndata));
        
        goodpt = find(any(dists>.5)==0);
        eyecalijgood = eyecalij(:, goodpt);
        
        Ngoodpt = length(goodpt);
        boundLOW = prctile(eyecalijgood', 0);
        boundHIGH = prctile(eyecalijgood', 100);
        
        InRange = eyecalijgood>boundLOW'*ones(1,Ngoodpt) & eyecalijgood<boundHIGH'*ones(1,Ngoodpt);
        Eyept = find((InRange(1,:) & InRange(2,:) & InRange(3,:) & InRange(4,:))>0);
        
        DataPts = [DataPts; reshape(inds(goodpt(Eyept)), [], 1)];
        
        CenterEst = zeros(1,4);
        NoiseEst = zeros(1,4);
        for k=1:4
            data = eyecalijgood(k,:);
            
            pts = Eyept;
            CenterEst(k) = mean(data(pts));
            NoiseEst(k) = std(data(pts));
        end
        if ~isempty(Eyept)
            EyeLdatas = [EyeLdatas; [eyecalijgood([1 2], Eyept)' ones(size(Eyept))'*[FixXu(i) FixYu(j)] ] ];
            EyeRdatas = [EyeRdatas; [eyecalijgood([3 4], Eyept)' ones(size(Eyept))'*[FixXu(i) FixYu(j)] ]];
            FixCenter = [FixCenter; [FixXu(i) FixYu(j)]];
            
            EyeCluster = [EyeCluster; CenterEst];
            EyeClusterStd = [EyeClusterStd; NoiseEst];
        end
    end
end
fprintf('Using %d data pts from %d fixations\n',length(DataPts), size(eyecal,2))

% figure,
% subplot(2,2,1);imagesc(reshape(EyeClusterStd(:,1), 5,5))
% subplot(2,2,2);imagesc(reshape(EyeClusterStd(:,2), 5,5))
% subplot(2,2,3);imagesc(reshape(EyeClusterStd(:,3), 5,5))
% subplot(2,2,4);imagesc(reshape(EyeClusterStd(:,4), 5,5))

% figure,subplot(2,2,1);plot(eyecal(1,:),eyecal(2,:),'.');
% subplot(2,2,2);plot(eyecal(1,DataPts),eyecal(2,DataPts),'.');axis([-10 10 -10 10])
% subplot(2,2,3);plot(eyecal(1,:),eyecal(2,:),'.');axis([-1 1 -1 1])
% subplot(2,2,4);plot(eyecal(1,DataPts),eyecal(2,DataPts),'.');axis([-1 1 -1 1])

gains = zeros(4, 2);
offsets = zeros(4, 1);
highergain = zeros(4, 3);
% funObj = @(p, x) x*p(1:2)'+p(3);
% costObj = @(p, x, y) sum((x*p(1:2)+p(3) - y).^2);
switch FitMode
    case 'linear'
    costObj = @(p, x, y) sum(abs(x*p(1:2)+p(3) - y));
    p0 = [1 1 0]';
    case 'offsetonly'
    costObj = @(p, x, y) sum(abs(x+p(1) - y));
    p0 = 0;
    gains = [1 0; 0 1; 1 0; 0 1];
    case 'highergain'
        costObj = @(p, x, y) sum((x*p(1:2)+p(3) + p(4)*x(:,1).^2+ p(5)*x(:,2).^2+ p(6)*x(:,1).*x(:,2) - y).^2);
        p0 = [1 1 0 0 0 0]';
end


opts = optimset('GradObj','off','Display','off','MaxIter',10000,'MaxFunEvals',10000,'TolFun',1e-7);
for eyei=1:2
    for hv=1:2
        EyeMeasure = EyeCluster(:,[eyei eyei+2]);
        
        switch FitMode
            case 'linear'
                if eyei==1
                    beta = fminsearch(@(p)costObj(p, EyeLdatas(:,1:2),EyeLdatas(:,2+hv)),p0,opts);
                elseif eyei==2
                    beta = fminsearch(@(p)costObj(p, EyeRdatas(:,1:2),EyeRdatas(:,2+hv)),p0,opts);
                end
                
                gains(2*(eyei-1)+hv,:) = beta(1:2);
                offsets(2*(eyei-1)+hv,:) = beta(3);
            case 'offsetonly'
                if eyei==1
                    beta = fminsearch(@(p)costObj(p, EyeLdatas(:,hv),EyeLdatas(:,2+hv)),p0,opts);
                elseif eyei==2
                    beta = fminsearch(@(p)costObj(p, EyeRdatas(:,hv),EyeRdatas(:,2+hv)),p0,opts);
                end
                offsets(2*(eyei-1)+hv,:) = beta(1);
            case 'highergain'
                if eyei==1
                    beta = fminsearch(@(p)costObj(p, EyeLdatas(:,1:2),EyeLdatas(:,2+hv)),p0,opts);
                elseif eyei==2
                    beta = fminsearch(@(p)costObj(p, EyeRdatas(:,1:2),EyeRdatas(:,2+hv)),p0,opts);
                end
                gains(2*(eyei-1)+hv,:) = beta(1:2);
                offsets(2*(eyei-1)+hv,:) = beta(3);                
                highergain(2*(eyei-1)+hv,1:(length(beta)-3)) = beta(4:end);
        end
                
        %                 figure,plot(funObj(beta, EyeMeasure), FixCenter(:,hv),'r.')
    end
end


EyePara = struct('gain', gains, 'offset', offsets,'highergain',highergain);
EyePara.label = {'LeyeH','LeyeV','ReyeH','ReyeV'};

if ~DISPLAY return; end


% after calibration
EyeLdatascal = CalibrateEyeSignal(EyeLdatas(:,1:2), gains(1:2,:),offsets(1:2), highergain(1:2,:));
EyeRdatascal = CalibrateEyeSignal(EyeRdatas(:,1:2), gains(3:4,:),offsets(3:4), highergain(3:4,:));
% estimate error
rMSEs = zeros(4,length(FixXu), length(FixYu));
meanDIST = zeros(4,length(FixXu), length(FixYu));
Vergence = zeros(2,length(FixXu), length(FixYu));
for i=1:length(FixXu)
    for j=1:length(FixYu)
        % left
        Linds = find(EyeLdatas(:,3)==FixXu(i) & EyeLdatas(:,4)==FixYu(j));
        rMSEs(1,i,j) = sqrt(mean((EyeLdatascal(Linds,1) - FixXu(i)).^2));
        rMSEs(2,i,j) = sqrt(mean((EyeLdatascal(Linds,2) - FixYu(j)).^2));
        meanDIST(1,i,j) = (median(EyeLdatascal(Linds,1)) - FixXu(i));
        meanDIST(2,i,j) = (median(EyeLdatascal(Linds,2)) - FixYu(j));
        
        % right
        Rinds = find(EyeRdatas(:,3)==FixXu(i) & EyeRdatas(:,4)==FixYu(j));
        rMSEs(3,i,j) = sqrt(mean((EyeRdatascal(Rinds,1) - FixXu(i)).^2));
        rMSEs(4,i,j) = sqrt(mean((EyeRdatascal(Rinds,2) - FixYu(j)).^2));  
        meanDIST(3,i,j) = (median(EyeRdatascal(Rinds,1)) - FixXu(i));
        meanDIST(4,i,j) = (median(EyeRdatascal(Rinds,2)) - FixYu(j));
        % vergence
        Vergence(1:2,i,j) = mean(EyeLdatascal(Linds,1:2) - EyeRdatascal(Rinds,1:2));        
    end
end

% display clustering data
axissize = [FixXu(1) FixXu(end) FixYu(1) FixYu(end)]*1.2;
h=figure(1);clf
subplot(2,2,1),plot(EyeLdatas(:,1), EyeLdatas(:,2),'.');hold on;
plot(FixCenter(:,1), FixCenter(:,2),'r.','MarkerSize',12);axis square
title('Before Calibration-L');axis(axissize);

subplot(2,2,2),plot(EyeRdatas(:,1), EyeRdatas(:,2),'.');hold on;
plot(FixCenter(:,1), FixCenter(:,2),'r.','MarkerSize',12);axis square
title('Before Calibration-R');axis(axissize);
subplot(2,2,3),plot(EyeLdatascal(:,1), EyeLdatascal(:,2),'.');hold on;
plot(FixCenter(:,1), FixCenter(:,2),'r.','MarkerSize',12);axis square
title('After Calibration-L');axis(axissize);

subplot(2,2,4),plot(EyeRdatascal(:,1), EyeRdatascal(:,2),'.');hold on;
plot(FixCenter(:,1), FixCenter(:,2),'r.','MarkerSize',12);axis square
title('After Calibration-R');axis(axissize);


% display error
figure(2),clf; % root-mean square error
subplot(2,2,1);imagesc(squeeze(rMSEs(1,:,:)));colorbar;title('L Eye-H');
subplot(2,2,2);imagesc(squeeze(rMSEs(2,:,:)));colorbar;title('L Eye-V');
subplot(2,2,3);imagesc(squeeze(rMSEs(3,:,:)));colorbar;title('R Eye-H');
subplot(2,2,4);imagesc(squeeze(rMSEs(4,:,:)));colorbar;title('R Eye-V');
figure(3),clf; % average distance from the center
subplot(2,2,1);imagesc(squeeze(meanDIST(1,:,:)));colorbar;title('L Eye-H');
subplot(2,2,2);imagesc(squeeze(meanDIST(2,:,:)));colorbar;title('L Eye-V');
subplot(2,2,3);imagesc(squeeze(meanDIST(3,:,:)));colorbar;title('R Eye-H');
subplot(2,2,4);imagesc(squeeze(meanDIST(4,:,:)));colorbar;title('R Eye-V');
figure(4),clf; % vergence
subplot(2,2,1);imagesc(squeeze(Vergence(1,:,:)));colorbar;title('H: L-R')
subplot(2,2,2);imagesc(squeeze(Vergence(2,:,:)));colorbar;title('V: L-R')
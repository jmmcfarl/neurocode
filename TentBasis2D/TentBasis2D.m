classdef TentBasis2D
    % TentBasis2D is a 2D piecewise linear basis
    
    properties
        Xtick;
        Ytick;
        
        gridX; % tent basis center pt
        gridY; % tent basis center pt
        
        tentcoeff;  % structure array
        dNLx;       % derivative (Nx,Ny,2) the two entries are for lower and
                    % upper triangle of the bin(i,j)
        dNLy;       % derivatives along y direction
        
        Counts; % sample number per vertice, note that each sample will be used by 3
                % vertices, so sum(Counts(:)) is approximately 3 times of
                % the sample number                        
                
        prescale;   % scale matrix (2x2) also do whitening
    end
    
    methods
        function TB = TentBasis2D(Xtick, Ytick)
            Nx = length(Xtick);
            Ny = length(Ytick);
            
            TB.Xtick = Xtick;
            TB.Ytick = Ytick;
            TB.gridX = zeros(Nx,Ny);
            TB.gridY = zeros(Nx,Ny);
            TB.tentcoeff = zeros(Nx,Ny);
            TB.prescale = eye(2,2);
            TB.Counts = zeros(Nx,Ny); % sample per vertice

            for x=1:Nx
                for y=1:Ny
                    TB.gridX(x,y) = Xtick(x);
                    TB.gridY(x,y) = Ytick(y);                    
                end
            end
        end
        
        function TB = NLderivative(TB)
            Nx = length(TB.Xtick);
            Ny = length(TB.Ytick);
            TB.dNLx = zeros(Nx-1,Ny-1,2);
            TB.dNLy = zeros(Nx-1,Ny-1,2);
            NLx = TB.Xtick;
            NLy = TB.Ytick;
            for x=1:Nx-2
                for y=1:Ny-2
                    % vertices of lower triangle
                    LT = [NLx(x) NLy(y); NLx(x+1) NLy(y); NLx(x) NLy(y+1)];
                    LTf = [TB.tentcoeff(x,y) TB.tentcoeff(x+1,y) TB.tentcoeff(x,y+1)];
                    % vertices of upper triangle
                    UT = [NLx(x+1) NLy(y); NLx(x+1) NLy(y+1); NLx(x) NLy(y+1)];
                    UTf = [TB.tentcoeff(x+1,y) TB.tentcoeff(x+1,y+1) TB.tentcoeff(x,y+1)];
                    [TB.dNLx(x,y,1) TB.dNLy(x,y,1)]=TB.PieceSlope(LT,LTf);
                    [TB.dNLx(x,y,2) TB.dNLy(x,y,2)]=TB.PieceSlope(UT,UTf);
                end
            end
        end
        
        function output = Process(TB,stim)
                        
            
            NLmat = InputNL2D(TB,stim);
            output = NLmat*TB.tentcoeff(:);
%             keyboard
            return;
            % the following code gives same result, just for sanity check
%             stim = TB.CutOffset(stim);
%             NT = size(stim,1);
%             NNx = length(TB.Xtick);
%             NNy = length(TB.Ytick);
%             
%             output = zeros(NT,1);
%             NLx = TB.Xtick;
%             NLy = TB.Ytick;
%             
%             for t=1:NT
%                 x = find(TB.Xtick<stim(t,1),1,'last');
%                 y = find(TB.Ytick<stim(t,2),1,'last');
%                 if isempty(x) || isempty(y)
%                     continue
%                 end
%                 if x>NNx-2 || y>NNy-2
%                     continue;
%                 end
%                 
%                 % vertices of lower triangle
%                 LT = [NLx(x) NLy(y); NLx(x+1) NLy(y); NLx(x) NLy(y+1)];
%                 % vertices of upper triangle
%                 UT = [NLx(x+1) NLy(y); NLx(x+1) NLy(y+1); NLx(x) NLy(y+1)];
%                 flag = PointInTriangle(stim(t,:), LT(1,:),LT(2,:),LT(3,:)) ;
%                 
%                 if flag==1
%                     lambdasLT = BaryCentric(stim(t,:), LT(1,:),LT(2,:),LT(3,:)) ;
%                     output(t) = output(t)+ lambdasLT(:,1)*TB.tentcoeff(x,y);
%                     output(t) = output(t)+ lambdasLT(:,2)*TB.tentcoeff(x+1,y);
%                     output(t) = output(t)+ lambdasLT(:,3)*TB.tentcoeff(x,y+1);
%                 else
%                     lambdasUT = BaryCentric(stim(t,:), UT(1,:),UT(2,:),UT(3,:)) ;
%                     output(t) = output(t)+ lambdasUT(:,1)*TB.tentcoeff(x+1,y);
%                     output(t) = output(t)+lambdasUT(:,2)*TB.tentcoeff(x+1,y+1);
%                     output(t) = output(t)+ lambdasUT(:,3)*TB.tentcoeff(x,y+1);
%                 end
%             end            
        end
        
        function [stim outstimi] = CutOffset(TB,stim)
            % cut off stimulus at boundaries            
            
            NLx = TB.Xtick;
            NLy = TB.Ytick;
            if nargout>=2
                outstimi = find(stim(:,1)<=NLx(1) | stim(:,2)<=NLy(1) | ...
                    stim(:,1)>=NLx(end) | stim(:,2)>=NLy(end));
            end
            
            stim(stim(:,1)<=NLx(1),1) = NLx(1);
            stim(stim(:,2)<=NLy(1),2) = NLy(1);
            stim(stim(:,1)>=NLx(end),1) = NLx(end)-1e-7;
            stim(stim(:,2)>=NLy(end),2) = NLy(end)-1e-7;
            
        end
        
        function output = Process_derivative(TB,stim)
            % process stimulus using derivative
            % stim must be a 2D signal (NT,2)
            % output will also be 2D, corresponding to derivative along X
            % and Y direction
            [NT Ndim] = size(stim);
            % rescale stimulus to match NL operating domain
            stim = stim*TB.prescale;
            stim = TB.CutOffset(stim);
            if Ndim~=2
                disp('Input must be a 2D signal');
                return
            end
            
            NLx = TB.Xtick;
            NLy = TB.Ytick;
            % number of sample pts along x/y dimensions
            NNx = length(NLx);
            NNy = length(NLy);            
            
            % bin centers
            [n, binx] =histc(stim(:,1),NLx);
            [n, biny] =histc(stim(:,2),NLy);
            
            output = zeros(NT,2);
            
            % go through all bins
            for x=1:NNx-1
                for y=1:NNy-1                    
                        
                    % data points in the current bin
                    inds=find(binx(:)==x & biny(:)==y);
                    if isempty(inds)
                        continue;
                    end
                    
                    dfdxLT = TB.dNLx(x,y,1);
                    dfdyLT = TB.dNLy(x,y,1);
                    dfdxUT = TB.dNLx(x,y,2);
                    dfdyUT = TB.dNLy(x,y,2);
                    % vertices of lower triangle
                    LT = [NLx(x) NLy(y); NLx(x+1) NLy(y); NLx(x) NLy(y+1)];
                    % vertices of upper triangle
                    UT = [NLx(x+1) NLy(y); NLx(x+1) NLy(y+1); NLx(x) NLy(y+1)];
                    
                    % decide whether data pt lies in the upper/lower
                    % triangle of the current square
                    flag = PointInTriangle(stim(inds,:), LT(1,:),LT(2,:),LT(3,:)) ;
                    indLT = inds(flag==1);
                    indUT = inds(flag==0);
                    
                    if ~isempty(indLT)
                        output(indLT,1) = dfdxLT;
                        output(indLT,2) =  dfdyLT;
                    end
                    
                    if ~isempty(indUT)
                        output(indUT,1) = dfdxUT;
                        output(indUT,2) =  dfdyUT;
                    end
                end
            end
            output = output*(TB.prescale');
%             output(:,1) = output(:,1)*TB.prescale(1);
%             output(:,2) = output(:,2)*TB.prescale(2);
        end
        
        function DisplayDetail(TB, Mask)
            baseline = 0;
            if nargin<2
                Mask = ones(size(TB.tentcoeff));
            end
            [Nx Ny] = size(TB.tentcoeff);
%             centerX = find(abs(TB.Xtick)==min(abs(TB.Xtick)),1);
%             centerY = find(abs(TB.Ytick)==min(abs(TB.Ytick)),1);
%             baseline = TB.tentcoeff(centerX,centerY);
 
            subplot(2,2,1);
            imagesc(TB.Xtick,TB.Xtick,(TB.tentcoeff-baseline).*Mask);
            set(gca,'Ydir','normal');axis square
            colormap jet
            
            subplot(2,2,2);            
            cmap = jet(Nx);
            for x=1:Nx
                NL = TB.tentcoeff(x,:);
                inds = find(NL~=baseline);
                plot(TB.Ytick(inds), NL(inds),'Color',cmap(x,:));hold on;                
            end
%             colormap cool
            colorbar;caxis([TB.Xtick(1) TB.Xtick(end)]);
            xlabel('supression')
            
            subplot(2,2,4);            
            cmap = jet(Ny);
            for y=1:Ny
                NL = TB.tentcoeff(:,y);
                inds = find(NL~=baseline);
                plot(TB.Xtick(inds), NL(inds),'Color',cmap(y,:));hold on;
            end
%             colormap cool
            colorbar;caxis([TB.Xtick(1) TB.Xtick(end)]);
            xlabel('excitation')            
        end
        
        function Display(TB)
            baseline = 0;
            if nargin<2
                Mask = ones(size(TB.tentcoeff));
            end
            [Nx Ny] = size(TB.tentcoeff);
%             centerX = find(abs(TB.Xtick)==min(abs(TB.Xtick)),1);
%             centerY = find(abs(TB.Ytick)==min(abs(TB.Ytick)),1);
%             baseline = TB.tentcoeff(centerX,centerY);
            
            TB.tentcoeff(Mask==0) = nan;
            cmap = colormap;
            cmap(1,:) = [1 1 1];
            imagesc(TB.Xtick,TB.Xtick,(TB.tentcoeff-baseline));
            set(gca,'Ydir','normal');axis square
        end
        
        function TB = CalculateScaleMatrix(TB,gsig,Xrangefrac,decorr)
            if nargin<4
                decorr = 0;
            end
            
            if decorr==1
                % covariance of the signal
                sigma = (gsig'*gsig)/length(gsig);
                Wtrans = chol(sigma); % R'*R gives sigma
                Wtrans = inv(Wtrans);
                % Note: gsig*Wtrans gives I, the signal is decorrelated                
            else
                Wtrans = [1 0; 0 1];                
            end
            gsigW = gsig*Wtrans;
            
            eqscale = 1;            
            if Xrangefrac<100
                if eqscale                    
                    scale1 = TB.Xtick(end)/prctile(gsigW(:),Xrangefrac);
                    scale2 = scale1;
                else
                    scale1 = TB.Xtick(end)/prctile(gsigW(:,1),Xrangefrac);
                    scale2 = TB.Ytick(end)/prctile(gsigW(:,2),Xrangefrac);
                end
            else
                if eqscale                    
                    scale1 = TB.Xtick(end)/prctile(gsigW(:),100)*Xrangefrac/100;                
                    scale2 = scale1;
                else
                    scale1 = TB.Xtick(end)/prctile(gsigW(:,1),100)*Xrangefrac/100;
                    scale2 = TB.Ytick(end)/prctile(gsigW(:,2),100)*Xrangefrac/100;
                end
            end                        
            
            TB.prescale = Wtrans*[scale1 0; 0 scale2];
        end
        
        function [TB offset] = setK(TB,K,offset)
            % set kernel values
            if nargin<3
                offset=0;
            end
            
            TB.tentcoeff = reshape(K(offset+(1:length(TB.tentcoeff(:)))), ...
                length(TB.Xtick),length(TB.Ytick));
            offset = offset+length(TB.tentcoeff(:));            
        end
        
        function Counts = SignalDist(TB, stim)
            stim = stim*TB.prescale;
            stim = TB.CutOffset(stim);
            
            [NT Ndim] = size(stim);
            if Ndim~=2
                disp('Input must be a 2D signal');
                return
            end
            NLx = TB.Xtick;
            NLy = TB.Ytick;
            % number of sample pts along x/y dimensions
            NNx = length(NLx);
            NNy = length(NLy);
            
            % bin centers
            [n, binx] =histc(stim(:,1),NLx);
            [n, biny] =histc(stim(:,2),NLy);
                        
            Counts = zeros(NNx,NNy);
            
            % go through all bins
            for x=1:NNx-1
                for y=1:NNy-1
                    
                    % data points in the current bin
                    inds=find(binx(:)==x & biny(:)==y);
                    if isempty(inds)
                        continue;
                    end
                    
                    % vertices of lower triangle
                    LT = [NLx(x) NLy(y); NLx(x+1) NLy(y); NLx(x) NLy(y+1)];
                    % vertices of upper triangle
                    UT = [NLx(x+1) NLy(y); NLx(x+1) NLy(y+1); NLx(x) NLy(y+1)];
                    flag = PointInTriangle(stim(inds,:), LT(1,:),LT(2,:),LT(3,:)) ;
                    indLT = inds(flag==1);
                    indUT = inds(flag==0);
                                        
                    NLT = length(indLT);
                    Counts(x,y) = Counts(x,y)+NLT;
                    Counts(x+1,y) = Counts(x+1,y)+NLT;
                    Counts(x,y+1) = Counts(x,y+1)+NLT;
                    
                    NUT = length(indUT);
                    Counts(x+1,y) = Counts(x+1,y)+NUT;
                    Counts(x+1,y+1) = Counts(x+1,y+1)+NUT;
                    Counts(x,y+1) = Counts(x,y+1)+NUT;
                 end
            end                        
        end
        
        function [NLoutput Counts binx biny] = InputNL2D(TB,stim)
            % calculate output
            % Input:
            %       stim(NT,2): 2D signal
            %       NLx/NLy, bin boundaries on X/Y direction, if NLx is of length N, there
            %       will be (N-1) bins, the ith bin given by [NLx(i) NLx(i+1)]
            % Output:
            %       NLoutput(NT, Npar), where Npar=Nx*Ny
            %       Counts(NT, NS/2): number of samples per bin
            % Yuwei Cui, Created by Oct 20, 2012
            % Oct 28, 2012 Extend Pixel filter to 2D tent filter
            
%             stim(:,1) = stim(:,1)*TB.prescale(1);
%             stim(:,2) = stim(:,2)*TB.prescale(2);
            stim = stim*TB.prescale;
            [NT Ndim] = size(stim);
            stim = TB.CutOffset(stim);
            
            if Ndim~=2
                disp('Input must be a 2D signal');
                return
            end
            NLx = TB.Xtick;
            NLy = TB.Ytick;
            % number of sample pts along x/y dimensions
            NNx = length(NLx);
            NNy = length(NLy);
            
            % bin centers
            [n, binx] =histc(stim(:,1),NLx);
            [n, biny] =histc(stim(:,2),NLy);
            
            NLoutput=zeros(NT, NNx,NNy);
            if nargout>1
                Counts = zeros(NNx,NNy);
            end
            %                         keyboard
            % go through all bins
            for x=1:NNx-1
                for y=1:NNy-1                
                    % data points in the current bin
                    inds=find(binx(:)==x & biny(:)==y);
                    if isempty(inds)
                        continue;
                    end
                    
                    % vertices of lower triangle
                    LT = [NLx(x) NLy(y); NLx(x+1) NLy(y); NLx(x) NLy(y+1)];
                    % vertices of upper triangle
                    UT = [NLx(x+1) NLy(y); NLx(x+1) NLy(y+1); NLx(x) NLy(y+1)];
                    flag = PointInTriangle(stim(inds,:), LT(1,:),LT(2,:),LT(3,:)) ;
                    indLT = inds(flag==1);
                    indUT = inds(flag==0);
                    
                    lambdasLT = BaryCentric(stim(indLT,:), LT(1,:),LT(2,:),LT(3,:)) ;
                    
                    NLoutput(indLT,x,y) = NLoutput(indLT,x,y)+ lambdasLT(:,1);
                    NLoutput(indLT,x+1,y) = NLoutput(indLT,x+1,y)+ lambdasLT(:,2);
                    NLoutput(indLT,x,y+1) = NLoutput(indLT,x,y+1)+ lambdasLT(:,3);
                    
                    lambdasUT = BaryCentric(stim(indUT,:), UT(1,:),UT(2,:),UT(3,:)) ;
                    
                    NLoutput(indUT,x+1,y) = NLoutput(indUT,x+1,y)+ lambdasUT(:,1);
                    NLoutput(indUT,x+1,y+1) = NLoutput(indUT,x+1,y+1)+ lambdasUT(:,2);
                    NLoutput(indUT,x,y+1) = NLoutput(indUT,x,y+1)+ lambdasUT(:,3);
                    
                    if nargout>1
                    % count sample per vertice
                    NLT = length(indLT);
                    Counts(x,y) = Counts(x,y)+NLT;
                    Counts(x+1,y) = Counts(x+1,y)+NLT;
                    Counts(x,y+1) = Counts(x,y+1)+NLT;
                    
                    NUT = length(indUT);
                    Counts(x+1,y) = Counts(x+1,y)+NUT;
                    Counts(x+1,y+1) = Counts(x+1,y+1)+NUT;
                    Counts(x,y+1) = Counts(x,y+1)+NUT;
                    end
                end
            end
            NLoutput = reshape(NLoutput, NT, (NNx)*(NNy));
        end
        
        function [A B] = MakeConstr(TB, Xmono, Ymono, Mask)
            % Make constrained matrix, s.t.
            % A*TB.tentcoeff(:) < B
            % Xmono=1 if monotonic along X direction
            % Ymono=1 if monotonic along Y direction
            % Mask: apply the constraints only for pts within Mask (Mask=1)
            
            NLx = TB.Xtick;
            NLy = TB.Ytick;
            
            % number of sample pts along x/y dimensions
            NNx = length(NLx);
            NNy = length(NLy);
            ND = 0;
            if nargin<4
                Mask = ones(NNx,NNy);
            end
            
            if Xmono~=0
                ND = ND+(NNx-1)*NNy;
            end
            if Ymono~=0
                ND = ND+NNx*(NNy-1);
            end
            
            Npar = NNx*NNy;
            A = zeros(ND, Npar);
            B = zeros(ND,1);
            
            i=1; % index for constraint
            for x=1:NNx
                for y=1:NNy
                    % enforcing mono constraint along X axis
                    if x<NNx && Xmono~=0
                        if Mask(x,y)==0 || Mask(x+1,y) ==0
                            continue;
                        end
                        p1 = XY2pixN(x,y,NNx,NNy);
                        p2 = XY2pixN(x+1,y,NNx,NNy);
                        if Xmono<0 % decreasing along X
                            A(i,p1) = -1; A(i,p2) = 1;
                        elseif Xmono>0 % increasing along X
                            A(i,p1) = 1; A(i,p2) = -1;
                        end
                        i=i+1;
                    end
                    
                    % enforcing constraint along Y axis
                    if y<NNy && Ymono~=0
                        if Mask(x,y)==0 || Mask(x,y+1) ==0
                            continue;
                        end                        
                        p1 = XY2pixN(x,y,NNx,NNy);
                        p2 = XY2pixN(x,y+1,NNx,NNy);
                        if Ymono<0 % decreasing along Y
                            A(i,p1) = -1; A(i,p2) = 1;
                        elseif Ymono>0 % increasing along X
                            A(i,p1) = 1; A(i,p2) = -1;
                        end
                        i=i+1;
                    end                    
                end
            end
            A = A(1:(i-1),:);
            B = B(1:(i-1));
%             assert((i-1)==ND);       
            
%             testMat = TB.tentcoeff*0;            
%             x = 10; y=12;
%             p1 = XY2pixN(x,y,NNx,NNy);
%             p2 = XY2pixN(x+1,y,NNx,NNy);
%             testMat(p1) = 1; testMat(p2)=-1;
%             figure, imagesc(testMat);
        end
        
        function [NL1Dx NL1Dy sigdistx sigdisty Offset] = DecomposeNL2D(TB, weight,NoDisp)
            % Decompose the 2D nonlinearity f(x,y) into two 1D NL's
            % s.t f(x,y)~ fX(x) * fY(y)
            % Weighted MSE is minimized by the 1D NLs
            % disp('Decompose 2D into two 1D NLs');
            if nargin<3
                NoDisp=0;
            end
            opts = optimset('GradObj','on','Algorithm','active-set','Display','off','MaxIter',1000,'MaxFunEvals',1000,'TolFun',1e-6);
            Nx = length(TB.Xtick);
            Ny = length(TB.Ytick);
            NLexcInit = zeros(Nx,1);
            for i=1:Nx
                NLexcInit(i) = sum(TB.tentcoeff(i,:))/sum(TB.tentcoeff(i,:)~=0);
            end
            NLexcInit = NLexcInit-min(NLexcInit);
            
            NL1D = ones(Nx+Ny+1,1);
            NL1D(1:Nx) = (1:Nx)-1;
%             Mask = (TB.Counts>5);
%             Pmono = [1 Nx; -(Nx+1) -(Nx+Ny)];
            Pmono = [1 Nx]; % monotonic constraint on 1st filter
            Pcon = 1:(Nx+Ny); % postive constraint on everything
            
            [A B] = MakeConstrainMat(Nx+Ny+1, Pcon, Pmono);
            if nargin<2
                weight = ones(Nx,Ny);
            else
                weight = weight/sum(weight(:));
            end            
            Mask = weight>0;
            MSEthresh=0.001;
            MSEdiff = 100;
            MSEtrack = [];
            
            % Shift the Function to be all positive
            baseline = min(TB.tentcoeff(:));
            TB.tentcoeff = TB.tentcoeff-baseline;
            TB.tentcoeff(Mask==0) = 0;
            
            sigdistx = sum(TB.Counts,2); sigdistx=sigdistx/max(sigdistx);
            sigdisty = sum(TB.Counts,1); sigdisty=sigdisty/max(sigdisty);
            converge=0;
            iter=0;
            MaxIter = 20;
            while converge~=1
                % hold NLx, optimize NLy and offset
                hold_const =[];                
                hold_const(1,:) = (1:Nx)';
                hold_const(2,:) = NL1D(hold_const(1,:));
                [NL1Dopt MSEnow] = fmincon(@(k)MSE2DNL(k,TB.tentcoeff, weight, hold_const), NL1D, A, B, [],[],[],[],[], opts );
                NL1D = NL1Dopt;
                MSEtrack(end+1) = MSEnow;
                NLapprox = NL1D(1:Nx)*NL1D(Nx+(1:Ny))'+NL1D(end);
                
                % hold NLy, optimize NLx and offset
                hold_const =[];
                hold_const(1,:) = Nx+(1:Ny)';
                hold_const(2,:) = NL1D(hold_const(1,:));
                [NL1Dopt MSEnow]= fmincon(@(k)MSE2DNL(k,TB.tentcoeff, weight, hold_const), NL1D, A, B, [],[],[],[],[], opts );
                NL1D = NL1Dopt;
                MSEtrack(end+1) = MSEnow;
                NLapprox = NL1D(1:Nx)*NL1D(Nx+(1:Ny))'+NL1D(end);
                
                NL1Dx = NL1D(1:Nx);
                NL1Dy = NL1D(Nx+(1:Ny));
                Offset = NL1D(end);
                
                scale = max(NL1Dy(:)); 
                NL1Dy = NL1Dy/scale; NL1Dx = NL1Dx*scale;
                if ~NoDisp
                figure(7);clf;
                subplot(2,2,1);imagesc(NLapprox.*Mask);set(gca,'Ydir','normal');axis square;title('Approx');
                ca = caxis();
                subplot(2,2,2);imagesc(TB.tentcoeff);set(gca,'Ydir','normal');axis square;title('Original')               
                caxis(ca);
                subplot(2,2,3); plot(TB.Xtick,NL1Dx);hold on;plot(TB.Xtick,sigdistx*max(NL1Dx),'r--'); 
                title('f_X(x)');axis tight;axis square
                subplot(2,2,4); plot(TB.Ytick,NL1Dy);hold on;plot(TB.Ytick,sigdisty*max(NL1Dy),'r--'); 
                title('f_Y(y)');axis tight; axis square
                end
                fprintf('Last MSE %2.5f MSE now %2.5f\n',MSEtrack(end-1),MSEtrack(end));
                iter = iter+1;
                if length(MSEtrack)>2
                    MSEdiff = MSEtrack(end-2)-MSEtrack(end);
                    converge = MSEdiff/MSEtrack(end)<MSEthresh || iter>MaxIter;
                end
            end
            
            Offset = Offset+baseline;
        end
        
        
        function [TBr RotM Scales] = RotateNL(TB, alpha, beta)
            % rotate the 2D nonlinearity
            % Rotate x-axis by alpha
            % Rotate Y-axis by beta
            % to use the new nonlinearity, the signals should
            % (1) be rotated by RotM, e.g. gs = RotM * gs0(:,1:2)'
            % (2) be re-scaled by Scales 
            % gs(1,:)=gs(1,:)*Scales(1); gs(2,:)=gs(2,:)*Scales(1)
            % (3) TBr.Process(gs'); gives similar result as 
            % TB.Process(gs);
            
            RotM = [cos(alpha) sin(alpha);
                    -sin(beta) cos(beta)];
            RotMi = inv(RotM);
            
            Nx = length(TB.Xtick);
            Ny = length(TB.Ytick);
            
            % find value of
            TBr = TB;
            TBc = RotM * [TB.gridX(:) TB.gridY(:)]';
            Mask = TB.Counts>10;
            goodpt = find(Mask(:)==1);
%             figure, 
%             circle(0,0,1); hold on;
%             plot(TBc(1,:), TBc(2,:),'b.'); axis equal; 
%             plot(TBc(1,goodpt), TBc(2,goodpt),'r.'); axis equal; 
 
            MaxX= max(abs(TBc(1,goodpt)));
            MaxY= max(abs(TBc(2,goodpt)));
%             MaxX = MaxX*1.5;
%             MaxY = MaxY*1.5;
%             MaxX = 1;
%             MaxY = 1;
            % change scale
%             TBr.prescale(1,1) = TB.prescale(1,1)/MaxX;
%             TBr.prescale(2,2) = TB.prescale(2,2)/MaxY;
            
            % re-scale
            TBcmapX = TBr.gridX*MaxX;
            TBcmapY = TBr.gridY*MaxY;            
            % rotate back
            TBbackproj0 = RotMi * [TBcmapX(:) TBcmapY(:)]';
            % process with original NL
            TBbackproj(1,:) = TBbackproj0(1,:)/TB.prescale(1,1);
            TBbackproj(2,:) = TBbackproj0(2,:)/TB.prescale(2,2);
            TBr.tentcoeff = TB.Process(TBbackproj');
            
            [dum outstimi] = CutOffset(TB,TBbackproj0');
            TBr.tentcoeff(outstimi) = 0;
            
            TBr.tentcoeff = reshape(TBr.tentcoeff,Nx,Ny);
            
            
            
            Scales = 1./[MaxX MaxY];
%             
%             figure(5);
%             subplot(2,2,1); TB.Display; colorbar
%             subplot(2,2,2); TBr.Display; colorbar
%             
%             figure(10);clf
%             subplot(2,2,1);plot(TBr.gridX(:),TBr.gridY(:),'.'); axis equal
%             subplot(2,2,3);plot(TBcmapX(:),TBcmapY(:),'.'); axis equal
%             subplot(2,2,4); 
%             plot(TBbackproj(1,:), TBbackproj(2,:),'.'); axis equal; hold on;
%             plot(TBbackproj(1,outstimi), TBbackproj(2,outstimi),'g.'); axis equal; hold on;
%             plot(TB.gridX(goodpt)/TB.prescale(1,1),TB.gridY(goodpt)/TB.prescale(2,2),'r.');axis equal
        end
        
        function bin_areas = GetBinAreas(TB)
           bin_areas = ones(length(TB.Xtick),length(TB.Ytick));
           bin_areas(:,1) = 0.5;
           bin_areas(:,end) = 0.5;
           bin_areas(1,:) = 0.5;
           bin_areas(end,:) = 0.5;
           bin_areas(1,1) = 0.25;
           bin_areas(1,end) = 0.25;
           bin_areas(end,1) = 0.25;
           bin_areas(end,end) = 0.25;
        end
        
    end
    methods (Static)
        
        function Output=SingleTentBasis(stim, TBcenter, TBbase1, TBbase2)
            % stim is the 2D input signal
            % each 2D tent is defined by three 2D points
            % TBcenter is the tent center (top)
            % TBbase1 and TBbase2 are the two base point
            % Assume TBcenter is of height 1 and the two base of height 0
            if length(TBcenter)==2
                TBcenter(end+1) = 1;
            end
            if length(TBbase1)==2
                TBbase1(end+1) = 0;
            end
            if length(TBbase2)==2
                TBbase2(end+1) = 0;
            end
            NT = size(stim,1);
            % normal vector of the tent plane
            n=cross(TBcenter-TBbase1,TBcenter-TBbase2);
            Output = TBcenter(3) - ((stim(:,1)-TBcenter(1))*n(1)+...
                (stim(:,2)-TBcenter(2))*n(2))/n(3);
            
            flag = PointInTriangle(stim, TBcenter(1:2),TBbase1(1:2),TBbase2(1:2));
            Output(flag==0) =0;
        end
        
        function [dX dY] = PieceSlope(T,f)
            x = T(:,1);
            y = T(:,2);
            
            M = ((y(2)-y(3))*(x(1)-x(3))+(x(3)-x(2))*(y(1)-y(3)));
            
            dX = (f(1)*(y(2)-y(3)) + f(2)*(y(3)-y(1)) + ...
                f(3)*(-(y(2)-y(3))-(y(3)-y(1))) )/M;
            
            dY = (f(1)*(x(3)-x(2)) + f(2)*(x(1)-x(3))+...
                f(3)*(-(x(3)-x(2))-(x(1)-x(3))) )/M;
        end                       
    end
    
end
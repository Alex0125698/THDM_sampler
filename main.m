% 2HDM sampling density tests
% by A.S. Woodock
% JUN/2022
% License: GPL2
% 
% file contents: sampling density tests for the softly-broken Z2-aligned 2HDMs
%                tree-level only!

% setup
tic
clearvars
close all
clc
maxNumCompThreads(32)

% fill in parameters
p.nPoints = 10000000; % number of samples
p.v       = 246;      % Higgs VEV (GeV)
p.mode    = "Z2_symmetric"; % 

% run the sampler & spectrum generator to generate rest
% p = Sampler.generic(p);
p = Sampler.physical(p);
% p = Sampler.Higgs(p);

% calculate the likelihood
LL = ones(1,p.nPoints) == 1;

% LL = LL & TConstraint.vacuumStability(p);
% LL = LL & TConstraint.vacuumMetaStability(p);
% LL = LL & TConstraint.LOunitarity(p);
% LL = LL & TConstraint.massPositivity(p);

figure
histogram(sin(2*p.beta), 80, 'Normalization','probability');
% set(gca,'YScale','log')
xlabel('xx');
ylabel('sampling density');

figure
histogram(2*p.m122./sin(2*p.beta), 80, 'Normalization','probability');
% set(gca,'YScale','log')
xlabel('xx');
ylabel('sampling density');

%  make the plots
% allPlots(p, LL, "no constraints"); % full lambda range

outRange("lam1",p.lam1,0,4*pi);
outRange("lam2",p.lam2,0,4*pi);
outRange("lam3",p.lam3,-4*pi,4*pi);
outRange("lam4",p.lam4,-4*pi,4*pi);
outRange("lam5",p.lam5,-4*pi,4*pi);

toc

function outRange(name, lam, low, high)
    F_low = 100*sum(lam<low)/length(lam);
    F_high = 100*sum(lam>high)/length(lam);
    F = F_low+F_high;
    fprintf(1,"%s tot: %f low: %f high: %f\n", name, F, F_low, F_high);
end


function allPlots(p,avail,name)

    validCount = length(avail(avail~=0));
    totalCount = length(avail);
    per = (100*validCount) / totalCount;

    fprintf(1,"(%s) valid/total = %d / %d  (%f %%)\n", name, validCount, totalCount, per);
    
    if validCount > 0,

%         doSamplingDensityPlots(p,avail);
%         sgtitle(name);

        lambdaSamplingDensityPlots(p,avail)
        sgtitle(name);

        % massSamplingDensityPlots(p,avail);
        % sgtitle(name);

        % massDiffPlots(p,avail);
        % sgtitle(name);

    end
end

% --- main plotting functions ---

function doSamplingDensityPlots(p, avail)

    figure
    
    cosba = p.cosba(avail==1);
    tanb = p.tanb(avail==1);
    m122 = p.m122(avail==1);

    % don't plot m122 above 800 GeV
    m122(m122 > 1e6 | m122 < -1e6) = nan;
    
    % dont plot too many points otherwise matlab will lagg
    nKeep = min(length(cosba),1000000);
    cosba = cosba(1:nKeep); tanb = tanb(1:nKeep); m122 = m122(1:nKeep);

    % put m122 in mass dimensions
    m122 = sign(m122).*sqrt(abs(m122));

    % sampling density for cos(alpha-beta)
    subplot(2,3,1);
    histogram(cosba, 80, 'Normalization','probability');
    set(gca,'YScale','log')
    xlabel('cos(\beta-\alpha)');
    ylabel('sampling density');
    
    % sampling density for tan(beta)
    subplot(2,3,2);
    histogram(tanb, 80, 'Normalization','probability');
    set(gca,'YScale','log')
    xlabel('tan(\beta)');
    ylabel('sampling density');

    % sampling density for m122
    subplot(2,3,3);
    edges = linspace(-1e3,+1e3,20);
    histogram(m122,edges,'Normalization','probability')
    set(gca,'YScale','log')
    xlabel('m_1_2^2');
    ylabel('sampling density');
    
    % hist for tanb vs cosba plane 
    subplot(2,3,4);
    [N,c] = hist3([tanb;cosba]','Nbins',[40,40]);
    [X,Y] = meshgrid(c{1},c{2});
    pcolor(X,Y,N');
    set(gca,'ColorScale','log')
    colorbar
    xlabel('tan(\beta)');
    ylabel('cos(\alpha-\beta)');
    
    % hist for tanb vs cosba plane 
    subplot(2,3,5);
    [N,c] = hist3([m122;cosba]','Nbins',[40,40]);
    [X,Y] = meshgrid(c{1},c{2});
    pcolor(X,Y,N');
    set(gca,'ColorScale','log')
    colorbar
    xlabel('m_1_2^2');
    ylabel('cos(\alpha-\beta)');
    
    % hist for tanb vs cosba plane 
    subplot(2,3,6);
    [N,c] = hist3([tanb;m122]','Nbins',[40,40]);
    [X,Y] = meshgrid(c{1},c{2});
    pcolor(X,Y,N');
    set(gca,'ColorScale','log')
    colorbar
    xlabel('tan(\beta)');
    ylabel('m_1_2^2');

end

function massDiffPlots(p,avail)

    mH2 = p.mH2; mH2(mH2<0) = nan;
    mA2 = p.mA2; mA2(mA2<0) = nan;
    mHp2 = p.mHp2; mHp2(mHp2<0) = nan;

    mHA = sqrt(mH2) - sqrt(mA2); mHA(mHA < -500 | mHA > 500) = nan;
    mHHp = sqrt(mH2) - sqrt(mHp2); mHHp(mHHp < -500 | mHHp > 500) = nan;
    mAHp = sqrt(mA2) - sqrt(mHp2); mAHp(mAHp < -500 | mAHp > 500) = nan;
    
    mHA(avail == 0) = nan; mHHp(avail == 0) = nan; mAHp(avail == 0) = nan;

    lam1 = p.lam1; lam2 = p.lam2; lam3 = p.lam3; lam4 = p.lam4; lam5 = p.lam5;
    lam1(avail==0) = nan;
    lam2(avail==0) = nan;
    lam3(avail==0) = nan;
    lam4(avail==0) = nan;
    lam5(avail==0) = nan;

    R1 = 4*pi;
    R2 = -4*pi;
    R3 = R2/10;

    figure
    makeplot2(mHA,lam1,'m_H-m_A','\lambda_1',1,1,3,5,-500,500,R3,R1,20);
    makeplot2(mHA,lam2,'m_H-m_A','\lambda_2',1,2,3,5,-500,500,R3,R1,20);
    makeplot2(mHA,lam3,'m_H-m_A','\lambda_3',1,3,3,5,-500,500,R2,R1,20);
    makeplot2(mHA,lam4,'m_H-m_A','\lambda_4',1,4,3,5,-500,500,R2,R1,20);
    makeplot2(mHA,lam5,'m_H-m_A','\lambda_5',1,5,3,5,-500,500,R2,R1,20);

    makeplot2(mHHp,lam1,'m_H-m_H_+','\lambda_1',2,1,3,5,-500,500,R3,R1,20);
    makeplot2(mHHp,lam2,'m_H-m_H_+','\lambda_2',2,2,3,5,-500,500,R3,R1,20);
    makeplot2(mHHp,lam3,'m_H-m_H_+','\lambda_3',2,3,3,5,-500,500,R2,R1,20);
    makeplot2(mHHp,lam4,'m_H-m_H_+','\lambda_4',2,4,3,5,-500,500,R2,R1,20);
    makeplot2(mHHp,lam5,'m_H-m_H_+','\lambda_5',2,5,3,5,-500,500,R2,R1,20);

    makeplot2(mAHp,lam1,'m_A-m_H_+','\lambda_1',3,1,3,5,-500,500,R3,R1,20);
    makeplot2(mAHp,lam2,'m_A-m_H_+','\lambda_2',3,2,3,5,-500,500,R3,R1,20);
    makeplot2(mAHp,lam3,'m_A-m_H_+','\lambda_3',3,3,3,5,-500,500,R2,R1,20);
    makeplot2(mAHp,lam4,'m_A-m_H_+','\lambda_4',3,4,3,5,-500,500,R2,R1,20);
    makeplot2(mAHp,lam5,'m_A-m_H_+','\lambda_5',3,5,3,5,-500,500,R2,R1,20);
    

end

function lambdaSamplingDensityPlots(p,avail)

    % apply constraints by making points invisible
    lam1 = p.lam1(avail); lam2 = p.lam2(avail); lam3 = p.lam3(avail); lam4 = p.lam4(avail); lam5 = p.lam5(avail);
    nKeep = min(length(lam1),500000);
    lam1 = lam1(1:nKeep); lam2 = lam2(1:nKeep); lam3 = lam3(1:nKeep); lam4 = lam4(1:nKeep); lam5 = lam5(1:nKeep); 

    % upper and lower limits of plots
    H = +4*pi*3;
    L = -1*H*3;
    LL = L/10;
    
    % create lambda sampling density plots
    figure
    makePlot3(lam1,'\lambda_1',1,1,5,5,LL,H,20);
    makePlot3(lam2,'\lambda_2',2,2,5,5,LL,H,20);
    makePlot3(lam3,'\lambda_3',3,3,5,5,L,H,20);
    makePlot3(lam4,'\lambda_4',4,4,5,5,L,H,20);
    makePlot3(lam5,'\lambda_5',5,5,5,5,L,H,20);
    
    makeplot2(lam1,lam2,'\lambda_1','\lambda_2',1,2,5,5,LL,H,LL,H,20);
    makeplot2(lam1,lam3,'\lambda_1','\lambda_3',1,3,5,5,LL,H,L,H,20);
    makeplot2(lam1,lam4,'\lambda_1','\lambda_4',1,4,5,5,LL,H,L,H,20);
    makeplot2(lam1,lam5,'\lambda_1','\lambda_5',1,5,5,5,LL,H,L,H,20);
    makeplot2(lam2,lam3,'\lambda_2','\lambda_3',2,3,5,5,LL,H,L,H,20);
    makeplot2(lam2,lam4,'\lambda_2','\lambda_4',2,4,5,5,LL,H,L,H,20);
    makeplot2(lam2,lam5,'\lambda_2','\lambda_5',2,5,5,5,LL,H,L,H,20);
    makeplot2(lam3,lam4,'\lambda_3','\lambda_4',3,4,5,5,L,H,L,H,20);
    makeplot2(lam3,lam5,'\lambda_3','\lambda_5',3,5,5,5,L,H,L,H,20);
    makeplot2(lam4,lam5,'\lambda_4','\lambda_5',4,5,5,5,L,H,L,H,20);

end

function massSamplingDensityPlots(p,avail)

    % apply constraints by making points invisible
    ltanb = log10(p.tanb(avail==1)); m12 = p.m122(avail==1); mHp = p.mHp2(avail==1); mA = p.mA2(avail==1); mH = p.mH2(avail==1);
    % ltanb(avail==0) = nan; m12(avail==0) = nan; mHp(avail==0) = nan; mA(avail==0) = nan; mH(avail==0) = nan;

    % ignore imaginary masses
    m12(m12 < 0) = nan; mHp(mHp < 0) = nan; mA(mA < 0) = nan; mH(mH < 0) = nan;
    
    m12 = sqrt(m12); mHp = sqrt(mHp); mA = sqrt(mA); mH = sqrt(mH);
    % m12 = m12.^2;
    
    perInvalid = max(0.01,length(m12(isnan(m12))) / length(m12));
    nKeep = min(length(m12),round(500000/perInvalid));
    mH = mH(1:nKeep); mA = mA(1:nKeep); mHp = mHp(1:nKeep); m12 = m12(1:nKeep); ltanb = ltanb(1:nKeep); 
    
    % plot limits
    H = 10000;
    H2 = 10000;
    
    figure
    makePlot3(mH,'m_H',1,1,5,5,0,H,20);
    makePlot3(mA,'m_A',2,2,5,5,0,H,20);
    makePlot3(mHp,'m_H^+',3,3,5,5,0,H,20);
    makePlot3(m12,'sqrt(m_1_2^2)',4,4,5,5,0,H2,20);
    makePlot3(ltanb,'log(tan(\beta))',5,5,5,5,0,2.0,20);

    makeplot2(mH,mA,'m_H','m_A',1,2,5,5,0,H,0,H,20);
    makeplot2(mH,mHp,'m_H','m_H^+',1,3,5,5,0,H,0,H,20);
    makeplot2(mH,m12,'m_H','sqrt(m_1_2^2)',1,4,5,5,0,H,0,H2,20);
    makeplot2(mH,ltanb,'m_H','log(tan(\beta))',1,5,5,5,0,H,0,2.0,20);
    makeplot2(mA,mHp,'m_A','m_H^+',2,3,5,5,0,H,0,H,20);
    makeplot2(mA,m12,'m_A','sqrt(m_1_2^2)',2,4,5,5,0,H,0,H2,20);
    makeplot2(mA,ltanb,'m_A','log(tan(\beta))',2,5,5,5,0,H,0,2.0,20);
    makeplot2(mHp,m12,'m_H^+','sqrt(m_1_2^2)',3,4,5,5,0,H,0,H2,20);
    makeplot2(mHp,ltanb,'m_H^+','log(tan(\beta))',3,5,5,5,0,H,0,2.0,20);
    makeplot2(m12,ltanb,'sqrt(m_1_2^2)','log(tan(\beta))',4,5,5,5,0,H2,0,2.0,20);

end

% --- plot helper function ---

function makePlot3(A,Alabel,x,y,xWidth,yWidth,xMin,xMax,nBins)

    % A(A<xMin | A>xMax) = nan;
    A = A(A>=xMin & A<=xMax);
    
    subplot(yWidth,xWidth,(x-1)+(y-1)*xWidth+1);
    edges = linspace(xMin,xMax,nBins);
    % N = min(1,N); % turn into on/off avaliability plot 
    % m = max(N(:)); N(N==0) = -m; % compromise plot
    histogram(A,edges,'Normalization','probability')
    
    % set(gca,'YScale','log') % also good for avaliability plots

    if y-1 == yWidth-1, 
        xlabel(Alabel); 
    end
    if x-1 == 0, 
        ylabel(Alabel); 
    end

    % ylim([0.0001,1.0]);

end

function makeplot(A,B,Alabel,Blabel,x,y,xWidth,yWidth)

    % A = [A,-4*pi,-4*pi,4*pi,4*pi];
    % B = [B,-4*pi,-4*pi,4*pi,4*pi];

    subplot(yWidth,xWidth,(x-1)+(y-1)*xWidth+1)
    [N,c] = hist3([A;B]','Nbins',[20,20]);
    [X,Y] = meshgrid(c{1},c{2});
    % N = min(1,N); % turn into on/off avaliability plot 
    % m = max(N(:)); N(N==0) = -m; % compromise
    pcolor(X,Y,N');
    % set(gca,'ColorScale','log') % also good for avaliability plots

    if y-1 == yWidth-1, 
        xlabel(Alabel); 
    end
    if x-1 == 0, 
        ylabel(Blabel); 
    end

end

function makeplot2(A,B,Alabel,Blabel,x,y,xWidth,yWidth, xMin, xMax, yMin, yMax, nBins)

    A(A<xMin | A>xMax) = nan;
    B(B<yMin | B>yMax) = nan;
    
    subplot(yWidth,xWidth,(x-1)+(y-1)*xWidth+1);
    ranges = dataRange(xMin,xMax,yMin,yMax,nBins);
    [N,c] = hist3([A;B]','Edges',ranges);
    [X,Y] = meshgrid(c{1},c{2});
    % N = min(1,N); % turn into on/off avaliability plot 
    % m = max(N(:)); N(N==0) = -m; % compromise plot
    pcolor(X,Y,N');
    % set(gca,'ColorScale','log') % also good for avaliability plots

    if y-1 == yWidth-1, 
        xlabel(Alabel); 
    end
    if x-1 == 0, 
        ylabel(Blabel); 
    end

end

function binEdges = dataRange(xMin,xMax,yMin,yMax,nBins)

    xRange = (0:nBins)*((xMax-xMin)/nBins)+xMin;
    yRange = (0:nBins)*((yMax-yMin)/nBins)+yMin;

    binEdges = {xRange, yRange};

end

function printer(p, avail)

    ind = 1;

    for i=1:length(avail)
        if avail(i) == 1
            ind = i;
        end
    end

    % lambdas
    fprintf('lambda1: %f\n', p.lam1(ind));
    fprintf('lambda2: %f\n', p.lam2(ind));
    fprintf('lambda3: %f\n', p.lam3(ind));
    fprintf('lambda4: %f\n', p.lam4(ind));
    fprintf('lambda5: %f\n', p.lam5(ind));
    % fprintf('lambda6: %f\n', p.lam6(ind));
    % fprintf('lambda7: %f\n', p.lam7(ind));

    % masses, m122, m112, m222
    fprintf('m122: %f\n', p.m122(ind));
    fprintf('m112: %f\n', p.m112(ind));
    fprintf('m222: %f\n', p.m222(ind));
    fprintf('mH: %f\n', sqrt(p.mH2(ind)));
    fprintf('mA: %f\n', sqrt(p.mA2(ind)));
    fprintf('mh: %f\n', sqrt(p.mh2(ind)));
    fprintf('mHp: %f\n', sqrt(p.mHp2(ind)));
    
    % angles
    fprintf('beta: %f\n', p.beta(ind));
    fprintf('alpha: %f\n', p.alpha(ind));
    fprintf('tanb: %f\n', p.tanb(ind));
    fprintf('cosba: %f\n', p.cosba(ind));
    fprintf('sinba: %f\n', p.sinba(ind));
    
    % VEVs
    fprintf('v1: %f\n', p.v1(ind));
    fprintf('v2: %f\n', p.v2(ind));

end

function p = calcLoopCouplings(p,m)
    
    % top - bottom - charm - tau
    Nc = [3,     3,      3,     1]'; % colour factors
    Q = [+2/3,  -1/3,  +2/3,  -1]'; % charges
    m2 = ([173.07, 4.180, 1.275, 1.777]').^2; % masses
    mHp2 = p.mHp2;
    mW2 = (80.39)^2; % mass of W boson squared
    mZ2 = (91.19)^2; % mass of Z boson squared

    T = p.mh2 ./ (4*m2); TW = p.mh2 ./ (4*mW2); THp = p.mh2 ./ (4*mHp2);
    C = sqrt([m.k_htt2, m.k_hbb2, m.k_hcc2, m.k_htautau2]'); CWW = sqrt(m.k_hWW); 
    ChHpHp = 0.0;
    
    % get_v(container) * ((q[k][1]).real() * container.he->get(Par::mass1,"Lambda_3") + (q[k][2] * (double)sgn(container.he->get(Par::mass1,"Lambda_6")) * container.he->get(Par::mass1,"Lambda_7")).real());
    
    % v*(q[

    f = @(x) (x<=1).*(asin(sqrt(x)).^2) + ...
              (x>1).*(-1/4*(  log( (1+sqrt(1-x.^-1))/(1-sqrt(1-x.^-1)) ) - 1i*pi  ));

    A12 = @(x) 2*x.^(-2).*(x+(x-1).*f(x));
    A1 = @(x) -1*x.^(-2).*(2*x.^2+3*x+3*(2*x-1).*f(x));
    A0 = @(x) x.^(-2).*(f(x)-x);

    tmpA = A12(T);
    p.k_hgg2 = sum(C(1:3,:).*tmpA(1:3),1) ./ sum(tmpA(1:3),1);

    tmpC = A1(TW);
    tmpB = CWW.*tmpC-ChHpHp.*A0(THp);
    p.k_hgaga2 = (sum( Nc.*Q.^2.*C.*tmpA,1)+tmpB) ./ (sum(Nc.*Q.^2.*tmpA,1)+tmpC);

    % see ref 128,41,129
    
    p.k_hZga2 = 0;
    
    I3 = [0,0,0,0];
    Lam = p.mZ2 ./ (4*m2); LamW = p.mZ2 ./ (4*mW2); LamHp = p.mZ2 ./ (4*mHp2);
    
    I1 = @(a,b) 0.0;
    I2 = @(a,b) 0.0;
    g = @(x) 0.0;
    A12Z = @(a,b) 0.0;
    A1Z = @(a,b) 0.0;
    A0Z = 0.0;
    
    tmpA = A12Z(T,Lam);
    tmpC = A1(TW,LamW);
    tmpB = CWW.*tmpC-ChHpHp.*A0(THp,WHp);
    
    % p.k_Zga2 = sum(C.*tmpA+tmpB,1) / sum(tmpA+tmpC,1);
    
    p.k_hgg2 = abs(p.k_hgg2).^2;
    p.k_hgaga2 = abs(p.k_hgaga2).^2;
    p.k_hZga2 = abs(p.k_hZga2).^2;

end


% 2HDM sampling density tests (parameter samplers)
% by A.S. Woodock
% JUN/2022
% License: GPL2
% 
% file contents: parameter samplers (and spectrum generation)
% for the generic, physical and Higgs bases at tree-level

% unfortunately MATLAB doesn't allow multiple functions 
% in one file so we create a dummy class
classdef Sampler
methods(Static)

% INPUTS: p = the parameter struct,
%         nPoints = number of samples

% OUTPUTS: p = the full set of parameters

function r = randAB(A,B,n) 
    r = A + (B-A)*rand(1,n);
end

function [y, abc] = cubicspline(x, xj, fj)

    % Document: name, date, description, function inputs/outputs
    xj = xj(:);
    fj = fj(:);
    x = x(:);

    N = length(xj);
    A = [ones(N,1) xj abs(xj - xj(2:N-1)').^1];
    abc = A\fj;
    y = abc(1) + abc(2)*x + abs(x - xj(2:N-1)').^1*abc(3:N);

end

% Generic (aka coupling) basis
function p = generic(p)

    % Lambda range
    lam_max = 4*pi;

    % basis params
    p.lam1 = Sampler.randAB(-lam_max/10,lam_max,p.nPoints);
    p.lam2 = Sampler.randAB(-lam_max/10,lam_max,p.nPoints);
    p.lam3 = Sampler.randAB(-lam_max,lam_max,p.nPoints);
    p.lam4 = Sampler.randAB(-lam_max,lam_max,p.nPoints);
    p.lam5 = Sampler.randAB(-lam_max,lam_max,p.nPoints);
    p.lam6 = 0;
    p.lam7 = 0;
    p.m122  = Sampler.randAB(-1e6,1e6,p.nPoints);
    p.tanb = Sampler.randAB(0.3,50,p.nPoints);

    % derived parameters
    p.lam345 = p.lam3+p.lam4+p.lam5;
    p.beta = atan(p.tanb);
    p.sinb = sin(p.beta);
    p.cosb = cos(p.beta);
    
    % !!CHECK
    p.vsq  = p.v^2;
    p.lam345 = p.lam3+p.lam4+p.lam5;
    m122v2 = (2/p.vsq)*p.m122;
    cot2b = cot(2*p.beta);
    tan2a = (m122v2 - p.lam345.*sin(2*p.beta)) ./ (m122v2.*cot2b-(p.lam1.*(cos(p.beta).^2)-p.lam2.*(sin(p.beta).^2)));
    p.alpha = atan(tan2a)/2;
    p.cosba = cos(p.beta-p.alpha);
    p.sinba = sin(p.beta-p.alpha);

    p = Sampler.extras(p);
    p = Sampler.calcMasses(p);

end

% Physical (aka mass) basis
function p = physical(p)
    
    % all masses in GeV^2

    % basis params
    p.mh2   = 125.35.^2;
    p.mH2   = Sampler.randAB(0,2000,p.nPoints).^2;
    p.mA2   = Sampler.randAB(0,2000,p.nPoints).^2;
    p.mHp2  = Sampler.randAB(0,2000,p.nPoints).^2;
    p.m122  = Sampler.randAB(-1e6,+1e6,p.nPoints);
    p.cosba = Sampler.randAB(-0.01,+0.01,p.nPoints);
    p.tanb  = Sampler.randAB(0.3,50,p.nPoints);
    p.lam6  = 0;
    p.lam7  = 0;

    % derived params
    p.beta = atan(p.tanb);
    p.sinb = sin(p.beta);
    p.cosb = cos(p.beta);
    p.alpha = p.beta - acos(p.cosba);
    p.sinba = sin(p.beta-p.alpha);
    
    % !!CHECK
    p.vsq  = p.v^2;
    cosa = cos(p.alpha);
    sina = sin(p.alpha);
    p.lam1 = (p.mH2.*cosa.^2+p.mh2.*sina.^2-p.m122.*p.tanb) ./ (p.vsq.*p.cosb.^2);
    p.lam2 = (p.mH2.*sina.^2+p.mh2.*cosa.^2-p.m122./p.tanb) ./ (p.vsq.*p.sinb.^2);
    p.lam3 = ((p.mH2-p.mh2).*cosa.*sina+2*p.mHp2.*p.cosb.*p.sinb-p.m122) ./ (p.vsq.*p.sinb.*p.cosb);
    p.lam4 = ((p.mA2-2*p.mHp2).*p.sinb.*p.cosb+p.m122)  ./ (p.vsq.*p.sinb.*p.cosb);
    p.lam5 = (p.m122-p.mA2.*p.sinb.*p.cosb)./(p.vsq.*p.sinb.*p.cosb);
    p.lam345 = p.lam3+p.lam4+p.lam5;

    p = Sampler.extras(p);

end

% Higgs basis (TODO)
function p = Higgs(p)
    
end

function p = calcMasses(p)
    
    % Higgs Masses (general - not just Z2 symmetric)
    % p.mA2 = p.m122 ./ (p.sinb.*p.cosb) - (p.vsq/2)*(2*p.lam5+p.lam6.*p.cotb+p.lam7.*p.tanb); % ---
    % p.mHp2 = p.mA2 + (0.5*p.vsq)*(p.lam5-p.lam4); % +--

    Vbar = -2*p.m122+(p.lam4+p.lam5).*p.v1.*p.v2+p.lam6.*p.v1.^2+p.lam7.*p.v2.^2;
    p.mHp2 = -(p.vsq./(2*p.v1.*p.v2)) .* Vbar;
    p.mA2 = p.mHp2 + (0.5*p.vsq)*(p.lam4-p.lam5);

    % sM112 = p.mA2.*p.sinb.^2 + p.vsq*(p.lam1.*p.cosb.^2+2*p.lam6.*p.sinb.*p.cosb+p.lam5.*p.sinb.^2); % ---
    % sM222 = p.mA2.*p.cosb.^2 + p.vsq*(p.lam2.*p.sinb.^2+2*p.lam7.*p.sinb.*p.cosb+p.lam5.*p.cosb.^2); % ---
    % sM122 = -p.mA2.*p.sinb.^2.*p.cosb.^2 + p.vsq*((p.lam3+p.lam4).*p.sinb.*p.cosb+p.lam6.*p.cosb.^2+p.lam7.*p.sinb.^2); % ---

    % tmpA = sM112 + sM222; % ---
    % tmpB = sqrt((sM112-sM222).^2 + 4*sM122.^2); % ---
    % p.mh2 = 0.5*(tmpA - tmpB); % --- !!!!!
    % p.mH2 = 0.5*(tmpA + tmpB); % ---

    % expressions for P.114(Cheng)

    tmpA = 2*p.m122./sin(2*p.beta);
    % p.mH2 = tmpA-0.5*(p.lam4+p.lam5)*p.vsq;
    % p.mA2 = tmpA-p.lam5*p.vsq;
    tmpMassA = tmpA.*cos(p.beta-p.alpha).^2 + p.vsq*(p.lam1.*sin(p.alpha).^2.*p.cosb.^2+p.lam2.*cos(p.alpha).^2.*p.sinb.^2-0.5*p.lam345.*sin(2*p.alpha).*sin(2*p.beta));
    tmpMassB = tmpA.*sin(p.beta-p.alpha).^2 + p.vsq*(p.lam1.*cos(p.alpha).^2.*p.cosb.^2+p.lam2.*sin(p.alpha).^2.*p.sinb.^2+0.5*p.lam345.*sin(2*p.alpha).*sin(2*p.beta));

    p.mh2 = min(tmpMassA,tmpMassB);
    p.mH2 = max(tmpMassA,tmpMassB);

end

function p = extras(p)

    p.cotb = 1 ./ p.tanb;
    p.v2   = p.v .* p.sinb;
    p.v1   = p.v .* p.cosb;
    p.vsq  = p.v^2;

    % calc m122, m222 using tree-level tadpole conditions
    p.m112 = p.m122.*p.tanb - 0.5*(p.lam1.*p.v1.^2 + p.lam345.*p.v2.^2);
    p.m222 = p.m122.*p.cotb - 0.5*(p.lam2.*p.v2.^2 + p.lam345.*p.v1.^2);

    % general

    % p.m112 = (real(p.m122).*p.v2 - 0.5*p.lam1.*p.v1.^3 - 0.5*p.lam345.*p.v1.*p.v2.^2 ...
    %        - 0.5*(3*real(p.lam6).*p.v1.^2.*p.v2 + real(p.lam7).*p.v2.^3)) ./ p.v1;
    % 
    % p.m112 = (real(p.m122).*p.v2 - 0.5*p.lam2.*p.v2.^3 - 0.5*p.lam345.*p.v2.*p.v1.^2 ...
    %        - 0.5*(real(p.lam6).*p.v1.^3 + 3*real(p.lam7).*p.v2.*p.v1.^2)) ./ p.v2;

end


end
end

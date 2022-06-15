% 2HDM sampling density tests (theoretical constraints)
% by A.S. Woodock
% JUN/2022
% License: GPL2

% unfortunately MATLAB doesn't allow multiple functions 
% in one file so we create a dummy class
classdef TConstraint
methods(Static)

% (VERIFIED) vacuum stability
function LL = vacuumStability(p)

    % https://arxiv.org/pdf/1303.5098.pdf (pages 6)

    sqrt_lam12 = sqrt(p.lam1.*p.lam2);

    % ERROR: i think it should be l3+l4-abs(l5) !!!!
    LL = p.lam1 >= 0 & p.lam2 >= 0 & p.lam3 + sqrt_lam12 >= 0 & ...
         p.lam3+p.lam4-abs(p.lam5) + sqrt_lam12 >= 0 ; % & ...
         % 2*abs(p.lam6+p.lam7) < 0.5*(p.lam1+p.lam2)+p.lam345;
     
end

% (VERIFIED) vacuum meta-stability
function LL = vacuumMetaStability(p)

    % https://arxiv.org/pdf/1303.5098.pdf (pages 5,7)

    k = (p.lam1./p.lam2).^(1/4);
    k2m222 = (k.^2).*p.m222;
    % root_12 = sqrt(p.lam1.*p.lam2);
    % denom = p.m112 + k2m222;
    
    % x = (4*k.*p.m122 ./ denom) .* (root_12 ./ (p.lam345 - root_12));
    % y = ((p.m112-k2m222) ./ denom) .* ((root_12+p.lam345)./(root_12-p.lam345));

    % the following condition determines if the potential has only one minimum
    % (note that for a normal vacuum there can be at most two minima)
    % R = x.^(2/3) + y.^(2/3); % TODO: double-check this
    
    % the vacuum meta-stability discriminant
    D = p.m122.*(p.m112-k2m222).*(p.tanb - k);
    
    % not (panic vacuum condition)
    LL = D>0;
    
end

% (VERIFIED) leading-order S-matrix unitarity
function LL = LOunitarity(p)

    % http://arxiv.org/pdf/hep-ph/0312374.pdf (page 5)

    % leading-order tree-level unitarity constraints
    tmp = sqrt((p.lam1 - p.lam2).^2 + 4*abs(p.lam5).^2);
    a11p_even = (1/2)*(p.lam1 + p.lam2 + tmp);
    a11m_even = (1/2)*(p.lam1 + p.lam2 - tmp);

    tmp = sqrt((p.lam1 - p.lam2).^2 + 4*p.lam4.^2);
    a01p_even = (1/2)*(p.lam1 + p.lam2 + tmp);
    a01m_even = (1/2)*(p.lam1 + p.lam2 - tmp);

    tmp = sqrt(9*(p.lam1 - p.lam2).^2 + 4*(2*p.lam3+p.lam4).^2);
    a00p_even = (1/2)*(3*(p.lam1 + p.lam2) + tmp);
    a00m_even = (1/2)*(3*(p.lam1 + p.lam2) - tmp);
    
    a11_odd = (1/2) * (p.lam3+p.lam4);
    a10_odd = (1/2) * (p.lam3-p.lam4);
    a01p_odd = (1/2) * (p.lam3+abs(p.lam5));
    a01m_odd = (1/2) * (p.lam3-abs(p.lam5));
    a00p_odd = (1/2) * (p.lam3+2*p.lam4+3*abs(p.lam5));
    a00m_odd = (1/2) * (p.lam3+2*p.lam4-3*abs(p.lam5));
    
    lim = 8*pi;
    
    LL = (abs(a11p_even) < lim & abs(a11m_even) < lim & abs(a01p_even) < lim) ...
       & (abs(a01m_even) < lim & abs(a00p_even) < lim & abs(a00m_even) < lim) ...
       & (abs(a11_odd)   < lim & abs(a10_odd)   < lim & abs(a01p_odd)  < lim) ...
       & (abs(a01m_odd)  < lim & abs(a00p_odd)  < lim & abs(a00m_odd)  < lim);
    
end

% (VERIFIED) perturbativity of the generic couplings
function LL = perturbativityGeneric(p)
    
    cutL = -4*pi;
    cutH = +4*pi;
    
    LL = p.lam1 > 0    & p.lam1 < cutH ...
       & p.lam2 > 0    & p.lam2 < cutH ...
       & p.lam3 > cutL & p.lam3 < cutH ...
       & p.lam4 > cutL & p.lam4 < cutH ...
       & p.lam5 > cutL & p.lam5 < cutH;
end

% perturbativity of all scalar couplings
% TODO

% perturbativity of the Yukawa couplings
% TODO

% (VERIFIED) mass positivity constraint
function LL = massPositivity(p)

    LL = imag(p.mh2) == 0 & imag(p.mH2) == 0 & imag(p.mHp2) == 0 ...
       & imag(p.mA2) == 0 & real(p.mh2) > 0 & real(p.mH2) > 0 ...
       & real(p.mHp2) > 0 & real(p.mA2) > 0;
   
end

end
end

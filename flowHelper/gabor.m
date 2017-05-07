function [ res ] = gabor( x, y, sigma, fx0, fy0 )
% c.f. On event based optical flow Tobias Brosch et al. 2015 and
% c.f. Diploma Thesis Tobias Brosch (with DC compensation)
    res = 2*pi/sigma^2 * (exp(2*pi*1i*(fx0*x+fy0*y)) ... % modulation
        - exp(-.5*hypot(fx0,fy0)^2*sigma^2) ... %DC compensation
        ).*exp(-2*pi^2*(x.^2+y.^2)/sigma^2); % Gaussian

    
%     2*pi/sigma^2 *exp(-.5*hypot(fx0,fy0)^2*sigma^2) % DC-compensation
end


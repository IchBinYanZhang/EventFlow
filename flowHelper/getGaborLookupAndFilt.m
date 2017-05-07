function [gaborLookup, gabIm, evenGaussX, evenMod, oddGauss, oddMod] = ...
    getGaborLookupAndFilt(theta, sigma, f0, plotTrue)
%GETGABORFILT theta in degree

    fS = round(sigma/(2*pi)*4); if ~mod(fS, 2), fS = fS + 1; end
    fShalf = floor(fS/2);
%     theta = 90; % in degree % 0 degree results in filter selective for downward movement
    fx0 = cosd(90-theta)*f0; fy0=sind(90-theta)*f0; % 0 degree is horizontal with imshow
    
    
    xs = -fShalf:fShalf; ys = -fShalf:fShalf;
    [X, Y] = meshgrid(xs, ys); gabIm = gabor(X, Y, sigma, fx0, fy0);
    [X, Y] = meshgrid(-xs, -ys); gaborLookup = gabor(X, Y, sigma, fx0, fy0);
    
    
    
    %% Test for separable construction
    % c.f. neuromorphic coding project under cpe/tobi
    evenGaussX = []; evenMod = []; oddGauss = []; oddMod = [];
    if fx0 == 0 || fy0 == 0
        DCcompFac = exp(-.5*hypot(fx0,fy0)^2*sigma^2);
        fac = 2*pi/sigma^2;
        %fx0 = 0;
        
        evenGaussX = exp(-2*pi^2*xs.^2/sigma^2);
        evenMod = cos(2*pi*f0*ys).*exp(-2*pi^2*ys.^2/sigma^2);
        evenGaussY = DCcompFac.*exp(-2*pi^2*ys.^2/sigma^2); % DC comp
        
        evenMod = evenMod - evenGaussY; % DC comp
        
        % real DC comp
        even = fac*( ...
            repmat(evenGaussX, numel(ys), 1) .* repmat(evenMod', 1, numel(xs)));
        
        
        oddGauss = exp(-2*pi^2*(xs.^2)/sigma^2);
        oddMod = sin(2*pi*(f0*ys)).*exp(-2*pi^2*(ys.^2)/sigma^2);
        
        odd = fac * repmat(oddGauss, numel(ys), 1) .* repmat(oddMod', 1, numel(xs));
        if fx0 ~= 0, even = even'; end
        if fx0 ~= 0, odd = odd'; end
        
        tmp = imag(gabIm) - odd; if max(abs(tmp(:))) > 1E-10, error('too much difference'); end
        tmp = real(gabIm) - even; if max(abs(tmp(:))) > 1E-10, error('too much difference'); end
        
%         gabIm = evenApprox + 1i * odd;
    end
%     evenGaussX
%     evenMod
%     evenGaussY
    
%     oddGauss
%     oddMod
    
    %%
    
    
%     gabIm = gabIm - mean(real(gabIm(:)));
%     gaborLookup = gaborLookup - mean(real(gaborLookup(:)));
    
%     sum(real(gabIm(:))), sum(imag(gabIm(:))),

    if plotTrue
        % show Gabor filter
        figure(30);
        subplot 221; imshow(real(gabIm), [], 'InitialMagnification', 1000); colorbar; title('Gabor even');
        subplot 222; imshow(imag(gabIm), [], 'InitialMagnification', 1000); colorbar; title('Gabor odd');
        
        % Gabor lookup (point symmetric to original Gabor
        subplot 223; imshow(real(gaborLookup), [], 'InitialMagnification', 1000); colorbar; title('Gabor lookup even');
        subplot 224; imshow(imag(gaborLookup), [], 'InitialMagnification', 1000); colorbar; title('Gabor lookup odd');
        
        drawnow;
    end
end


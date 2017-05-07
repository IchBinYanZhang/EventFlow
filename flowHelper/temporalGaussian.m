function [ res ] = temporalGaussian( t, sigma, mu )
%TEMPORALGAUSSIAN processes t in mus i.e. 1s = 1E6
    res = exp((t*1E-6-mu).^2/(-2*sigma^2));
end


function [ res ] = tMono( t, sigMo, muMo )
%TMONO processes t in mus i.e. 1s = 1E6
    res = temporalGaussian(t, sigMo, muMo);
end


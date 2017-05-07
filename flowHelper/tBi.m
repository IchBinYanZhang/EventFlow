function res = tBi( t, s1, s2, sigBi1, muBi1, sigBi2, muBi2 )
%TBI processes t in mus i.e. 1s = 1E6
    res = -s1*temporalGaussian(t, sigBi1, muBi1) + s2*temporalGaussian(t, sigBi2, muBi2);
end


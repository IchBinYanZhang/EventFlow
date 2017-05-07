function [ filtResp ] = getFiltRespFromTempBuff( ...
    curTempFilteredBuffRespBi, curTempFilteredBuffRespMono, gabor, borderCond)
%GETFILTRESPFROMTEMPBUFF 
    biGabFilteredEven = imfilter(curTempFilteredBuffRespBi, real(gabor), borderCond);
    monoGabFilteredOdd = imfilter(curTempFilteredBuffRespMono, imag(gabor), borderCond);
    
    filtResp = abs(biGabFilteredEven+monoGabFilteredOdd) - ...
        abs(biGabFilteredEven-monoGabFilteredOdd);
end


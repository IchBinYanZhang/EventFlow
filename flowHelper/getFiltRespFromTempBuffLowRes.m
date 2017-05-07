function [ filtResp ] = getFiltRespFromTempBuffLowRes( ...
    curTempFilteredBuffRespBi, curTempFilteredBuffRespMono, gabor, borderCond)
%GETFILTRESPFROMTEMPBUFF 
    respBi = imfilter(imfilter(curTempFilteredBuffRespBi, [1 4 6 4 1]/16), [1 4 6 4 1]'/16);
    respBi = respBi(1:2:end, 1:2:end);
%     respBi = imfilter(imfilter(respBi, [1, 2, 1]/4), [1; 2; 1]/4);
%     respBi = respBi(1:2:end, 1:2:end);
    
    
    respMo = imfilter(imfilter(curTempFilteredBuffRespMono, [1 4 6 4 1]/16), [1 4 6 4 1]'/16);
    respMo = respMo(1:2:end, 1:2:end);
%     respMo = imfilter(imfilter(respMo, [1, 2, 1]/4), [1; 2; 1]/4);
%     respMo = respMo(1:2:end, 1:2:end);
    
    

    biGabFilteredEven = imfilter(respBi, real(gabor), borderCond);
    monoGabFilteredOdd = imfilter(respMo, imag(gabor), borderCond);
    
    filtResp = abs(biGabFilteredEven+monoGabFilteredOdd) - ...
        abs(biGabFilteredEven-monoGabFilteredOdd);
end


close all;
clear;
clc;

restoredefaultpath;
addpath('flowHelper');
addpath('loadData');

baseOutDir = 'output';
baseDir = '/Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/eventBasedWeizman/weizmannEventCropped';
personNum = 2;
% curSeqName = sprintf('person%02d_bend', personNum);
curSeqName = sprintf('person%02d_jack', personNum);
% % curSeqName = sprintf('person%02d_jump', personNum);
% % curSeqName = sprintf('person%02d_pjump', personNum);
% % curSeqName = sprintf('person%02d_run', personNum);
% curSeqName = sprintf('person%02d_side', personNum);
% % curSeqName = sprintf('person%02d_skip', personNum);
% curSeqName = sprintf('person%02d_walk', personNum);
% % curSeqName = sprintf('person%02d_wave01', personNum);
% % curSeqName = sprintf('person%02d_wave02', personNum);
% 
% 
% 
% baseDir = '';
% curSeqName = 'person03_walk';
S = load(fullfile(baseDir, [curSeqName, '.mat']));
x_coords = S.x_coords; y_coords = S.y_coords;
ts = S.ts; on_offs = S.on_offs; clear S;

% % add noise to the event time stamps
% [ts, inds] = sort(ts + randn(size(ts))*20E3); % 40E3, i.e. 40ms
% y_coords = y_coords(inds); x_coords = x_coords(inds); on_offs = on_offs(inds);

curSeqName = 'freeway.aedat';
[allAddr,allTs]=loadaerdat(curSeqName,2000000);
[x_coords, y_coords, ts, on_offs] = getDataTobi(allTs, allAddr); clear allTs allAddr;

x_coords = x_coords - min(x_coords(:)) + 1;
y_coords = y_coords - min(y_coords(:)) + 1;
maxY = max(y_coords(:));
maxX = max(x_coords(:));



borderCond = 0;
%% Play sequence and 3D plot
% %     dtPlot = 0.05E6; % mu s
% dtPlot = 50E3; % mu s
% 
% figure(10);
% for curT = 1:dtPlot:ts(end)-dtPlot-1
%     im = zeros(128, 128, 3);
%     fromPlot = find(ts > curT, 1, 'first');
%     toPlot = find(ts > curT + dtPlot, 1, 'first');
%     for k = fromPlot:toPlot
%         if on_offs(k)
%             im(y_coords(k), x_coords(k), 1) = im(y_coords(k), x_coords(k), 1) + (toPlot-k)/(toPlot-fromPlot);
%         else
%             im(y_coords(k), x_coords(k), 3) = im(y_coords(k), x_coords(k), 3) + (toPlot-k)/(toPlot-fromPlot);
%         end
%     end
%     imshow(im, 'InitialMagnification', 500);
% %     ylim([35, 119]);
%     title(sprintf('%s t=%4.2fs', curSeqName, curT*1E-6), 'Interpreter', 'none');
%     drawnow; %pause;
% end
% 
% onInds = logical(on_offs == 1);
% offInds = logical(on_offs == 0);
% 
% figure(20);
% plot3(x_coords(onInds), ts(onInds), -y_coords(onInds), '.r', ...
%     x_coords(offInds), ts(offInds), -y_coords(offInds), '.b')
% xlabel('x'); ylabel('t'); zlabel('y');
% ylim([1E6, 3E6]);

%% Parameters

muBi1 = 0.015; %0.05;
sigma = 23; %23; 
f0 = 0.06; %0.06; 


s1 = 0.5; s2 = 0.75;
muMo = 0.2*muBi1*(1+sqrt(36+10*log(s1/s2)));
muBi2 = 2*muBi1;
sigMo = muMo/3.; sigBi1 = muBi1/3.;
sigBi2 = 0.5*muBi1;

tFiltBound = round((muBi2+3*sigBi2)*1E6); %0.2E6; % adjust if parameters change!


sigmaPlotBlurr = 3;
HmotBlurr = fspecial('Gaussian', 4*sigmaPlotBlurr+1, sigmaPlotBlurr);

%% Gabor and temporal filters

theta = 90;
plotTrue = false;
[~, gabor0] = getGaborLookupAndFilt(0, sigma, f0, plotTrue);
[~, gabor90] = getGaborLookupAndFilt(90, sigma, f0, plotTrue);
% fS = size(gabor0, 1); fShalf = floor(fS/2);


%% subplot 222; % show color wheel
% [X, Y] = meshgrid(-20:20, -20:20);
% imshow(getColorWheelCoded(Y, X, 20));
fprintf('Determining flow\n');


dt = 13.4E3; % in mus % sample rate of temporal response function
bufferTs = 0:dt:tFiltBound;


valsForBufferTmonoTemplate = tMono(bufferTs, sigMo, muMo);
valsForBufferTbiTemplate = tBi(bufferTs, s1, s2, sigBi1, muBi1, sigBi2, muBi2);

figure(10)
subplot 411;
% plotTs = linspace(0, tFiltBound, 6);
plot(bufferTs, valsForBufferTmonoTemplate, '-og', ...
    bufferTs, valsForBufferTbiTemplate, '--ob', ...
    'LineWidth', 2);
xlabel('t [ms]')
title(sprintf('Temporal filters (muBi1=%6.3f)', muBi1));

subplot 412;
plot(real(gabor90(round(size(gabor0,1)/2), :)));

%%
figure(11);

% ring buffer idea: buffer(y, x, nextPlotBufferInd) determines responses used at
% time nextPlotTime
% -> when event arrives, determine temporal filter response and add it to
% previous entries
% -> when plot time is reached, filter with Gabors, combine etc. and set
% buffer of this time step to zero and advance the buffer ind
% nextPlotBufferInd
bufferMonoExact = zeros(maxY, maxX, numel(bufferTs));
bufferBiExact = zeros(maxY, maxX, numel(bufferTs));
bufferMonoApprox = zeros(maxY, maxX, numel(bufferTs));
bufferBiApprox = zeros(maxY, maxX, numel(bufferTs));
nextBuffTime = dt;
nextBufferInd = 1;

dtPlot = 40E3; % in mus (for 40E3mus you get 25frames/s)
nextPlotTime = dtPlot;


respPRange = [-1, 1]; %[-.15, .15]; %[-.5, .5];
maxPlotSpeed = 1; %1/10; % velocities plot scaling

% subplot 221; % show color wheel
% [X, Y] = meshgrid(-20:20, -20:20);
% imshow(getColorWheelCoded(Y, X, 20));

for tInd = 1:numel(ts)
    curT = ts(tInd);
    curY = y_coords(tInd); curX = x_coords(tInd); curOnOff = ((on_offs(tInd)-.5)*2);
    
    %% progress by one delay (=move buffer pointer foward in ring)
    while curT > nextBuffTime,
        
        % progress by one delay (=move buffer pointer foward in ring)
        bufferMonoExact(:, :, nextBufferInd) = 0; % not needed anymore must be set to zero!
        bufferBiExact(:, :, nextBufferInd) = 0; % otherwise interferes with new values
        
        
        bufferMonoApprox(:, :, nextBufferInd) = 0; % not needed anymore must be set to zero!
        bufferBiApprox(:, :, nextBufferInd) = 0; % otherwise interferes with new values
        nextBuffTime = nextBuffTime + dt;
        nextBufferInd = mod(nextBufferInd+1-1, numel(bufferTs))+1; 
    end
    
    %% temporal filtering and storage in buffer
    % which corresponds to delayed propagation in neuron models
    %% update buffers
    offset = nextBuffTime-curT; % exact
    valsForBufferTmonoExact = curOnOff*tMono(offset+bufferTs, sigMo, muMo);
    valsForBufferTbiExact = curOnOff*tBi(offset+bufferTs, s1, s2, sigBi1, muBi1, sigBi2, muBi2);
    
    valsForBufferTmonoApprox = curOnOff*valsForBufferTmonoTemplate;
    valsForBufferTbiApprox = curOnOff*valsForBufferTbiTemplate;
    
    
    bufferMonoApprox(curY, curX, :) = ...
        squeeze(bufferMonoApprox(curY, curX, :)) + circshift(valsForBufferTmonoApprox, [0, nextBufferInd-1])';
    bufferBiApprox(curY, curX, :) = ...
        squeeze(bufferBiApprox(curY, curX, :)) + circshift(valsForBufferTbiApprox, [0, nextBufferInd-1])';
    
    bufferMonoExact(curY, curX, :) = ...
        squeeze(bufferMonoExact(curY, curX, :)) + circshift(valsForBufferTmonoExact, [0, nextBufferInd-1])';
    bufferBiExact(curY, curX, :) = ...
        squeeze(bufferBiExact(curY, curX, :)) + circshift(valsForBufferTbiExact, [0, nextBufferInd-1])';
    
    
    %% plot
    if curT > nextPlotTime
        
        nextPlotTime = nextPlotTime + dtPlot/2;
        %% determine Gabor filter response of mono and bi filtered responses
        
        curTempFilteredBuffRespMonoExact = bufferMonoExact(:, :, nextBufferInd);
        curTempFilteredBuffRespBiExact = bufferBiExact(:, :, nextBufferInd);
        
        filtResp0Exact = getFiltRespFromTempBuff( ...
            curTempFilteredBuffRespBiExact, curTempFilteredBuffRespMonoExact, gabor0, borderCond);
        filtResp90Exact = getFiltRespFromTempBuff( ...
            curTempFilteredBuffRespBiExact, curTempFilteredBuffRespMonoExact, gabor90, borderCond);
        
        
        
        curTempFilteredBuffRespMonoApprox = bufferMonoApprox(:, :, nextBufferInd);
        curTempFilteredBuffRespBiApprox = bufferBiApprox(:, :, nextBufferInd);
        
        filtResp0Approx = getFiltRespFromTempBuff( ...
            curTempFilteredBuffRespBiApprox, curTempFilteredBuffRespMonoApprox, gabor0, borderCond);
        filtResp90Approx = getFiltRespFromTempBuff( ...
            curTempFilteredBuffRespBiApprox, curTempFilteredBuffRespMonoApprox, gabor90, borderCond);
        
%         fprintf('%8.5f, %8.5f\n', min(filtResp0Approx(:)), max(filtResp0Approx(:)));
        
        %% actual plot
        % exact
        respPlot90Exact = imfilter(filtResp90Exact, HmotBlurr);
        respPlot0Exact = imfilter(filtResp0Exact, HmotBlurr);
        colorWheelCodedPlotExact = getColorWheelCoded(respPlot0Exact, respPlot90Exact, maxPlotSpeed);
        subSamp = 5; quiverScale = 10*subSamp;
        
        subplot(3,3,1);
        imshow(colorWheelCodedPlotExact); hold on;
        quiver(1:subSamp:maxX, 1:subSamp:maxY, ...
            respPlot90Exact(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
            respPlot0Exact(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
            'Autoscale', 'off');
        axis ij equal; xlim([1, maxX]); ylim([1, maxY]);
        title(sprintf('Mot vecs (sigma=%4.1f)', sigmaPlotBlurr)); hold off;
        
        
        subplot(3, 3, 4);
        imshow(respPlot0Exact, respPRange);
        title('Vertical motion response');
        
        subplot(3, 3, 7);
        imshow(respPlot90Exact, respPRange);
        title('Horizontal motion response');
        
        
        
        % approx
        respPlot90Approx = imfilter(filtResp90Approx, HmotBlurr);
        respPlot0Approx = imfilter(filtResp0Approx, HmotBlurr);
        colorWheelCodedPlotApprox = getColorWheelCoded(respPlot0Approx, respPlot90Approx, maxPlotSpeed);
        
        subplot(3,3,2);
        imshow(colorWheelCodedPlotApprox); hold on;
        quiver(1:subSamp:maxX, 1:subSamp:maxY, ...
            respPlot90Approx(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
            respPlot0Approx(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
            'Autoscale', 'off');
        axis ij equal; xlim([1, maxX]); ylim([1, maxY]);
        title(sprintf('Mot vecs (sigma=%4.1f)', sigmaPlotBlurr)); hold off;
        
        
        subplot(3, 3, 5);
        imshow(respPlot0Approx, respPRange);
        title('Vertical motion response');
        
        subplot(3, 3, 8);
        imshow(respPlot90Approx, respPRange);
        title('Horizontal motion response');
        
        
        
        % diff
        respPlot90Diff = respPlot90Exact-respPlot90Approx;
        respPlot0Diff = respPlot0Exact-respPlot0Approx;
        colorWheelCodedPlotDiff = getColorWheelCoded(respPlot0Diff, respPlot90Diff, maxPlotSpeed);
        
        subplot(3,3,3);
        imshow(colorWheelCodedPlotDiff); hold on;
        quiver(1:subSamp:maxX, 1:subSamp:maxY, ...
            respPlot90Diff(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
            respPlot0Diff(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
            'Autoscale', 'off');
        axis ij equal; xlim([1, maxX]); ylim([1, maxY]);
        title(sprintf('Mot vecs (sigma=%4.1f)', sigmaPlotBlurr)); hold off;
        
        
        subplot(3, 3, 6);
        imshow(respPlot0Diff, respPRange);
        title('Vertical motion response');
        
        subplot(3, 3, 9);
        imshow(respPlot90Diff, respPRange);
        title('Horizontal motion response');
        
        drawnow; %pause;
%         print('-dpdf', '-r600', fullfile(baseOutDir, 'imsComp', sprintf('im%07d.pdf', tInd)));
    end
end





















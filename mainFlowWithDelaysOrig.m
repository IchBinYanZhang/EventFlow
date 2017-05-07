close all;
clear;
clc;

restoredefaultpath;
addpath('flowHelper');
addpath('loadData');
vn = @(x) inputname(1);


baseOutDir = 'output';
outDir = fullfile(baseOutDir, 'flowCPU');
% % baseDir = '/Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/eventBasedWeizman/weizmannEventCropped';
% personNum = 2;
% % curSeqName = sprintf('person%02d_bend', personNum);
% curSeqName = sprintf('person%02d_jack', personNum);
% % % curSeqName = sprintf('person%02d_jump', personNum);
% % % curSeqName = sprintf('person%02d_pjump', personNum);
% % % curSeqName = sprintf('person%02d_run', personNum);
% % curSeqName = sprintf('person%02d_side', personNum);
% % % curSeqName = sprintf('person%02d_skip', personNum);
% % curSeqName = sprintf('person%02d_walk', personNum);
% % % curSeqName = sprintf('person%02d_wave01', personNum);
% % % curSeqName = sprintf('person%02d_wave02', personNum);
% % 
% % 
% % 
% % % baseDir = '';
% % % curSeqName = 'person03_walk';
% % S = load(fullfile(baseDir, [curSeqName, '.mat']));
% % x_coords = S.x_coords; y_coords = S.y_coords;
% % ts = S.ts; on_offs = S.on_offs; clear S;
% 
% % % add noise to the event time stamps
% % [ts, inds] = sort(ts + randn(size(ts))*20E3); % 40E3, i.e. 40ms
% % y_coords = y_coords(inds); x_coords = x_coords(inds); on_offs = on_offs(inds);

preciseFold = fullfile('../dvs128/precise');
curSeqName = fullfile(preciseFold, 'halfCircleNotBlurred24V.aedat'); maxNumEvents = 2E6; fInd = 1;
% curSeqName = fullfile('stimuli', 'freeway.aedat'); maxNumEvents = 2E6; fInd = 1;
curSeqName = fullfile('stimuli', 'straight_slanted_bar2.aedat'); maxNumEvents = .045E6; fInd = 0.01E6;
[allAddr,allTs]=loadaerdat(curSeqName, maxNumEvents);
[x_coords, y_coords, ts, on_offs] = getDataTobi(allTs, allAddr); clear allTs allAddr;

x_coords = x_coords(fInd:end);
y_coords = y_coords(fInd:end);
ts = ts(fInd:end);
ts = ts - ts(1);
on_offs = on_offs(fInd:end);

% x_coords = x_coords - min(x_coords(:)) + 1;
% y_coords = y_coords - min(y_coords(:)) + 1;
% maxY = max(y_coords(:));
% maxX = max(x_coords(:));
maxY = 128; maxX = 128;



borderCond = 0;
%% Play sequence and 3D plot
% %     dtPlot = 0.05E6; % mu s
% dtPlot = 5E3; % mu s
% 
% figure(10);
% for curT = 1:dtPlot:ts(end)-dtPlot-1
%     im = zeros(128, 128, 3);
%     fromPlot = find(ts > curT, 1, 'first');
%     toPlot = find(ts > curT + dtPlot, 1, 'first');
%     for k = fromPlot:toPlot
%         if on_offs(k)
%             im(y_coords(k), x_coords(k), 1) = im(y_coords(k), x_coords(k), 1) + (k-fromPlot)/(toPlot-fromPlot);
%         else
%             im(y_coords(k), x_coords(k), 3) = im(y_coords(k), x_coords(k), 3) + (k-fromPlot)/(toPlot-fromPlot);
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

plotTrue = false;
[~, gabor0] = getGaborLookupAndFilt(0, sigma, f0, plotTrue);
[~, gabor90] = getGaborLookupAndFilt(90, sigma, f0, plotTrue);
% fS = size(gabor0, 1); fShalf = floor(fS/2);


%% subplot 222; % show color wheel
% [X, Y] = meshgrid(-20:20, -20:20);
% imshow(getColorWheelCoded(Y, X, 20));
fprintf('Determining flow\n');


dt = .1E3; % in mus % sample rate of temporal response function
bufferTs = 0:dt:tFiltBound;


valsForBufferTmonoTemplate = tMono(bufferTs, sigMo, muMo);
valsForBufferTbiTemplate = tBi(bufferTs, s1, s2, sigBi1, muBi1, sigBi2, muBi2);

% ring buffer idea: buffer(y, x, nextPlotBufferInd) determines responses used at
% time nextPlotTime
% -> when event arrives, determine temporal filter response and add it to
% previous entries
% -> when plot time is reached, filter with Gabors, combine etc. and set
% buffer of this time step to zero and advance the buffer ind
% nextPlotBufferInd
bufferMono = zeros(maxY, maxX, numel(bufferTs));
bufferBi = zeros(maxY, maxX, numel(bufferTs));
nextBuffTime = dt;
nextBufferInd = 1;

dtPlot = 40E3; % in mus (for 40E3mus you get 25frames/s)
nextPlotTime = dtPlot;




save(fullfile(outDir, 'params.mat'), vn(x_coords), vn(y_coords), vn(ts), vn(on_offs), ...
    vn(muBi1), vn(sigma), vn(f0), vn(s1), vn(s2), vn(muMo), vn(muBi2), ...
    vn(sigMo), vn(sigBi1), vn(sigBi2), ...
    vn(tFiltBound), vn(dt), vn(dtPlot));




respPRange = [-1, 1]; %[-.15, .15]; %[-.5, .5];
maxPlotSpeed = 1; %1/10; % velocities plot scaling

[X, Y] = meshgrid(-20:20, -20:20);
colorWheel = imresize(getColorWheelCoded(Y, X, 20), [128, 128]);
% subplot 221; % show color wheel
% imshow(colorWheel);
if ishandle(20); close(20); end; figure(20);

for tInd = 1:numel(ts)
    curT = ts(tInd);
    curY = y_coords(tInd); curX = x_coords(tInd); curOnOff = ((on_offs(tInd)-.5)*2);
    
    %% progress by one delay (=move buffer pointer foward in ring)
    while curT > nextBuffTime,
        % progress by one delay (=move buffer pointer foward in ring)
        bufferMono(:, :, nextBufferInd) = 0; % not needed anymore must be set to zero!
        bufferBi(:, :, nextBufferInd) = 0; % otherwise interferes with new values
        nextBuffTime = nextBuffTime + dt;
        nextBufferInd = mod(nextBufferInd+1-1, numel(bufferTs))+1; 
    end
        
    %% temporal filtering and storage in buffer
    % which corresponds to delayed propagation in neuron models
    %% update buffers
    offset = nextBuffTime-curT; % exact
%     offset = 0; % approximation (suitable for e.g. Truenorth)
    if offset == 0
        valsForBufferTmono = curOnOff*valsForBufferTmonoTemplate;
        valsForBufferTbi = curOnOff*valsForBufferTbiTemplate;
    else
        valsForBufferTmono = curOnOff*tMono(offset+bufferTs, sigMo, muMo);
        valsForBufferTbi = curOnOff*tBi(offset+bufferTs, s1, s2, sigBi1, muBi1, sigBi2, muBi2);
    end
    bufferMono(curY, curX, :) = ...
        squeeze(bufferMono(curY, curX, :)) + circshift(valsForBufferTmono, [0, nextBufferInd-1])';
    
    bufferBi(curY, curX, :) = ...
        squeeze(bufferBi(curY, curX, :)) + circshift(valsForBufferTbi, [0, nextBufferInd-1])';
    
    %% plot
    if curT > nextPlotTime
        %% determine Gabor filter response of mono and bi filtered responses
        
        curTempFilteredBuffRespMono = bufferMono(:, :, nextBufferInd);
        curTempFilteredBuffRespBi = bufferBi(:, :, nextBufferInd);
        
        filtResp0 = getFiltRespFromTempBuff( ...
            curTempFilteredBuffRespBi, curTempFilteredBuffRespMono, gabor0, borderCond);
        filtResp90 = getFiltRespFromTempBuff( ...
            curTempFilteredBuffRespBi, curTempFilteredBuffRespMono, gabor90, borderCond);
        
        %% Input
        im = zeros(128, 128, 3);
        fromPlot = find(ts > curT-dtPlot, 1, 'first');
        toPlot = find(ts > curT, 1, 'first');
        for k = fromPlot:toPlot
            if on_offs(k)
                im(y_coords(k), x_coords(k), 1) = im(y_coords(k), x_coords(k), 1) + (k-fromPlot)/(toPlot-fromPlot);
            else
                im(y_coords(k), x_coords(k), 3) = im(y_coords(k), x_coords(k), 3) + (k-fromPlot)/(toPlot-fromPlot);
            end
        end
        
        
        %% actual plot
        nextPlotTime = nextPlotTime + dtPlot/2;
        respPlot90 = imfilter(filtResp90, HmotBlurr);
        respPlot0 = imfilter(filtResp0, HmotBlurr);
        colorWheelCodedPlot = getColorWheelCoded(respPlot0, respPlot90, maxPlotSpeed);
        subSamp = 5; quiverScale = 10*subSamp;
        
        
        combPlot = [im, colorWheelCodedPlot; ...
            repmat((respPlot0-respPRange(1))/diff(respPRange), 1, 1, 3), ...
            repmat((respPlot90-respPRange(1))/diff(respPRange), 1, 1, 3)];
        imshow(combPlot, 'InitialMagnification', 200);
        
        %             subplot(2,2,2);
        %             imshow(colorWheelCodedPlot);
        hold on;
        quiver(128+(1:subSamp:maxX), 1:subSamp:maxY, ...
            respPlot90(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
            respPlot0(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
            'Autoscale', 'off');
        %             axis ij equal; xlim([1, maxX]); ylim([1, maxY]);
        %             title(sprintf('Mot vecs (sigma=%4.1f)', sigmaPlotBlurr));
        title(sprintf('Input, Optical flow field; vertical resp, horizontal resp (t=%dms)', ...
            nextPlotTime/1E3));
        hold off;
        %
        %
        %             subplot(2, 2, 3);
        %             imshow(respPlot0, respPRange);
        %             title('Vertical motion response');
        %
        %             subplot(2, 2, 4);
        %             imshow(respPlot90, respPRange);
        %             title('Horizontal motion response');
        
        drawnow; %pause;
        baseOutFileName = sprintf('%09d', nextPlotTime);
%         save(fullfile(outDir, sprintf('t%s.mat', baseOutFileName)), vn(filtResp0), vn(filtResp90));
%         print('-dpng', '-r300', fullfile(outDir, sprintf('im%s.png', baseOutFileName)));
    end
end





















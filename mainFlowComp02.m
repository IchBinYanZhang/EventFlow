close all;
clear;
clc;

restoredefaultpath;
addpath('flowHelper');
addpath('loadData');

baseOutDir = 'output';
% baseDir = '/Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/eventBasedWeizman/weizmannEventCropped';
% personNum = 1;
% curSeqName = sprintf('person%02d_bend', personNum);
% curSeqName = sprintf('person%02d_jack', personNum);
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
baseDir = '';
curSeqName = 'person03_walk';
S = load(fullfile(baseDir, [curSeqName, '.mat']));
x_coords = S.x_coords; y_coords = S.y_coords;
ts = S.ts; on_offs = S.on_offs; clear S;

% add noise to the event time stamps
[ts, inds] = sort(ts + randn(size(ts))*20E3); % 40E3, i.e. 40ms
y_coords = y_coords(inds); x_coords = x_coords(inds); on_offs = on_offs(inds);

% curSeqName = 'freeway.aedat';
% [allAddr,allTs]=loadaerdat(curSeqName,2000000);
% [x_coords, y_coords, ts, on_offs] = getDataTobi(allTs, allAddr); clear allTs allAddr;

x_coords = x_coords - min(x_coords(:)) + 1;
y_coords = y_coords - min(y_coords(:)) + 1;
maxY = max(y_coords(:));
maxX = max(x_coords(:));

%% Play sequence and 3D plot
% %     dtPlot = 0.05E6; % mu s
dtPlot = 50E3; % mu s
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

tFiltBound = 0.2E6; %0.2E6; % adjust if parameters change!
muBi1 = 0.05; %0.05;
sigma = 23; %23; 
f0 = 0.06; %0.06; 


s1 = 0.5; s2 = 0.75;
muMo = 0.2*muBi1*(1+sqrt(36+10*log(s1/s2)));
muBi2 = 2*muBi1;
sigMo = muMo/3.; sigBi1 = muBi1/3.;
sigBi2 = 0.5*muBi1;

sigmaPlotBlurr = 3;
HmotBlurr = fspecial('Gaussian', 4*sigmaPlotBlurr+1, sigmaPlotBlurr);

%% Gabor and temporal filters

theta = 90;
plotTrue = false;
[gaborLookup0, ~] = getGaborLookupAndFilt(0, sigma, f0, plotTrue); gaborLookup0 = single(gaborLookup0);
[gaborLookup90, ~] = getGaborLookupAndFilt(90, sigma, f0, plotTrue); gaborLookup90 = single(gaborLookup90);
gaborLookup0Real = real(gaborLookup0); gaborLookup0Imag = imag(gaborLookup0);
gaborLookup90Real = real(gaborLookup90); gaborLookup90Imag = imag(gaborLookup90);
clear gaborLookup0 gaborLookup90;
fS = size(gaborLookup0Real, 1); fShalf = floor(fS/2);

%%
% figure(40);
bufferSize = 60000;
if bufferSize > 10000; 
    fprintf('WARNING: You should have more than 8Gb RAM!!!\n'); %error('Remove if you are sure'); 
end
% if numel(ts) > 10000, error('Buffer would be too big'); end
buffer0Real = zeros(maxY, maxX, bufferSize, 'single'); buffer0Imag = zeros(maxY, maxX, bufferSize, 'single');
buffer90Real = zeros(maxY, maxX, bufferSize, 'single'); buffer90Imag = zeros(maxY, maxX, bufferSize, 'single');
resp0 = zeros(maxY, maxX); resp90 = zeros(maxY, maxX); respPlotDecay = 0.99995;

respPRange = [-.15, .15]; %[-.5, .5];
bufPRange = respPRange; %[-.5, .5];
maxPlotSpeed = 1/10; % velocities plot scaling

helperTforPlot = 0;


%% subplot 222; % show color wheel
% [X, Y] = meshgrid(-20:20, -20:20);
% imshow(getColorWheelCoded(Y, X, 20));
fprintf('Determining flow\n');

buffer_tInd1 = 1;
for tInd = 1:numel(ts)
    curY = y_coords(tInd); curX = x_coords(tInd); curT = ts(tInd); curOnOff = ((on_offs(tInd)-.5)*2);
    
    %% event evokes response in local neighborhood
    %  weight with gabor lookup and store in buffer
    
    % determine neighborhood indices in image and filter space
    yFrom = max(1, curY-fShalf); filtYFrom = yFrom - (curY-fShalf) + 1;
    yTo = min(maxY, curY+fShalf); filtYTo = fS - (curY+fShalf-yTo);
    xFrom = max(1, curX-fShalf); filtXFrom = xFrom - (curX-fShalf) + 1;
    xTo = min(maxX, curX+fShalf); filtXTo = fS - (curX+fShalf-xTo);
    buffer0Real(yFrom:yTo, xFrom:xTo, tInd-buffer_tInd1+1) ... % write spatial filter response into buffer
        = curOnOff * gaborLookup0Real(filtYFrom:filtYTo, filtXFrom:filtXTo);
    buffer0Imag(yFrom:yTo, xFrom:xTo, tInd-buffer_tInd1+1) ... % write spatial filter response into buffer
        = curOnOff * gaborLookup0Imag(filtYFrom:filtYTo, filtXFrom:filtXTo);
    
    buffer90Real(yFrom:yTo, xFrom:xTo, tInd-buffer_tInd1+1) ...
        = curOnOff * gaborLookup90Real(filtYFrom:filtYTo, filtXFrom:filtXTo);
    buffer90Imag(yFrom:yTo, xFrom:xTo, tInd-buffer_tInd1+1) ...
        = curOnOff * gaborLookup90Imag(filtYFrom:filtYTo, filtXFrom:filtXTo);
   
    
    %% determine temporal filter weights
    % first find entries in buffer at this location
    tFrom = find(ts > curT-tFiltBound, 1, 'first');
    temporalRange = (tFrom:tInd)-buffer_tInd1+1;
    bufferInds = temporalRange(...
        buffer0Real(curY, curX, temporalRange) ~= 0 | buffer0Imag(curY, curX, temporalRange) ~= 0 | ...
        buffer90Real(curY, curX, temporalRange) ~= 0 | buffer90Imag(curY, curX, temporalRange) ~= 0);
    clear temporalRange;
    relTs = curT - ts(bufferInds+buffer_tInd1-1); % relative times serve as input to temporal filter function
    
    % determine temporal filter weights
    tBiWeights = tBi(relTs, s1, s2, sigBi1, muBi1, sigBi2, muBi2);
    tMonoWeights = tMono(relTs, sigMo, muMo);
    
    %% filter orientation theta=0
    biTimesEven0 = tBiWeights * squeeze(buffer0Real(curY, curX, bufferInds));
    monoTimesOdd0 = tMonoWeights * squeeze(buffer0Imag(curY, curX, bufferInds));
    
    % filter selective for exactly one direction
%     filtResp0 =  biTimesEven0 + monoTimesOdd0; 

    % filter selective for both directions 
    % (preferred one white, other black)
    % note that addition results in filter in one direction, subtraction in
    % filter selective for other direction
    filtResp0 =  abs(biTimesEven0 + monoTimesOdd0) - abs(biTimesEven0 - monoTimesOdd0); 
    
    %% filter orientation theta=90
    biTimesEven90 = tBiWeights * squeeze(buffer90Real(curY, curX, bufferInds));
    monoTimesOdd90 = tMonoWeights * squeeze(buffer90Imag(curY, curX, bufferInds));
%     filtResp90 =  biTimesEven90 + monoTimesOdd90;
    filtResp90 =  abs(biTimesEven90 + monoTimesOdd90) - abs(biTimesEven90 - monoTimesOdd90);
    
    %% resp0 integrates activities, response fades to 0 with respPlotDecay per mu s
    if tInd > 1, 
        resp0 = resp0 * respPlotDecay^(ts(tInd)-ts(tInd-1)); 
        resp90 = resp90 * respPlotDecay^(ts(tInd)-ts(tInd-1)); 
    end
    resp0(curY, curX) = resp0(curY, curX) + filtResp0;
    resp90(curY, curX) = resp90(curY, curX) + filtResp90;
    
    
    %% shift buffer if at end
    if tInd - buffer_tInd1 + 1 == bufferSize,
        tic;
        tFrom = find(ts > curT-tFiltBound, 1, 'first'); 
        if tFrom <= buffer_tInd1, error('Buffer too small'); end
        firstIndToShiftToBeginning = tFrom-buffer_tInd1+1;
        indsToShift = firstIndToShiftToBeginning:bufferSize;
        
        nnzBuffer0Real = nnz(buffer0Real);
        buffer0RealNew = zeros(maxY, maxX, bufferSize, 'single'); 
        buffer0RealNew(:, :, 1:numel(indsToShift)) = buffer0Real(:, :, indsToShift);
        buffer0Real = buffer0RealNew; clear buffer0RealNew; % not needed anymore
        
        nnzBuffer0Imag = nnz(buffer0Real);
        buffer0ImagNew = zeros(maxY, maxX, bufferSize, 'single');
        buffer0ImagNew(:, :, 1:numel(indsToShift)) = buffer0Imag(:, :, indsToShift);    
        buffer0Imag = buffer0ImagNew; clear buffer0RealNew;
        
        nnzBuffer90Real = nnz(buffer0Real);
        buffer90RealNew = zeros(maxY, maxX, bufferSize, 'single'); 
        buffer90RealNew(:, :, 1:numel(indsToShift)) = buffer90Real(:, :, indsToShift);
        buffer90Real = buffer90RealNew; clear buffer0RealNew; % not needed anymore
        
        nnzBuffer90Imag = nnz(buffer0Real);
        buffer90ImagNew = zeros(maxY, maxX, bufferSize, 'single');
        buffer90ImagNew(:, :, 1:numel(indsToShift)) = buffer90Imag(:, :, indsToShift);
        buffer90Imag = buffer90ImagNew; clear buffer0RealNew;
        
        buffer_tInd1 = tFrom;
        fprintf('Shifted buffer in %6.3fs buffer_tInd1=%6.0f, sparsity=%7.5f\n', toc, buffer_tInd1, ...
            (nnzBuffer0Real + nnzBuffer0Imag + nnzBuffer90Real + nnzBuffer90Imag)/(4*maxY*maxX*bufferSize));
    end
    
    %% plot
    if helperTforPlot < curT
        helperTforPlot = helperTforPlot + 40E3;
        
        tFromPlot = find(ts > curT-tFiltBound, 1, 'first');
        
        bufferPlot0 = sum(buffer0Imag(:, :, (tFromPlot:tInd)-buffer_tInd1+1), 3);
        bufferPlot90 = sum(buffer90Imag(:, :, (tFromPlot:tInd)-buffer_tInd1+1), 3);
        
        % blurr responses for plot
        respPlot0 = imfilter(resp0, HmotBlurr);
        respPlot90 = imfilter(resp90, HmotBlurr);
        
%         imOut = [(bufferPlot0-bufPRange(1))/(bufPRange(2)-bufPRange(1)), ...
%             (resp0-respPRange(1))/(respPRange(2)-respPRange(1)); ...
%             (bufferPlot90-bufPRange(1))/(bufPRange(2)-bufPRange(1)), ...
%             (resp90-respPRange(1))/(respPRange(2)-respPRange(1))];
%         
%         subplot(2,2,3,'replace');
%         imshow(imOut); title(sprintf('upper: theta=0 lower: theta=90; l/r: buffer/response'));
%         
%%

        integralImagePlot = getIntegralImage(y_coords(tFromPlot:tInd), x_coords(tFromPlot:tInd), ...
            on_offs(tFromPlot:tInd), maxY, maxX);
        colorWheelCodedPlot = getColorWheelCoded(respPlot0, respPlot90, maxPlotSpeed);
        subSamp = 5; quiverScale = 30*subSamp;
        
%         subplot(2,2,1);
%         imshow(integralImagePlot);
%         title(sprintf('Event display t=%5dms (Delta t=%3dms)', round(curT*1E-4)*10, tFiltBound*1E-3));
%         
%         subplot(2,2,2);
%         imshow(colorWheelCodedPlot); hold on;
%         quiver(1:subSamp:maxX, 1:subSamp:maxY, ...
%             respPlot90(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
%             respPlot0(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
%             'Autoscale', 'off'); 
%         axis ij equal; xlim([1, maxX]); ylim([1, maxY]);
%         title(sprintf('Mot vecs (sigma=%4.1f)', sigmaPlotBlurr)); hold off;
%         
%         
%         subplot(2, 2, 3);
%         imshow(respPlot0, respPRange);
%         title('Vertical motion response');
%         
%         subplot(2, 2, 4);
%         imshow(respPlot90, respPRange);
%         title('Horizontal motion response');
%         
%         drawnow;
%         saveas(gcf, fullfile(baseOutDir, 'ims', sprintf('%s_t_%05d.pdf', curSeqName, round(curT*1E-4)*10)));
        save(fullfile(baseOutDir, 'mats', sprintf('%s_t_%05d.mat', curSeqName, round(curT*1E-4)*10)), ...
            'integralImagePlot', 'curT', 'tFiltBound', ...
            'colorWheelCodedPlot', 'subSamp', 'quiverScale', 'maxX', 'maxY', 'respPlot90', 'respPlot0', 'sigmaPlotBlurr',  ...
            'resp0', 'resp90', ...
            'respPRange', 'baseOutDir', 'curSeqName');


        fprintf('.'); 
    end
end
fprintf('\n');





%% 
% figure(50);
% t0Files = dir(fullfile(baseOutDir, 'mats', 'theta_000_t_*.mat'));
% t90Files = dir(fullfile(baseOutDir, 'mats', 'theta_090_t_*.mat'));
% 
% sigmaPlotBlurr = 3;
% H = fspecial('Gaussian', 4*sigmaPlotBlurr+1, sigmaPlotBlurr);
% respPlotPlotRange = [-.1, .1];
% 
% for k = 1:numel(t0Files)
%     S0 = load(fullfile(baseOutDir, 'mats', t0Files(k).name));
%     t0BuffPlot = S0.bufferPlot; t0Resp = S0.respPlot;
%     S90 = load(fullfile(baseOutDir, 'mats', t90Files(k).name));
%     t90BuffPlot = S90.bufferPlot; t90Resp = S90.respPlot;
%     
%     t0Resp = imfilter(t0Resp, H);
%     t90Resp = imfilter(t90Resp, H);
%     
%     subplot 321;
%     imshow(t0BuffPlot, bufferPlotRange); title('buffer theta=0');
%     subplot 322;
%     imshow(t0Resp, respPlotPlotRange); title('filter response theta=0');
%     
%     subplot 323;
%     imshow(t90BuffPlot, bufferPlotRange); title('buffer theta=90');
%     subplot 324;
%     imshow(t90Resp, respPlotPlotRange); title('filter response theta=90');
%     
%     subplot 326;
% %     motionDirs = atan2(t0Resp, t90Resp); motionDirsPlot = ind2rgb(round((motionDirs + pi)/(2*pi)*256), hsv);
% %     speeds = hypot(t0Resp, t90Resp); motionDirsPlot = motionDirsPlot .* repmat(speeds, [1, 1, 3]) * 10;
% %     imshow(motionDirsPlot);
%     subSamp = 5; quiverScale = 30*subSamp;
%     quiver(1:subSamp:maxX, 1:subSamp:maxY, t90Resp(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, t0Resp(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
%         'Autoscale', 'off'); axis ij equal; xlim([1, maxX]); ylim([1, maxY]);
%     drawnow;
% end
% 
% 












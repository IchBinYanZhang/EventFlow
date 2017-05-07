close all;
clear;
clc;


baseDir = fullfile('output', 'mats', 'person03_walk');
baseImDir = fullfile('output', 'ims');


files = dir(fullfile(baseDir, '*.mat'));
nFiles = numel(files);

for k = 1:nFiles
    load(fullfile(baseDir, files(k).name));
    [~, curSeqName, ~] = fileparts(curSeqName);
%     outFile = fullfile(baseImDir, sprintf('%s_t_%05d_%03d.png', curSeqName, round(curT*1E-4)*10, k));
    outFile = fullfile(baseImDir, sprintf('%s%03d.png', curSeqName, k));
    if exist(outFile, 'file'), 
        fprintf('Skipping %d/%d\n', k, nFiles);
        continue; 
    end
    fprintf('Processing %d/%d\n', k, nFiles);
    
    subplot(2,2,1);
    imshow(integralImagePlot);
    title(sprintf('Events t=%5dms (\\Delta t=%3dms)', round(curT*1E-4)*10, tFiltBound*1E-3));
    
    subplot(2,2,2);
    imshow(colorWheelCodedPlot); hold on;
    quiver(1:subSamp:maxX, 1:subSamp:maxY, ...
        respPlot90(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
        respPlot0(1:subSamp:maxY,1:subSamp:maxX)*quiverScale, ...
        'Autoscale', 'off');
    axis ij equal; xlim([1, maxX]); ylim([1, maxY]);
    title(sprintf('Motion (sigma=%4.1f)', sigmaPlotBlurr)); hold off;
    
    
    subplot(2, 2, 3);
    imshow(respPlot0, respPRange);
    title('Vertical motion response');
    
    subplot(2, 2, 4);
    imshow(respPlot90, respPRange);
    title('Horizontal motion response');
    
    drawnow;
%     saveas(gcf, outFile);
    print('-dpng', '-r300', outFile);
    
end
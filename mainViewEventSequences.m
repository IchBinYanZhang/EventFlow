close all;
clear;
clc;

baseDir = '/Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/eventBasedWeizman/weizmannEventCropped';
sequences = dir(fullfile(baseDir, '*.mat'));


%%
for seq_k = 1:numel(sequences)
    curSeqName = sequences(seq_k).name;
    fprintf('%s\n', curSeqName);
    S = load(fullfile(baseDir, curSeqName));
    x_coords = S.x_coords; y_coords = S.y_coords;
    ts = S.ts; on_offs = S.on_offs;
    
%     dt = 0.05E6; % mu s
    dt = 0.3E6; % mu s
    
    figure(10);
    for curT = 1:dt:ts(end)-dt-1
        im = zeros(128, 128, 3);
        fromPlot = find(ts > curT, 1, 'first');
        toPlot = find(ts > curT + dt, 1, 'first');
        for k = fromPlot:toPlot
            if on_offs(k)
                im(y_coords(k), x_coords(k), 1) = im(y_coords(k), x_coords(k), 1) + (toPlot-k)/(toPlot-fromPlot);
            else
                im(y_coords(k), x_coords(k), 3) = im(y_coords(k), x_coords(k), 3) + (toPlot-k)/(toPlot-fromPlot);
            end
        end
        imshow(im, 'InitialMagnification', 500); 
        ylim([35, 119]);
        title(sprintf('%s t=%4.2fs', curSeqName, curT*1E-6), 'Interpreter', 'none');
        drawnow;
    end
    
    onInds = logical(on_offs == 1);
    offInds = logical(on_offs == 0);
    
    figure(20);
    plot3(x_coords(onInds), y_coords(onInds), ts(onInds), '.r', ...
        x_coords(offInds), y_coords(offInds), ts(offInds), '.b')
    xlabel('x'); ylabel('y'); zlabel('t');
    pause;
end






%% 3D plot of events
% borderMask = logical(y_coords>29) & logical(x_coords>5);
% 
% onInds = logical(on_offs == 1) & borderMask;
% offInds = logical(on_offs == 0) & borderMask;
% 
% figure(20);
% plot3(x_coords(onInds), y_coords(onInds), allTs(onInds), '.r', ...
%     x_coords(offInds), y_coords(offInds), allTs(offInds), '.b')
% xlabel('x'); ylabel('y'); zlabel('t');
% 
% 
% zlim([.1E6, 3.5E6]);
% 
% %% 2D plot
% 
% im = zeros(128, 128, 3);
% fromPlot = 55000;
% toPlot = fromPlot + 3000;
% for k = fromPlot:toPlot
%     if on_offs(k)
%         im(y_coords(k), x_coords(k), 1) = im(y_coords(k), x_coords(k), 1) + (toPlot-k)/(toPlot-fromPlot);
%     else
%         im(y_coords(k), x_coords(k), 3) = im(y_coords(k), x_coords(k), 3) + (toPlot-k)/(toPlot-fromPlot);
%     end
% end
% imshow(im, 'InitialMagnification', 500); 
% ylim([30, 128]);
% 
% 
% 

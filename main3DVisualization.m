close all;
clear;
clc;

restoredefaultpath;
addpath('flowHelper');

baseOutDir = 'output';
baseDir = '/Users/bloodbreaker/Desktop/experiment2B/subject4_search_right/';
personNum = 1;
% curSeqName = sprintf('person%02d_bend', personNum);
% curSeqName = sprintf('person%02d_jack', personNum);
% curSeqName = sprintf('person%02d_jump', personNum);
% curSeqName = sprintf('person%02d_pjump', personNum);
% curSeqName = sprintf('person%02d_run', personNum);
% curSeqName = sprintf('person%02d_side', personNum);
% curSeqName = sprintf('person%02d_skip', personNum);
% curSeqName = sprintf('person%02d_walk', personNum);
% curSeqName = sprintf('person%02d_wave01', personNum);
% curSeqName = sprintf('person%02d_wave02', personNum);
curSeqName = 'DVS128-2016-06-16T17-45-16+0200-0220-0'



S = load(fullfile(baseDir, [curSeqName, '.mat']));
x_coords = S.values{1}; y_coords = S.values{2};
ts = S.values{4}; on_offs = S.values{3}; 
trigger = S.values{5};
clear S;

% % add noise to the event time stamps
% [ts, inds] = sort(ts + randn(size(ts))*20E3); % 40E3, i.e. 40ms
% y_coords = y_coords(inds); x_coords = x_coords(inds); on_offs = on_offs(inds);

x_coords = x_coords - min(x_coords(:)) + 1;
y_coords = y_coords - min(y_coords(:)) + 1;
maxY = max(y_coords(:));
maxX = max(x_coords(:));
idx_trigger = find(trigger==1);
%% Play sequence and 3D plot
% %     dtPlot = 0.05E6; % mu s
dtPlot = 100E3; % mu s

% red - on events
% off - off events
% purple - on and off events both happened during that time interval

%% plot moving events 
fig_frame = figure(10);
fig_events = figure(11);
% video = VideoFileReader('/Users/bloodbreaker/Desktop/experiment2B/rec_2016_6_16_17_45.mp4');
for tt = 1:length(idx_trigger)-1
    curT = ts(idx_trigger(tt));
    im = zeros(128, 128, 3);
    fromPlot = idx_trigger(tt);
    toPlot = idx_trigger(tt+1);
    
    for k = fromPlot:toPlot
        if on_offs(k)
            im(129-y_coords(k), 129-x_coords(k), 1) = im(129-y_coords(k), 129-x_coords(k), 1) + (toPlot-k)/(toPlot-fromPlot);
        else
            im(129-y_coords(k), 129-x_coords(k), 3) = im(129-y_coords(k), 129-x_coords(k), 3) + (toPlot-k)/(toPlot-fromPlot);
        end
    end
%     imshow((im), 'InitialMagnification', 500);
    set(fig_events, 'cdata',im);
%     ylim([35, 119]);
    title(sprintf('%s t=%4.2fs', curSeqName, curT*1E-6), 'Interpreter', 'none');
    drawnow; %pause;
end


%% plot the 3D cube
yMask = y_coords > 30;

onInds = logical(on_offs == 1) & yMask;
offInds = logical(on_offs == 0) & yMask;

% specify the plot range
tsRange = ts(end)-ts(1);
lower = 0.75;
upper = 0.85;
mid = (lower + upper)/2;
% plotRange = [ts(1) + lower*tsRange, ts(1) + mid*tsRange, ts(1) + upper*tsRange];
plotRange = [ts(72898), ts(83842), ts(96875)] % frame idx = 500, 550, 600
imm = cell(length(plotRange),1);
for tt = 1:numel(plotRange)
    curT = plotRange(tt);
    imm{tt} = zeros(128, 128, 3);
    fromPlot = find(ts > plotRange(tt)-dtPlot, 1,'first');
    toPlot = find(ts > plotRange(tt)+dtPlot, 1,'first');
    
    for k = fromPlot:toPlot
        if on_offs(k)
            imm{tt}(129-y_coords(k), x_coords(k), 1) = imm{tt}(129-y_coords(k), x_coords(k), 1) + (toPlot-k)/(toPlot-fromPlot);
        else
            imm{tt}(129-y_coords(k), x_coords(k), 3) = imm{tt}(129-y_coords(k), x_coords(k), 3) + (toPlot-k)/(toPlot-fromPlot);
        end
    end
%     im{tt} = uint8(im{tt});
%     imshow((im), 'InitialMagnification', 500);
%     ylim([35, 119]);
    title(sprintf('%s t=%4.2fs', curSeqName, curT*1E-6), 'Interpreter', 'none');
    drawnow; %pause;
end




coloredPlotRange = [plotRange(1), plotRange(1)+2*dtPlot, ...
                    plotRange(2)-dtPlot, plotRange(2)+dtPlot,...
                    plotRange(3)-2*dtPlot, plotRange(3)];
grayPlotRange = [plotRange(1)+dtPlot, plotRange(2)-dtPlot,...
                 plotRange(2)+dtPlot, plotRange(3)-dtPlot];

% cPRfromInd = find(ts>coloredPlotRange(1), 1, 'first');
% cPRtoInd = find(ts>coloredPlotRange(2), 1, 'first');
% gPRfromInd = find(ts>grayPlotRange(1), 1, 'first');
% gPRtoInd = find(ts>grayPlotRange(2), 1, 'first');

cPRIndex = [find((ts > coloredPlotRange(1) & ts < coloredPlotRange(2)));...
            find((ts > coloredPlotRange(3) & ts < coloredPlotRange(4)));...
            find((ts > coloredPlotRange(5) & ts < coloredPlotRange(6)))];...
            
gPRIdx = [find(ts>grayPlotRange(1) & ts<grayPlotRange(2));...
          find(ts>grayPlotRange(3) & ts<grayPlotRange(4))];
    

x_coordsCol = x_coords(cPRIndex); x_coordsGray = x_coords(gPRIdx);
y_coordsCol = y_coords(cPRIndex); y_coordsGray = y_coords(gPRIdx);
ts_Col = ts(cPRIndex); ts_Gray = ts(gPRIdx); yMaskGray = yMask(gPRIdx);
trigger_Gray = trigger(gPRIdx);
onIndsCol = onInds(cPRIndex);
offIndsCol = offInds(cPRIndex);


figure(20);
% myAxes = plot3(x_coordsCol(onIndsCol), ts_Col(onIndsCol), y_coordsCol(onIndsCol),  '.r', ...
%     x_coordsCol(offIndsCol), ts_Col(offIndsCol), y_coordsCol(offIndsCol), '.b', ...
    myAxes = plot3(x_coordsGray(yMaskGray), ts_Gray(yMaskGray), y_coordsGray(yMaskGray), '.k');
xlabel('x'); ylabel('time'); zlabel('y');
set(myAxes, 'Color', 0.75*[1 1 1]);
% myAxes(3).color(4) = 0.5;
% alpha(myAxes(3),0.2);
% zlim([0, 128]);

for tt = 1:length(plotRange)
    hold on;
    surface('XData',[1 128; 1 128],'ZData',[1 1; 128 128],...
            'YData',[plotRange(tt), plotRange(tt); plotRange(tt) plotRange(tt)],'CData',flipdim((imm{tt}),1),...
        'FaceColor','texturemap','EdgeColor','none');
end
% colormap(gray);
















clear;
close all;
clc;

addpath('loadData');
basePath = fullfile('/Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/eventBasedWeizman');
fileName = '2015_06_08WeizmannEvents02.aedat';
filePath = fullfile(basePath, fileName);


%% load and convert data
maxNumEvents = 100E6;
tic; [allAddr,allTs] = loadaerdat(filePath, maxNumEvents);
fprintf('%s: Loaded file "%s" in %6.3fs\n', ...
    datestr(now, 'yyyy_mm_dd_HH_MM_SS'), fileName, toc);

tic; [x_coords, y_coords, allTs, on_offs] = getDataTobi(allTs, allAddr);
fprintf('%s: Converted location information in: %6.3fs\n', ...
    datestr(now, 'yyyy_mm_dd_HH_MM_SS'), toc);

% save('weizmanEB_preproc.mat', ...
%     'x_coords', 'y_coords', 'allTs', 'on_offs');



%% Option 1: with cumsum
eventsPerTime = zeros(allTs(end), 1);
eventsPerTime(allTs) = 1;

tic; evPTcumsum = cumsum(eventsPerTime);
fprintf('%s: Cumsum in: %6.3fs\n', datestr(now, 'yyyy_mm_dd_HH_MM_SS'), toc);

%%
pS = 100; % plot sample
plot(1:pS:numel(evPTcumsum), evPTcumsum(1:pS:end), '-k');
drawnow;

thres = 40000; 
deltaT = 200000; % mu s
stepFw = 1.4E6; % jump forward 1s when a marker is found
episodeMarkers = 1;
k = 1;
testUntil = numel(evPTcumsum)-deltaT-1;
while k < testUntil
    if evPTcumsum(k+deltaT) - evPTcumsum(k) > thres, 
        episodeMarkers(end+1) = k; 
        k = k + stepFw;
        fprintf('Found a fw marker at %d\n', episodeMarkers(end));
        hold on; plot(ones(2, 1)*episodeMarkers(end), [0, 3.4e7], '-g'); drawnow;
    end
    k = k+1;
end


k = numel(evPTcumsum);
while k > deltaT
    if evPTcumsum(k) - evPTcumsum(k-deltaT) > thres, 
        episodeMarkers(end+1) = k; 
        k = k - stepFw;
        fprintf('Found a bw marker at %d\n', episodeMarkers(end));
        hold on; plot(ones(2, 1)*episodeMarkers(end), [0, 3.4e7], '-b'); drawnow;
    end
    k = k-1;
end
hold off;


save('weizmanEB02_withEpisodeMarkers.mat', ...
    'episodeMarkers', 'x_coords', 'y_coords', 'allTs', 'on_offs', 'thres', 'deltaT', 'stepFw');


%% Episode markers are shifted by 130000 mu s
episodeMarkers = sort(episodeMarkers);
plot(1:pS:numel(evPTcumsum), evPTcumsum(1:pS:end), '-k');
finalEpisodeMarkers = zeros(2, 93); count = 1;
hold on;
for k = 1:2:numel(episodeMarkers)-1
    from = episodeMarkers(k) - 130000; 
    if k == 1, from = episodeMarkers(k) + 400000; end
    if k == 103, from = from + 900000; end
    to = episodeMarkers(k+1)+130000; 
    if k == 97, to = to-1400000; end
    if k == 99, to = to-200000; end
    if k == 103, to = to-1600000; end
    plot(ones(2, 1)*from, [0,  3.4e7], '-g', ...
        ones(2, 1)*to, [0,  3.4e7], '-b');
    text(from, evPTcumsum(from), sprintf('%d', k));
    finalEpisodeMarkers(:, count) = [from, to]; count = count + 1;
end
hold off;


%% 
str = '/Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_bend.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_jack.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_jump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_pjump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_run.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_side.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_skip.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_walk.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_wave01.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person01_wave02.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_bend.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_jack.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_jump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_pjump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_run.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_side.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_skip.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_walk.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_wave01.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person02_wave02.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_bend.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_jack.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_jump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_pjump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_run.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_side.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_skip.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_walk.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_wave01.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person03_wave02.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_bend.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_jack.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_jump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_pjump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_run.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_side.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_skip.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_walk.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_wave01.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person04_wave02.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_bend.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_jack.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_jump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_pjump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_run.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_side.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_skip.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_walk.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_wave01.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person05_wave02.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_bend.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_jack.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_jump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_pjump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_run1.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_run2.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_side.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_skip1.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_skip2.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_walk1.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_walk2.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_wave01.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person06_wave02.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_bend.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_jack.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_jump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_pjump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_run.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_side.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_skip.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_walk.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_wave01.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person07_wave02.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_bend.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_jack.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_jump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_pjump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_run.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_side.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_skip.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_walk.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_wave01.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person08_wave02.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_bend.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_jack.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_jump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_pjump.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_run.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_side.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_skip.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_walk.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_wave01.avi /Users/tbhunderbird/Daten/TestSets/Weizmann_HumanActions/movie/person09_wave02.avi';
episodeNames = textscan(str, '%s'); episodeNames = episodeNames{1};
for k = 1:numel(episodeNames),
    parts = strsplit(episodeNames{k}, '/');
    name = strsplit(parts{end}, '.'); name = name{1};
    
    fromTo = finalEpisodeMarkers(:, k);
    fromInd = find(allTs > fromTo(1), 1, 'first');
    toInd = find(allTs > fromTo(2), 1, 'first');
    
    fprintf('%2d: %20s (%8d, %8d)\n', k, name, fromInd, toInd);
    
    xCoords = x_coords(fromInd:toInd);
    yCoords = y_coords(fromInd:toInd);
    ts = allTs(fromInd:toInd)-allTs(fromInd)+1;
    onOffs = on_offs(fromInd:toInd);
    
    save(fullfile('weizmannEvent', [name, '.mat']), 'xCoords', 'yCoords', 'ts', 'onOffs');
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

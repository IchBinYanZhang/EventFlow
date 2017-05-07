function [x_coords, y_coords, allTsnew, on_offs] = getDataTobi(allTs, allAddr)
%GETDATATOBI 
% make sure that timestamps are really monotonically increasing
[allTs, inds] = sort(allTs);
allAddr = allAddr(inds);

% preallocate
numEvents=numel(allTs);

x_coords = zeros(1, numEvents);
y_coords = zeros(1, numEvents);
allTsnew = zeros(1, numEvents);
on_offs = zeros(1, numEvents);

% convert address to readable coordinates
count=0;
for currentEvent=1:numEvents(1)
    curUint = allAddr(currentEvent);
    onOffBit = mod(curUint, 2); % last byte
%     curAddr = double(bitshift(curUint, -1, 'uint32')); % first bytes except last
    curAddr = double(bitshift(curUint, -1)); % first bytes except last

    column = mod(curAddr, 128);
    row = (curAddr - column) / 128;
    column = 128-column; row = 128-row;
    
    if (column>0) && (column<129) && (row>0) && (row<129),
        count = count + 1;
        x_coords(count) = column;
        y_coords(count) = row;
        allTsnew(count) = allTs(currentEvent);        
        on_offs(count) = ~onOffBit;
    end
end
if numEvents-count>0, 
    fprintf('WARNING: Sorted %d/%d events out\n', numEvents-count, numEvents); 
end

% strip to valid length count
x_coords = x_coords(1:count);
y_coords = y_coords(1:count);
allTsnew = allTsnew(1:count) - allTsnew(1)+1;
on_offs = on_offs(1:count);

% % only on or offs
% validInds = (on_offs == 0);
% x_coords = x_coords(validInds);
% y_coords = y_coords(validInds);
% allTsnew = allTsnew(validInds);
% on_offs = on_offs(validInds);

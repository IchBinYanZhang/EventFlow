function [x_coord, y_coord, allTsnew, on_off] = getDataStephan(allTs, allAddr)
%GETDATASTEPHAN

ON='1'; % OF='0';
total_t=numel(allTs);
count=1;

x_coord = zeros(1, total_t);
y_coord = zeros(1, total_t);
allTsnew = zeros(1, total_t);
on_off = zeros(1, total_t);

for currentEvent=1:total_t(1)
    
    %% Stephan
    % Extract address from allAddr format(thanks to Stephan Tschechne)
    bs1 = dec2bin( allAddr(currentEvent), 15 );
    d=bs1(15);
    currentAddress = bs1(1:14);
    currentAddressBinary = bin2dec(currentAddress);
    column = mod(currentAddressBinary ,128);
    row = (currentAddressBinary-column) / 128;
    dx = 128-column;
    dy = 128-row;
    
    if (dx>0) && (dx<129) && (dy>0) && (dy<129),
        x_coord(count) = dx;
        y_coord(count) = dy;
        allTsnew(count)= allTs(currentEvent);
        
        % Extract ON and OFF events
        if  strcmp(d,ON)    %ON
            on_off(count)=1;
        else                %OFF
            on_off(count)=0;
        end
        count=count+1;
    end
    
end

% strip to valid length count
x_coord = x_coord(1:count-1);
y_coord = y_coord(1:count-1);
allTsnew = allTsnew(1:count-1) - allTsnew(1);
on_off = on_off(1:count-1);




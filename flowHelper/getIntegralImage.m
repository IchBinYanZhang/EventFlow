function [ intIm ] = getIntegralImage(y_coords, x_coords, on_offs, maxY, maxX)
%GETINTEGRALIMAGE
    intIm = zeros(maxY, maxX, 3);
    for k = 1:numel(y_coords)
        if on_offs(k) == 1, 
            intIm(y_coords(k), x_coords(k), 1) = intIm(y_coords(k), x_coords(k), 1) + 1;
        else
            intIm(y_coords(k), x_coords(k), 3) = intIm(y_coords(k), x_coords(k), 3) + 1;
        end
    end
end


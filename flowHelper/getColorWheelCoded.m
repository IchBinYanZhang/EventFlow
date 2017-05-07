function anglesPlot = getColorWheelCoded(vecY, vecX, maxSpeed)
    % get color coded output of vector map vecY, vecX
    % using quiver with axis ij (i.e. in imshow direction),
    % vector (1, 1) points into lower right corner
    %  \
    %   \
    %    \ |
    %    __|
    if maxSpeed <= 0, error('maxSpeed must be positive'); end
    angles = atan2(-vecY, vecX); % "-" because of plot dir of imshow
    velocities = hypot(vecY, vecX);
    anglesPlot = ind2rgb(ceil((angles+pi)/(2*pi)*256), hsv(256));
    scaling = min(1, velocities * 1/maxSpeed);
%     anglesPlot = anglesPlot .* repmat(scaling, [1, 1, 3]); % slow speeds are black
    anglesPlot = anglesPlot .* repmat(scaling, [1, 1, 3]) + repmat(1-scaling, [1, 1, 3]); % slow speeds are white
end
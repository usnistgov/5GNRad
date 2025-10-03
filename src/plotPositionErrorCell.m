function [] = plotPositionErrorCell(x,y,errorValues,x1,y1, varargin)
% PLOTPOSITIONERRORCELL Visualize position error on a 2D spatial map
%   PLOTPOSITIONERRORCELL(X, Y, ERRORVALUES, X1, Y1) plots a 2D scatter
%   plot of position errors given by ERRORVALUES at coordinates (X,Y).
%   Points (X1,Y1) are highlighted in red to indicate miss-detections or
%   other events of interest.
%
%   PLOTPOSITIONERRORCELL(..., 'cellSize', [XMIN XMAX YMIN YMAX]) allows
%   specification of plot bounds.
%
%   Parameters:
%     X, Y            - Coordinates of the error data points
%     ERRORVALUES     - Scalar error values to visualize (e.g., localization error in meters)
%     X1, Y1          - Coordinates of points to highlight (e.g., miss-detections)
%
%   Optional Name-Value Pair:
%     'cellSize'      - [XMIN XMAX YMIN YMAX], bounding box for the plot area
%
%   The function supports two display modes:
%     * Interpolated contour map (if plotInterp = 1)
%     * Scatter plot with color-coded errors (default)

%   2025 NIST/CTL Steve Blandino

%   This file is available under the terms of the NIST License.

% Define input parser
p = inputParser;

% Default values for cellSize (min(x), max(x), min(y), max(y))
defaultCellSize = [min(x), max(x), min(y), max(y)];

% Add parameter for cellSize
addParameter(p, 'cellSize', defaultCellSize, @(v) isnumeric(v) && numel(v) == 4);

% Parse inputs
parse(p, varargin{:});
cellSize = p.Results.cellSize;

plotInterp = 0;

if plotInterp
    % Define grid for interpolation based on cellSize
    xq = linspace(cellSize(1), cellSize(2), 200);  % X grid points
    yq = linspace(cellSize(3), cellSize(4), 200);  % Y grid points

    % Define grid for interpolation
    % xq = linspace(min(x), max(x), 50);  % X grid points
    % yq = linspace(min(y), max(y), 50);  % Y grid points
    [Xq, Yq] = meshgrid(xq, yq);

    % Interpolate scattered data onto a grid
    ErrorGrid = griddata(x, y, errorValues, Xq, Yq, 'cubic'); % 'cubic' for smooth interpolation

    % Plot the contour map
    figure;
    contourf(Xq, Yq, ErrorGrid, 3); % 20 levels, no line edges
    c = colorbar;
    c.Label.String = 'Position Error (m)';
    xlabel('X (m)');
    ylabel('Y (m)');
    hold on;
    % scatter(x, y, 30, errorValues, 'filled', 'MarkerEdgeColor', 'k'); % Show data points
    % hold off;
    % plot(x, y, 'ro', 'MarkerFaceColor','r')
    % title('2D Contour Plot of Error');
    p = plot(x1, y1, 'ro', 'MarkerFaceColor','r');
    legend(p, 'Miss-Detection')

else
    figure
    scatter(x, y, 40, errorValues, 'filled', 'MarkerEdgeColor', 'k','MarkerEdgeAlpha',.2)
    c = colorbar;
    c.Label.String = 'Position Error (m)';
    xlabel('X (m)');
    ylabel('Y (m)');
    hold on;
    p = plot(x1, y1, 'ro', 'MarkerFaceColor','r');
    legend(p, 'Miss-Detection')
end
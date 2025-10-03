%% Analyze 5GNRad Simulation Results
% This helper script analyzes and visualizes the output of 5GNRad
% simulations across one or more scenarios. It performs the following:
%
%   - Loads simulation output from /Output/error.csv for each scenario
%   - Computes and visualizes 2D spatial error maps
%   - Identifies miss-detections based on a position error threshold 
%   - Plots miss-detection rates and cumulative distribution of position errors
%
% Example:
%   Run after simulation completes and error.csv is generated in:
%       examples/UMi-Av25/Output/error.csv

%   2025 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.


clear
close all
addpath(genpath(pwd))

folder = 'examples3GPP';
scenarioFolderNames = dir(folder);
scenarioFolderNames(1:2) = [];
scenarioNameStrVector = {scenarioFolderNames.name};
% scenarioNameStrVector = { 'UMa-Av200-8x8-60'};

nscen = length(scenarioNameStrVector);
for i = 1:nscen
    scenarioPath = fullfile(pwd, folder,scenarioNameStrVector{i});
    errorFile = fullfile(scenarioPath, 'Output/error.csv');
    tgtFile = fullfile(scenarioPath, 'Input/targetConfig.txt');
    tgt = readmatrix(tgtFile);
    tgt = tgt(:,1:3);
    tgtMat(:,:, i) = tgt;
    T = readtable(errorFile);
    nz = T.positionError~=0;
    T = table2array(T);
    errorMat(:,:,i) = T(nz,:);

end

cellsize= 500;
eta = 3.4;

%% Cell plot
for scen = 1:nscen
    trueMissedDetection = errorMat(:,7,scen) < eta;
    errorValues = squeeze(errorMat(~trueMissedDetection,1, scen));

    plotPositionErrorCell(tgtMat(~trueMissedDetection,1,scen),tgtMat(~trueMissedDetection,2,scen), errorValues, tgtMat(trueMissedDetection,1,scen),tgtMat(trueMissedDetection,2,scen), 'cellSize', (mod(0:3,2)*2-1)*cellsize)
    clim([0 10])
    
    axis equal


    % set(gca, "YTick", -cellsize(mod(ceil(scen/nscen), nscen))/2:cellsize(mod(ceil(scen/nscen), nscen))/nscen:cellsize(mod(ceil(scen/nscen), nscen))/2)
    % set(gca, "XTick", -cellsize(mod(ceil(scen/nscen), nscen))/2:cellsize(mod(ceil(scen/nscen), nscen))/nscen:cellsize(mod(ceil(scen/nscen), nscen))/2)
    % axis([-cellsize(mod(ceil(scen/nscen), nscen))/(2*sind(60)) cellsize(mod(ceil(scen/nscen), nscen))/(2*sind(60)) -cellsize(mod(ceil(scen/nscen), nscen))/2 cellsize(mod(ceil(scen/nscen), nscen))/2])


    % set(gca, 'Position', [ 0.2451    0.2932    0.4962    0.4962])

end


%% Other
for scen = 1:nscen
    mdMat(scen) = sum(errorMat(:,7,scen) < eta);
end

% mdMat = reshape(mdMat/size(tgtMat,1)*100, 4,[])';


% Define x-axis labels
xLabels = {'UMi-AV'};
xLabels = scenarioNameStrVector;
% Create bar plot
figure;
bar(mdMat);

% Set x-axis labels
set(gca, 'XTickLabel', xLabels);

% Add labels and title
xlabel('Scenario');
ylabel('Miss Detection Probability');

% Add legend

% Set font size for readability
set(gca, 'FontSize', 12);

% Display the plot
set(gca, 'YGrid', 'on', 'XGrid', 'off')
ylim([0 25])


%%
figure
for scen = 1:nscen
    % if mod(scen,4)==1
    %     figure
    % end
    trueMissedDetection =  errorMat(:,7,scen) < eta;
    errorValues = squeeze(errorMat(~trueMissedDetection,1, scen));
    p = cdfplot(errorValues);
    p.LineWidth = 4;
    hold on
end

l.FontSize = 18;
title('')
% set(gca, "XTick", 0:2:10)
xlabel('Position Error (m)')

%%
figure
for scen = 1:nscen
    % if mod(scen,4)==1
    %     figure
    % end
    trueMissedDetection =  errorMat(:,7,scen) < eta;
    errorValues = squeeze(errorMat(~trueMissedDetection,3, scen));
    p = cdfplot(errorValues);
    p.LineWidth = 4;
    hold on
end

l.FontSize = 18;
title('')
% set(gca, "XTick", 0:2:10)
xlabel('Velocity Error (m)')


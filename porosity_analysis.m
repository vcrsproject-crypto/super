%% porosity_analysis.m
% =====================================================================
% POROSITY ANALYSIS FROM A MICROGRAPH IMAGE  (single-file version)
% ---------------------------------------------------------------------
% Reads an image, REMOVES THE SCALE-BAR REGION first (so the burned-in
% scale bar / text is not mistaken for a pore), segments the dark pores
% from the lighter matrix, cleans up noise, and reports porosity (%).
%
%       Porosity (%) = (pore pixels / valid pixels) * 100
%
% HOW TO USE
%   1) Put your micrograph in the same folder as this file.
%   2) Set imageFile below to its name.
%   3) Press Run.
%
% Requires the Image Processing Toolbox.
% =====================================================================

clear; clc; close all;

% ----- 1. Single image -------------------------------------------------
imageFile = '1 100x a.jpg';      % <-- change to your image file name

[porosity, results] = porosityAnalysis(imageFile, ...
    'ScaleBarCorner', 'bottomright', ...  % where the scale bar sits
    'ScaleBarFrac',   [0.18 0.10],  ...   % size of the corner to remove
    'Threshold',      [],           ...   % [] = automatic (Otsu)
    'Invert',         false,        ...   % pores are dark
    'MinPoreArea',    15,           ...   % remove speckle smaller than this
    'Show',           true);

fprintf('\n>> Porosity of %s = %.2f %%\n', imageFile, porosity);


% ----- 2. (Optional) Batch process a whole folder ----------------------
% Uncomment to process every JPG/PNG/TIF in the current folder and write a
% summary table to porosity_results.csv
%
% files = [dir('*.jpg'); dir('*.png'); dir('*.tif'); dir('*.tiff')];
% names = strings(numel(files),1);
% vals  = zeros(numel(files),1);
% for k = 1:numel(files)
%     names(k) = string(files(k).name);
%     vals(k)  = porosityAnalysis(files(k).name, 'Show', false);
% end
% T = table(names, vals, 'VariableNames', {'Image','Porosity_pct'});
% writetable(T, 'porosity_results.csv');
% disp(T);


% =====================================================================
%                          LOCAL FUNCTIONS
% =====================================================================

function [porosityPct, results] = porosityAnalysis(imagePath, varargin)
%POROSITYANALYSIS  Measure area porosity (%) from a micrograph image.
%
%   p = porosityAnalysis('sample.jpg');
%   [p, results] = porosityAnalysis('sample.jpg', 'Name', Value, ...);
%
%   NAME-VALUE OPTIONS
%     'ScaleBarCorner' : Where the scale bar sits. One of
%                        'bottomright' (default) | 'bottomleft' |
%                        'topright' | 'topleft' | 'none'
%     'ScaleBarFrac'   : [widthFrac heightFrac] of the corner to remove,
%                        as fractions of the image size. Default [0.18 0.10].
%     'CropStrip'      : If true, physically crop the full bottom/top strip
%                        that contains the scale bar (keeps a clean rectangle).
%                        If false (default), only the corner rectangle is
%                        masked out and excluded from the calculation.
%     'Threshold'      : [] for automatic Otsu threshold (default), or a
%                        scalar in [0,1] to force a manual threshold.
%     'Invert'         : false (default) -> pores are DARK (typical metallo-
%                        graphic / SEM images). true -> pores are BRIGHT.
%     'MinPoreArea'    : Remove connected pore blobs smaller than this many
%                        pixels (speckle noise removal). Default 15.
%     'FillHoles'      : Fill enclosed holes inside detected pores. Default true.
%     'EnhanceContrast': Apply imadjust contrast stretch before thresholding.
%                        Default true.
%     'Show'           : Show the 6-panel result figure. Default true.
%
%   OUTPUTS
%     porosityPct : Scalar porosity percentage.
%     results     : Struct with intermediate data (mask, threshold, counts).

% --- Parse inputs ----------------------------------------------------
p = inputParser;
p.FunctionName = 'porosityAnalysis';
addRequired(p,  'imagePath', @(x) ischar(x) || isstring(x));
addParameter(p, 'ScaleBarCorner', 'bottomright', ...
    @(x) any(strcmpi(x, {'bottomright','bottomleft','topright','topleft','none'})));
addParameter(p, 'ScaleBarFrac', [0.18 0.10], ...
    @(x) isnumeric(x) && numel(x)==2 && all(x>=0 & x<1));
addParameter(p, 'CropStrip', false, @(x) islogical(x) || isscalar(x));
addParameter(p, 'Threshold', [], @(x) isempty(x) || (isscalar(x) && x>=0 && x<=1));
addParameter(p, 'Invert', false, @(x) islogical(x) || isscalar(x));
addParameter(p, 'MinPoreArea', 15, @(x) isscalar(x) && x>=0);
addParameter(p, 'FillHoles', true, @(x) islogical(x) || isscalar(x));
addParameter(p, 'EnhanceContrast', true, @(x) islogical(x) || isscalar(x));
addParameter(p, 'Show', true, @(x) islogical(x) || isscalar(x));
parse(p, imagePath, varargin{:});
opt = p.Results;

% --- 1. Read the image ----------------------------------------------
rgb = imread(char(opt.imagePath));
if ndims(rgb) == 3 && size(rgb,3) == 4
    rgb = rgb(:,:,1:3);                 % drop alpha channel if present
end

% --- 2. Remove / crop the scale-bar region FIRST --------------------
[H, W, ~] = size(rgb);
validMask = true(H, W);                 % TRUE = analyse, FALSE = ignore

if ~strcmpi(opt.ScaleBarCorner, 'none')
    wPix = max(1, round(opt.ScaleBarFrac(1) * W));   % corner width  (pixels)
    hPix = max(1, round(opt.ScaleBarFrac(2) * H));   % corner height (pixels)

    switch lower(opt.ScaleBarCorner)
        case 'bottomright', rows = (H-hPix+1):H; cols = (W-wPix+1):W;
        case 'bottomleft',  rows = (H-hPix+1):H; cols = 1:wPix;
        case 'topright',    rows = 1:hPix;       cols = (W-wPix+1):W;
        case 'topleft',     rows = 1:hPix;       cols = 1:wPix;
    end

    if opt.CropStrip
        % Physically crop the full-width strip and rebuild the mask.
        if any(strcmpi(opt.ScaleBarCorner, {'bottomright','bottomleft'}))
            keepRows = 1:(H - hPix);             % drop bottom strip
        else
            keepRows = (hPix + 1):H;             % drop top strip
        end
        rgb       = rgb(keepRows, :, :);
        [H, W, ~] = size(rgb);
        validMask = true(H, W);
    else
        validMask(rows, cols) = false;          % exclude the corner
    end
end

% --- 3. Grayscale (+ optional contrast enhancement) -----------------
if ndims(rgb) == 3
    gray = rgb2gray(rgb);
else
    gray = rgb;
end
gray = im2double(gray);                 % normalise to [0,1]
if opt.EnhanceContrast
    gray = imadjust(gray);              % stretch contrast
end

% --- 4. Threshold to binary -----------------------------------------
if isempty(opt.Threshold)
    level = graythresh(gray(validMask));    % Otsu on valid region only
else
    level = opt.Threshold;
end
bw = imbinarize(gray, level);           % TRUE = bright matrix, FALSE = dark pore
if opt.Invert
    poreMask = bw;                      % pores are bright
else
    poreMask = ~bw;                     % pores are dark (default)
end

% --- 5. Clean up the pore mask --------------------------------------
poreMask = poreMask & validMask;            % never count the scale-bar area
if opt.FillHoles
    poreMask = imfill(poreMask, 'holes');
end
if opt.MinPoreArea > 0
    poreMask = bwareaopen(poreMask, round(opt.MinPoreArea));
end
poreMask = poreMask & validMask;            % re-apply after cleaning

% --- 6. Compute porosity --------------------------------------------
porePixels  = nnz(poreMask);
validPixels = nnz(validMask);
porosityPct = 100 * porePixels / validPixels;
matrixPct   = 100 - porosityPct;
cc          = bwconncomp(poreMask);
numPores    = cc.NumObjects;

% --- 7. Package results ---------------------------------------------
results = struct( ...
    'porosityPct', porosityPct, ...
    'matrixPct',   matrixPct, ...
    'threshold',   level, ...
    'poreMask',    poreMask, ...
    'validMask',   validMask, ...
    'numPores',    numPores, ...
    'porePixels',  porePixels, ...
    'validPixels', validPixels, ...
    'gray',        gray);

fprintf('--- Porosity Analysis: %s ---\n', char(opt.imagePath));
fprintf('Threshold (Otsu/manual) : %.3f\n', level);
fprintf('Detected pores          : %d\n', numPores);
fprintf('Pore pixels / valid     : %d / %d\n', porePixels, validPixels);
fprintf('Porosity                : %.2f %%\n', porosityPct);
fprintf('Matrix                  : %.2f %%\n', matrixPct);

% --- 8. Visualise (optional) ----------------------------------------
if opt.Show
    showPorosityResults(rgb, gray, bw, poreMask, ...
        porosityPct, matrixPct, level, char(opt.imagePath));
end
end


function showPorosityResults(rgb, gray, bw, poreMask, ...
        porosityPct, matrixPct, level, name)
%SHOWPOROSITYRESULTS  Draw the 6-panel summary figure.

% Red overlay of detected pores on the original image.
if ndims(rgb) == 3
    base = im2double(rgb);
else
    base = repmat(im2double(rgb), 1, 1, 3);
end
R = base(:,:,1); G = base(:,:,2); B = base(:,:,3);
R(poreMask) = 1; G(poreMask) = 0; B(poreMask) = 0;
overlay = cat(3, R, G, B);

figure('Name', sprintf('Porosity Analysis - %s', name), 'Color', 'w');
subplot(2,3,1); imshow(rgb);     title('1. Original Image');
subplot(2,3,2); imshow(gray);    title('2. Grayscale (Enhanced)');
subplot(2,3,3); imshow(bw);      title(sprintf('3. Binary (thresh = %.3f)', level));
subplot(2,3,4); imshow(poreMask);title('4. Cleaned Pore Mask');
subplot(2,3,5); imshow(overlay); title('5. Pores Highlighted (Red)');
subplot(2,3,6);
pie([porosityPct, matrixPct], {sprintf('Porosity\n%.2f%%', porosityPct), ...
                               sprintf('Matrix\n%.2f%%', matrixPct)});
title('6. Composition');
sgtitle(sprintf('Porosity Analysis - %s', name), 'FontWeight', 'bold');
end

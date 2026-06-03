%% porosity_analysis.m
% =====================================================================
% ADVANCED POROSITY & DEFECT ANALYSIS FROM A MICROGRAPH IMAGE
% (single self-contained MATLAB file)
% ---------------------------------------------------------------------
% FEATURES
%   1. Advanced scale-bar detection (black & white, horizontal & vertical)
%      with connected-component analysis, rectangularity / aspect-ratio /
%      geometric filtering, confidence scoring, OCR text reading and
%      automatic micron-per-pixel calibration. The bar + annotation text
%      are removed (inpainted) without disturbing the microstructure.
%   2. Preprocessing: grayscale, illumination correction / background
%      normalization, bilateral filtering, CLAHE, adaptive contrast.
%   3. Segmentation: adaptive threshold, Otsu (optional), Canny edges,
%      hybrid (intensity + edge), automatic sensitivity optimization,
%      small-pore detection.
%   4. Morphological refinement: scratch / dust removal, hole filling,
%      opening / closing, boundary smoothing, watershed separation.
%   5. Per-pore characterization (area, perimeter, equivalent diameter,
%      Feret diameter, major / minor axis, circularity, solidity,
%      aspect ratio, eccentricity, orientation, centroid).
%   6. Automatic defect classification (gas, keyhole, lack-of-fusion,
%      crack, irregular) with user-adjustable thresholds.
%   7. Statistics (porosity %, count, mean / median / std / min / max,
%      95% confidence interval, area-fraction distribution).
%   8. Publication-quality visualization (image montage + statistics).
%
% USAGE
%   1) Put your micrograph in the same folder as this file.
%   2) Set imageFile below.
%   3) Adjust any options in the opt struct (all optional).
%   4) Press Run.
%
% Requires the Image Processing Toolbox. OCR uses the Computer Vision
% Toolbox if available (it degrades gracefully if not).
% =====================================================================

clear; clc; close all;

% ----- INPUT ----------------------------------------------------------
imageFile = '1 100x a.jpg';        % <-- change to your image file name

% ----- OPTIONS (override any default here) ----------------------------
opt = getDefaultOptions();
% Examples of common overrides:
%   opt.SegmentationMethod = 'hybrid';   % 'adaptive' | 'otsu' | 'edge' | 'hybrid'
%   opt.MicronsPerPixel    = 0.25;       % force calibration, skip auto
%   opt.ScaleBarText       = '75 um';    % force scale length if OCR fails
%   opt.SaveCSV            = true;        % export per-pore table to CSV

% ----- RUN ------------------------------------------------------------
results = analyzePorosity(imageFile, opt);


% =====================================================================
%                          LOCAL FUNCTIONS
% =====================================================================

function opt = getDefaultOptions()
%GETDEFAULTOPTIONS  All tunable parameters in one place.

% --- Scale-bar detection / calibration ---
opt.DetectScaleBar     = true;
opt.MicronsPerPixel    = [];         % [] = auto from scale bar; else forced
opt.ScaleBarText       = '';         % e.g. '50 um' to force length if OCR fails
opt.SearchBorderFrac   = 0.35;       % bars usually within this border fraction
opt.MinBarAspect       = 2.5;        % long thin object
opt.MinRectangularity  = 0.55;       % Area / boundingbox area (extent)
opt.ScaleBarPadFrac    = 0.04;       % padding (frac of image) around removed band
opt.UseOCR             = true;

% --- Preprocessing ---
opt.IlluminationCorrect = true;
opt.BilateralFilter     = true;
opt.CLAHE               = true;
opt.CLAHEClipLimit      = 0.01;
opt.AdaptiveContrast    = true;

% --- Segmentation ---
opt.SegmentationMethod  = 'adaptive'; % 'adaptive'|'otsu'|'edge'|'hybrid'
opt.AdaptiveSensitivity = 0.50;       % used if AutoOptimize=false
opt.AutoOptimize        = true;       % auto-tune adaptive sensitivity
opt.PoresAreDark        = true;       % dark pores on bright matrix
opt.CannyThreshold      = [];         % [] = automatic

% --- Morphology ---
opt.MinPoreArea         = 6;          % px; detect very small pores
opt.FillHoles           = true;
opt.OpenRadius          = 1;
opt.CloseRadius         = 1;
opt.SmoothBoundaries    = true;
opt.RemoveScratches     = true;
opt.ScratchMaxThickness = 2.5;        % px; thinner+long => scratch
opt.ScratchMinAspect    = 12;
opt.Watershed           = true;
opt.WatershedDepth      = 1.5;        % imhmin depth (px) – higher = fewer splits

% --- Classification thresholds (user-adjustable) ---
opt.Class.GasCircularity   = 0.80;
opt.Class.GasAspectMax     = 1.60;
opt.Class.GasSolidity      = 0.90;
opt.Class.KeyholeDiam      = 60;      % >= this (um if calibrated, else px) & round => keyhole
opt.Class.LoFAspect        = 2.00;
opt.Class.LoFCircularity   = 0.55;
opt.Class.CrackAspect      = 4.00;
opt.Class.CrackCircularity = 0.30;
opt.Class.CrackSolidity    = 0.60;

% --- Output ---
opt.ShowFigures = true;
opt.SaveCSV     = false;
opt.CSVName     = 'pore_results.csv';
end


function results = analyzePorosity(imagePath, opt)
%ANALYZEPOROSITY  Full pipeline. Returns a results struct.

muStr = char(181);   % micro sign for display

% ---------------- 1. READ IMAGE ----------------
rgb = imread(char(imagePath));
if ndims(rgb) == 3 && size(rgb,3) == 4, rgb = rgb(:,:,1:3); end
if ndims(rgb) == 3, gray0 = rgb2gray(rgb); else, gray0 = rgb; end
gray0 = im2double(gray0);
[H, W] = size(gray0);

% ---------------- 2. SCALE-BAR DETECTION & CALIBRATION ----------------
sb = detectScaleBar(gray0, opt);

% Establish calibration
if ~isempty(opt.MicronsPerPixel)
    mpp = opt.MicronsPerPixel; calibrated = true; calSrc = 'user';
elseif sb.found && ~isnan(sb.micronsPerPixel)
    mpp = sb.micronsPerPixel; calibrated = true; calSrc = 'scalebar';
else
    mpp = 1; calibrated = false; calSrc = 'none (pixels)';
end
unit     = ternary(calibrated, 'um', 'px');          % ASCII token (column names)
unitDisp = ternary(calibrated, [muStr 'm'], 'px');   % pretty label for display
aunit    = ternary(calibrated, [muStr 'm^2'], 'px^2');

% Build valid (analysis) mask excluding scale bar + annotation, and an
% inpainted image so the removed region does not create false edges.
validMask = true(H, W);
grayClean = gray0;
if sb.found
    validMask(sb.removalMask) = false;
    try
        grayClean = regionfill(gray0, sb.removalMask);   % inpaint band
    catch
        grayClean(sb.removalMask) = median(gray0(~sb.removalMask));
    end
end

% ---------------- 3. PREPROCESSING ----------------
enh = grayClean;
if opt.IlluminationCorrect
    sigma = max(15, round(0.15*min(H,W)));
    bg  = imgaussfilt(enh, sigma);
    enh = enh ./ (bg + eps);
    enh = mat2gray(enh);
end
if opt.BilateralFilter
    try,  enh = imbilatfilt(enh);
    catch, enh = imgaussfilt(enh, 1); end
end
if opt.CLAHE
    try,  enh = adapthisteq(enh, 'ClipLimit', opt.CLAHEClipLimit, ...
                            'Distribution', 'rayleigh');
    catch, enh = histeq(enh); end
end
if opt.AdaptiveContrast
    enh = imadjust(enh, stretchlim(enh, [0.005 0.995]));
end

% Image whose BRIGHT pixels correspond to pores (simplifies thresholding)
if opt.PoresAreDark, poreImg = imcomplement(enh); else, poreImg = enh; end

% ---------------- 4. SEGMENTATION ----------------
[rawMask, segInfo] = segmentPores(poreImg, enh, validMask, opt);

% ---------------- 5. MORPHOLOGICAL REFINEMENT ----------------
mask = refineMask(rawMask, validMask, opt);

% ---------------- 6./7. CHARACTERIZATION + CLASSIFICATION ----------------
[poreTable, L] = characterizePores(mask, mpp, unit, opt);
stats = computeStatistics(poreTable, mask, validMask, unit, unitDisp, aunit);

% ---------------- ASSEMBLE RESULTS ----------------
results = struct();
results.imagePath     = char(imagePath);
results.rgb           = rgb;
results.gray          = gray0;
results.enhanced      = enh;
results.scaleBar      = sb;
results.calibrated    = calibrated;
results.calSource     = calSrc;
results.micronsPerPixel = mpp;
results.unit          = unit;
results.validMask     = validMask;
results.binaryMask    = mask;
results.labelMatrix   = L;
results.segInfo       = segInfo;
results.poreTable     = poreTable;
results.stats         = stats;

printSummary(results);

if opt.SaveCSV && ~isempty(poreTable)
    writetable(poreTable, opt.CSVName);
    fprintf('Per-pore table written to %s\n', opt.CSVName);
end
if opt.ShowFigures
    visualizeResults(results, opt);
end
end


% ---------------------------------------------------------------------
% SCALE-BAR DETECTION
% ---------------------------------------------------------------------
function sb = detectScaleBar(gray, opt)
%DETECTSCALEBAR  Detect a black or white, horizontal or vertical scale bar.

[H, W] = size(gray);
sb = struct('found', false, 'bbox', [], 'lengthPix', NaN, ...
            'orientation', '', 'polarity', '', 'confidence', 0, ...
            'text', '', 'physicalLength', NaN, 'physicalUnit', '', ...
            'micronsPerPixel', NaN, 'removalMask', false(H,W));

if ~opt.DetectScaleBar, return; end

% Candidate binary maps for dark and bright bars.
candidates = struct('mask', {}, 'polarity', {});
candidates(end+1) = struct('mask', gray < 0.18, 'polarity', 'black');
candidates(end+1) = struct('mask', gray > 0.90, 'polarity', 'white');

borderF = opt.SearchBorderFrac;
best = sb; best.confidence = 0;

for c = 1:numel(candidates)
    bw = candidates(c).mask;
    bw = imopen(bw, strel('rectangle', [1 3]));   % keep thin horizontal-ish
    bw = bw | imopen(candidates(c).mask, strel('rectangle', [3 1])); % or vertical
    CC = bwconncomp(bw);
    if CC.NumObjects == 0, continue; end
    rp = regionprops(CC, 'BoundingBox', 'Extent', 'Area', 'Centroid', ...
                          'MajorAxisLength', 'MinorAxisLength');
    for k = 1:numel(rp)
        bbox = rp(k).BoundingBox;            % [x y w h]
        w = bbox(3); h = bbox(4);
        Lpix = max(w, h); thick = min(w, h);
        if thick < 1, continue; end
        aspect = Lpix / thick;
        ext    = rp(k).Extent;               % rectangularity
        cx = rp(k).Centroid(1); cy = rp(k).Centroid(2);

        % Geometric gates
        if aspect < opt.MinBarAspect,            continue; end
        if ext    < opt.MinRectangularity,       continue; end
        if Lpix   < 0.03*W || Lpix > 0.65*W,     continue; end
        if thick  > 0.06*max(H,W),               continue; end

        % Confidence scoring
        rectScore   = min(ext, 1);
        aspectScore = min(aspect/8, 1);
        sz          = Lpix / W;
        sizeScore   = exp(-((sz - 0.15)/0.12)^2);
        corners = [1 1; W 1; 1 H; W H];
        dmin = min(sqrt(sum(([cx cy] - corners).^2, 2)));
        posScore = 1 - dmin/sqrt(W^2 + H^2);
        nearBorder = (cx < borderF*W || cx > (1-borderF)*W || ...
                      cy < borderF*H || cy > (1-borderF)*H);
        borderBonus = 0.1*double(nearBorder);
        conf = 0.32*rectScore + 0.22*aspectScore + 0.18*sizeScore + ...
               0.18*posScore + borderBonus;

        if conf > best.confidence
            best.found       = true;
            best.confidence  = conf;
            best.bbox        = bbox;
            best.lengthPix   = Lpix;
            best.orientation = ternary(w >= h, 'horizontal', 'vertical');
            best.polarity    = candidates(c).polarity;
        end
    end
end

if ~best.found, sb = best; return; end
sb = best;

% --- OCR to read the bar text (length + unit) ---
pad = round(opt.ScaleBarPadFrac * max(H,W));
bx = sb.bbox;
x1 = max(1, floor(bx(1) - 2.5*pad));   y1 = max(1, floor(bx(2) - 2.5*pad));
x2 = min(W, ceil(bx(1) + bx(3) + 2.5*pad));
y2 = min(H, ceil(bx(2) + bx(4) + 2.5*pad));
roiBox = [x1 y1 (x2-x1) (y2-y1)];

[microns, txt, punit] = readScaleText(gray, roiBox, opt);
sb.text = txt; sb.physicalUnit = punit;
if ~isnan(microns)
    sb.physicalLength  = microns;
    sb.micronsPerPixel = microns / sb.lengthPix;
end

% --- Removal mask: bar + annotation band (padded), inpainted later ---
rm = false(H, W);
rm(max(1,floor(y1)):min(H,ceil(y2)), max(1,floor(x1)):min(W,ceil(x2))) = true;
rm = imdilate(rm, strel('disk', max(1, round(pad/2))));
sb.removalMask = rm;
end


function [microns, txt, unit] = readScaleText(gray, roiBox, opt)
%READSCALETEXT  OCR the annotation near the bar and parse "value unit".
microns = NaN; txt = ''; unit = '';

% Allow a manual override if OCR is disabled or fails.
manual = strtrim(opt.ScaleBarText);

if opt.UseOCR
    try
        sub = imcrop(gray, roiBox);
        sub = imresize(sub, 3);                  % upsample helps OCR
        subA = imbinarize(sub);                  % try both polarities
        cand = {sub, subA, imcomplement(subA)};
        for i = 1:numel(cand)
            r = ocr(cand{i});
            t = strtrim(r.Text);
            if ~isempty(t)
                [m, u] = parseScaleString(t);
                if ~isnan(m), microns = m; unit = u; txt = t; return; end
                if isempty(txt), txt = t; end
            end
        end
    catch
        % OCR (Computer Vision Toolbox) unavailable -> fall through
    end
end

if ~isempty(manual)
    [m, u] = parseScaleString(manual);
    if ~isnan(m), microns = m; unit = u; txt = manual; end
end
end


function [microns, unit] = parseScaleString(s)
%PARSESCALESTRING  Extract a length in microns from text like "50 µm".
microns = NaN; unit = '';
s = strrep(s, char(956), 'u');   % greek mu -> u
s = strrep(s, char(181), 'u');   % micro sign -> u
tok = regexpi(s, '(\d+\.?\d*)\s*(nm|um|mm|m)\b', 'tokens', 'once');
if isempty(tok)
    tok = regexpi(s, '(\d+\.?\d*)\s*(nm|um|mm)', 'tokens', 'once');
end
if isempty(tok), return; end
val = str2double(tok{1});
u   = lower(tok{2});
switch u
    case 'nm', microns = val/1000;  unit = 'nm';
    case 'um', microns = val;       unit = 'um';
    case 'mm', microns = val*1000;  unit = 'mm';
    otherwise, microns = val;       unit = 'um';
end
end


% ---------------------------------------------------------------------
% SEGMENTATION
% ---------------------------------------------------------------------
function [mask, info] = segmentPores(poreImg, enh, validMask, opt)
%SEGMENTPORES  Produce a raw binary pore mask (TRUE = pore).
info = struct();
method = lower(opt.SegmentationMethod);

% Optimize adaptive sensitivity if requested
sens = opt.AdaptiveSensitivity;
if opt.AutoOptimize && any(strcmp(method, {'adaptive','hybrid'}))
    sens = optimizeSensitivity(poreImg, validMask);
end
info.sensitivity = sens;

% Intensity-based masks
adaptT   = adaptthresh(poreImg, sens, 'NeighborhoodSize', ...
                       2*floor(size(poreImg)/16)+1);
maskAdapt = imbinarize(poreImg, adaptT) & validMask;

otsuLevel = graythresh(poreImg(validMask));
maskOtsu  = imbinarize(poreImg, otsuLevel) & validMask;
info.otsuLevel = otsuLevel;

% Edge-based mask (Canny on enhanced image)
if isempty(opt.CannyThreshold)
    edges = edge(enh, 'Canny');
else
    edges = edge(enh, 'Canny', opt.CannyThreshold);
end
edgesClosed = imclose(edges, strel('disk', 2));
maskEdge = imfill(edgesClosed, 'holes') & validMask;
maskEdge = maskEdge & ~edgesClosed;     % keep interior, drop the rim

switch method
    case 'adaptive', mask = maskAdapt;
    case 'otsu',     mask = maskOtsu;
    case 'edge',     mask = maskEdge;
    case 'hybrid'
        % Union of intensity detection with edge-enclosed regions, but only
        % keep edge regions that are also relatively dark (real pores).
        darkish = poreImg > (otsuLevel*0.6);
        mask = maskAdapt | (maskEdge & darkish);
    otherwise
        error('Unknown SegmentationMethod: %s', opt.SegmentationMethod);
end
info.method = method;
end


function sens = optimizeSensitivity(poreImg, validMask)
%OPTIMIZESENSITIVITY  Pick the most stable adaptive sensitivity (plateau).
cand = 0.30:0.05:0.70;
por  = zeros(size(cand));
vp   = nnz(validMask);
for i = 1:numel(cand)
    T  = adaptthresh(poreImg, cand(i));
    bw = imbinarize(poreImg, T) & validMask;
    por(i) = 100 * nnz(bw) / vp;
end
% Stability = smallest local change in porosity (knee/plateau)
d = abs(gradient(por));
d(1) = inf; d(end) = inf;            % avoid endpoints
[~, idx] = min(d);
sens = cand(idx);
end


% ---------------------------------------------------------------------
% MORPHOLOGICAL REFINEMENT
% ---------------------------------------------------------------------
function mask = refineMask(mask, validMask, opt)
%REFINEMASK  Clean, smooth, and separate touching pores.
mask = mask & validMask;

if opt.FillHoles,  mask = imfill(mask, 'holes'); end
if opt.OpenRadius  > 0, mask = imopen(mask,  strel('disk', opt.OpenRadius));  end
if opt.CloseRadius > 0, mask = imclose(mask, strel('disk', opt.CloseRadius)); end

if opt.RemoveScratches
    mask = removeScratches(mask, opt);
end

if opt.MinPoreArea > 0
    mask = bwareaopen(mask, round(opt.MinPoreArea));
end

if opt.SmoothBoundaries
    mask = bwmorph(mask, 'majority');     % light boundary smoothing
end

if opt.Watershed
    mask = watershedSeparate(mask, opt.WatershedDepth);
end

mask = mask & validMask;
if opt.MinPoreArea > 0
    mask = bwareaopen(mask, round(opt.MinPoreArea));
end
end


function mask = removeScratches(mask, opt)
%REMOVESCRATCHES  Drop very thin, very elongated objects (polishing marks).
CC = bwconncomp(mask);
if CC.NumObjects == 0, return; end
rp = regionprops(CC, 'MajorAxisLength', 'MinorAxisLength');
toRemove = false(CC.NumObjects, 1);
for k = 1:numel(rp)
    minor = max(rp(k).MinorAxisLength, eps);
    aspect = rp(k).MajorAxisLength / minor;
    if minor <= opt.ScratchMaxThickness && aspect >= opt.ScratchMinAspect
        toRemove(k) = true;
    end
end
mask(cell2mat(CC.PixelIdxList(toRemove)')) = false;
end


function mask = watershedSeparate(mask, depth)
%WATERSHEDSEPARATE  Split touching convex pores via distance-transform WS.
if ~any(mask(:)), return; end
D  = -bwdist(~mask);
D  = imhmin(D, depth);          % suppress shallow minima (reduce oversplit)
D(~mask) = Inf;
Ld = watershed(D);
mask(Ld == 0) = false;
end


% ---------------------------------------------------------------------
% CHARACTERIZATION + CLASSIFICATION
% ---------------------------------------------------------------------
function [T, L] = characterizePores(mask, mpp, unit, opt)
%CHARACTERIZEPORES  Measure every pore and classify it.
L = bwlabel(mask);
props = regionprops(L, 'Area', 'Perimeter', 'Centroid', ...
    'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', ...
    'Orientation', 'Solidity', 'EquivDiameter');

n = numel(props);
if n == 0
    T = emptyPoreTable();
    return;
end

% Feret diameters (version-safe)
maxFeret = nan(n,1); minFeret = nan(n,1);
try
    fp = regionprops(L, 'MaxFeretProperties', 'MinFeretProperties');
    maxFeret = [fp.MaxFeretDiameter]';
    minFeret = [fp.MinFeretDiameter]';
catch
    maxFeret = [props.MajorAxisLength]';
    minFeret = [props.MinorAxisLength]';
end

% Pull arrays
Area   = [props.Area]';
Perim  = [props.Perimeter]';
Major  = [props.MajorAxisLength]';
Minor  = [props.MinorAxisLength]';
Ecc    = [props.Eccentricity]';
Orient = [props.Orientation]';
Solid  = [props.Solidity]';
Deq    = [props.EquivDiameter]';
Cent   = cat(1, props.Centroid);

% Derived shape descriptors
Circ = 4*pi*Area ./ max(Perim.^2, eps);
Circ = min(Circ, 1);
Aspect = Major ./ max(Minor, eps);

% Unit conversion
Area_u   = Area  * mpp^2;
Perim_u  = Perim * mpp;
Deq_u    = Deq   * mpp;
Major_u  = Major * mpp;
Minor_u  = Minor * mpp;
maxFer_u = maxFeret * mpp;
minFer_u = minFeret * mpp;

% Classification
Class = strings(n,1);
for k = 1:n
    Class(k) = classifyDefect(Circ(k), Aspect(k), Solid(k), Ecc(k), ...
                              Deq_u(k), opt.Class);
end

ID = (1:n)';
T = table(ID, Area_u, Perim_u, Deq_u, maxFer_u, minFer_u, Major_u, Minor_u, ...
          Circ, Solid, Aspect, Ecc, Orient, Cent(:,1), Cent(:,2), Class, ...
    'VariableNames', {'ID', ['Area_' unit '2'], ['Perimeter_' unit], ...
        ['EquivDiameter_' unit], ['MaxFeret_' unit], ['MinFeret_' unit], ...
        ['MajorAxis_' unit], ['MinorAxis_' unit], 'Circularity', 'Solidity', ...
        'AspectRatio', 'Eccentricity', 'Orientation_deg', 'CentroidX', ...
        'CentroidY', 'Class'});
end


function cls = classifyDefect(circ, aspect, solid, ~, deq, th)
%CLASSIFYDEFECT  Rule-based classifier (first match wins).
if aspect >= th.CrackAspect && circ < th.CrackCircularity && solid < th.CrackSolidity
    cls = "Crack";
elseif aspect >= th.LoFAspect && circ < th.LoFCircularity
    cls = "Lack-of-fusion";
elseif circ >= th.GasCircularity && aspect <= th.GasAspectMax && solid >= th.GasSolidity
    if deq >= th.KeyholeDiam
        cls = "Keyhole";
    else
        cls = "Gas pore";
    end
else
    cls = "Irregular";
end
end


function T = emptyPoreTable()
T = table('Size', [0 16], 'VariableTypes', ...
    [repmat({'double'},1,15), {'string'}]);
end


% ---------------------------------------------------------------------
% STATISTICS
% ---------------------------------------------------------------------
function stats = computeStatistics(T, mask, validMask, unit, unitDisp, aunit)
%COMPUTESTATISTICS  Porosity and pore-size statistics.
%   unit     = ASCII token used to index table columns ('um'/'px')
%   unitDisp = pretty unit string stored for display
stats = struct();
stats.unit       = unitDisp;
stats.areaUnit   = aunit;
stats.porosityPct = 100 * nnz(mask) / max(nnz(validMask),1);
stats.numPores    = height(T);

if height(T) == 0
    [stats.meanSize, stats.medianSize, stats.stdSize, stats.minSize, ...
     stats.maxSize] = deal(NaN);
    stats.ci95 = [NaN NaN];
    stats.classCounts = struct();
    stats.diam = [];
    stats.area = [];
    return;
end

d = T.(['EquivDiameter_' unit]);
a = T.(['Area_' unit '2']);
stats.diam = d;  stats.area = a;
stats.meanSize   = mean(d);
stats.medianSize = median(d);
stats.stdSize    = std(d);
stats.minSize    = min(d);
stats.maxSize    = max(d);

% 95% confidence interval of the mean (normal approximation)
nse = stats.stdSize / sqrt(numel(d));
stats.ci95 = stats.meanSize + [-1.96 1.96]*nse;

% Defect-class counts
classes = ["Gas pore","Keyhole","Lack-of-fusion","Crack","Irregular"];
cc = struct();
for i = 1:numel(classes)
    key = char(matlab.lang.makeValidName(classes(i)));
    cc.(key) = sum(T.Class == classes(i));
end
stats.classCounts = cc;
stats.classNames  = classes;
stats.classVec    = arrayfun(@(c) sum(T.Class == c), classes);
end


% ---------------------------------------------------------------------
% SUMMARY PRINT
% ---------------------------------------------------------------------
function printSummary(r)
s = r.stats;
fprintf('\n========================================================\n');
fprintf(' POROSITY ANALYSIS: %s\n', r.imagePath);
fprintf('--------------------------------------------------------\n');
fprintf(' Calibration source : %s\n', r.calSource);
if r.calibrated
    fprintf(' Microns per pixel  : %.5f %s/px\n', r.micronsPerPixel, char(181));
else
    fprintf(' Microns per pixel  : not calibrated (results in pixels)\n');
end
if r.scaleBar.found
    fprintf(' Scale bar          : %s, %s, conf=%.2f, len=%.0f px, text="%s"\n', ...
        r.scaleBar.polarity, r.scaleBar.orientation, r.scaleBar.confidence, ...
        r.scaleBar.lengthPix, r.scaleBar.text);
else
    fprintf(' Scale bar          : not detected\n');
end
fprintf(' Segmentation       : %s (sensitivity=%.2f)\n', ...
        r.segInfo.method, r.segInfo.sensitivity);
fprintf('--------------------------------------------------------\n');
fprintf(' Total porosity     : %.2f %%\n', s.porosityPct);
fprintf(' Number of pores    : %d\n', s.numPores);
if s.numPores > 0
    fprintf(' Mean pore size     : %.3f %s\n', s.meanSize, s.unit);
    fprintf(' Median pore size   : %.3f %s\n', s.medianSize, s.unit);
    fprintf(' Std deviation      : %.3f %s\n', s.stdSize, s.unit);
    fprintf(' Min / Max size     : %.3f / %.3f %s\n', s.minSize, s.maxSize, s.unit);
    fprintf(' 95%% CI of mean     : [%.3f, %.3f] %s\n', s.ci95(1), s.ci95(2), s.unit);
    fprintf('--------------------------------------------------------\n');
    fprintf(' Defect classes:\n');
    for i = 1:numel(s.classNames)
        fprintf('   %-16s : %d\n', s.classNames(i), s.classVec(i));
    end
end
fprintf('========================================================\n');
end


% ---------------------------------------------------------------------
% VISUALIZATION
% ---------------------------------------------------------------------
function visualizeResults(r, opt)
%VISUALIZERESULTS  Two publication-quality figures.
muStr = char(181);
classColors = [ ...
    0.15 0.70 0.20;    % Gas pore   - green
    0.10 0.45 0.95;    % Keyhole    - blue
    0.95 0.60 0.10;    % Lack-of-fusion - orange
    0.90 0.10 0.10;    % Crack      - red
    0.70 0.20 0.85];   % Irregular  - magenta
classNames = ["Gas pore","Keyhole","Lack-of-fusion","Crack","Irregular"];

% ---- FIGURE 1: image pipeline ----
f1 = figure('Name','Porosity Analysis - Images','Color','w', ...
            'Units','normalized','Position',[0.04 0.1 0.92 0.8]);
tiledlayout(f1, 2, 3, 'TileSpacing','compact','Padding','compact');

nexttile; imshow(r.rgb);       title('1. Original Image');

nexttile; imshow(r.gray); hold on; title('2. Scale-bar Detection');
if r.scaleBar.found
    bx = r.scaleBar.bbox;
    rectangle('Position', bx, 'EdgeColor', [1 0 0], 'LineWidth', 2);
    txt = sprintf('%s/%s  c=%.2f', r.scaleBar.polarity, ...
                  r.scaleBar.orientation, r.scaleBar.confidence);
    text(bx(1), max(1,bx(2)-8), txt, 'Color','y', 'FontSize',9, ...
         'FontWeight','bold', 'BackgroundColor',[0 0 0 ]);
else
    text(10,20,'No scale bar detected','Color','y','FontWeight','bold');
end
hold off;

nexttile; imshow(r.enhanced);  title('3. Enhanced Image (CLAHE + bilateral)');

nexttile; imshow(r.binaryMask);title('4. Binary Pore Mask');

nexttile;
overlay = makeOverlay(r.rgb, r.binaryMask, [1 0 0], 0.5);
imshow(overlay); title('5. Detected Pores Overlay');

nexttile;
imshow(makeClassMap(r.rgb, r.labelMatrix, r.poreTable, classNames, classColors));
title('6. Defect Classification Map');
hold on;
for i = 1:numel(classNames)
    plot(NaN,NaN,'s','MarkerFaceColor',classColors(i,:), ...
        'MarkerEdgeColor','none','MarkerSize',10, ...
        'DisplayName',char(classNames(i)));
end
legend('Location','southoutside','Orientation','horizontal','FontSize',7);
hold off;

sgtitle(sprintf('Porosity = %.2f%%   |   %d pores   |   %s', ...
        r.stats.porosityPct, r.stats.numPores, r.imagePath), ...
        'FontWeight','bold');

% ---- FIGURE 2: statistics ----
if r.stats.numPores == 0, return; end
s = r.stats;
f2 = figure('Name','Porosity Analysis - Statistics','Color','w', ...
            'Units','normalized','Position',[0.06 0.08 0.88 0.84]);
tiledlayout(f2, 2, 3, 'TileSpacing','compact','Padding','compact');

% Pore size histogram
nexttile; histogram(s.diam, 'FaceColor',[0.2 0.5 0.9]);
xlabel(sprintf('Equivalent diameter (%s)', s.unit)); ylabel('Count');
title('Pore Size Histogram'); grid on;

% Circularity histogram
nexttile; histogram(r.poreTable.Circularity, 0:0.05:1, 'FaceColor',[0.2 0.7 0.3]);
xlabel('Circularity'); ylabel('Count');
title('Circularity Histogram'); grid on;

% Aspect ratio histogram
nexttile; histogram(r.poreTable.AspectRatio, 'FaceColor',[0.9 0.55 0.1]);
xlabel('Aspect ratio'); ylabel('Count');
title('Aspect Ratio Histogram'); grid on;

% Cumulative pore size distribution (by count and by area)
nexttile;
ds = sort(s.diam); cc = (1:numel(ds))'/numel(ds)*100;
[da, ia] = sort(s.diam); ca = cumsum(s.area(ia))/sum(s.area)*100;
plot(ds, cc, '-o', 'LineWidth',1.4, 'MarkerSize',3); hold on;
plot(da, ca, '-s', 'LineWidth',1.4, 'MarkerSize',3); hold off;
xlabel(sprintf('Equivalent diameter (%s)', s.unit)); ylabel('Cumulative (%)');
title('Cumulative Pore Size Distribution');
legend({'By count','By area'}, 'Location','southeast'); grid on;

% Area-fraction distribution (share of total pore area per size bin)
nexttile;
edges = linspace(min(s.diam), max(s.diam)+eps, 11);
bin = discretize(s.diam, edges);
af = accumarray(bin(~isnan(bin)), s.area(~isnan(bin)), [numel(edges)-1 1]);
af = 100 * af / sum(af);
ctr = edges(1:end-1) + diff(edges)/2;
bar(ctr, af, 'FaceColor',[0.55 0.35 0.75]);
xlabel(sprintf('Equivalent diameter (%s)', s.unit));
ylabel('Area fraction (%)');
title('Area-Fraction Distribution'); grid on;

% Pie chart of defect classes
nexttile;
present = s.classVec > 0;
if any(present)
    labels = arrayfun(@(c,n) sprintf('%s (%d)', c, n), ...
        classNames(present), s.classVec(present), 'UniformOutput', false);
    pp = pie(double(s.classVec(present)));
    % colour the wedges
    pc = pp(1:2:end);
    cidx = find(present);
    for i = 1:numel(pc), pc(i).FaceColor = classColors(cidx(i),:); end
    legend(labels, 'Location','eastoutside', 'FontSize',8);
end
title('Defect Class Composition');

sgtitle(sprintf(['Statistics  |  mean=%.2f %s  median=%.2f %s  ' ...
    'std=%.2f  max=%.2f  95%%CI=[%.2f, %.2f] %s'], ...
    s.meanSize, s.unit, s.medianSize, s.unit, s.stdSize, s.maxSize, ...
    s.ci95(1), s.ci95(2), s.unit), 'FontWeight','bold');
end


function out = makeOverlay(rgb, mask, color, alpha)
%MAKEOVERLAY  Blend a colored mask over an image.
base = im2double(rgb);
if size(base,3) == 1, base = repmat(base,1,1,3); end
out = base;
for c = 1:3
    ch = out(:,:,c);
    ch(mask) = (1-alpha)*ch(mask) + alpha*color(c);
    out(:,:,c) = ch;
end
end


function out = makeClassMap(rgb, L, T, classNames, classColors)
%MAKECLASSMAP  Colour each pore by its defect class.
base = im2double(rgb);
if size(base,3) == 1, base = repmat(base,1,1,3); end
out = base;
if isempty(T) || height(T) == 0, return; end
R = out(:,:,1); G = out(:,:,2); B = out(:,:,3);
for k = 1:height(T)
    ci = find(classNames == T.Class(k), 1);
    if isempty(ci), ci = numel(classNames); end
    col = classColors(ci,:);
    idx = (L == T.ID(k));
    R(idx) = 0.25*R(idx) + 0.75*col(1);
    G(idx) = 0.25*G(idx) + 0.75*col(2);
    B(idx) = 0.25*B(idx) + 0.75*col(3);
end
out = cat(3, R, G, B);
end


% ---------------------------------------------------------------------
% SMALL UTILITY
% ---------------------------------------------------------------------
function out = ternary(cond, a, b)
%TERNARY  Inline conditional helper.
if cond, out = a; else, out = b; end
end

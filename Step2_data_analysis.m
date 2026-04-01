% In previous code, we identify the cell mask for both MED1 and RPB1. and then
% calculate the distance between to cluster; Pearson Coefficient;
% Overlapped region of MED1; Size of MED1 cluster; and number of sub-dots

clear all;
close all;

% Parameters
width = 0.1;  % Gaussian width for rendering
a1 = 20;      % Image scale for Channel 1
a2 = 20;      % Image scale for Channel 2
b = 1;        % Thresholding factor for binarization

% Folder paths
folder_path = '***';
filename = dir(fullfile(strcat([folder_path '\'], '*.bin')));
saving_path = fullfile(folder_path, 'MaskAnalysis');
mkdir(saving_path);

% Initialize results storage
results = [];

% Column headers for Excel output
headers = {'File', 'Size_Ch1_Mask', 'Size_Ch2_Mask', 'PCC', 'Ch2/Ch1_Ratio', 'Ch1_Ch2_Distance', 'Overlap_Area', 'Num_Subdots_Ch1', 'Num_Subdots_Ch2'};

for ii = 1:length(filename)
    % Import data from .bin file
    [MList, ~] = ReadMasterMoleculeList([folder_path '\' filename(ii).name], 'fieldsToLoad', {'xc', 'yc', 'zc', 'c'}, 'ZScale', 160);

    % Extract localizations for Channel 1 & 2
    xc1 = MList.xc(MList.c == 1);
    yc1 = MList.yc(MList.c == 1);
    xc2 = MList.xc(MList.c == 2);
    yc2 = MList.yc(MList.c == 2);

    % Skip if either channel is empty
    if isempty(xc1) || isempty(xc2)
        fprintf('Skipping empty file: %s\n', filename(ii).name);
        continue;
    end

    % Apply DBSCAN to merge nearby puncta
    labels1 = dbscan([xc1, yc1], 0.25, 10);
    labels2 = dbscan([xc2, yc2], 0.25, 10);

    % Retain only clustered points
    xc1 = xc1(labels1 ~= -1);
    yc1 = yc1(labels1 ~= -1);
    xc2 = xc2(labels2 ~= -1);
    yc2 = yc2(labels2 ~= -1);

    % Render images
    ROI = [min(yc1) max(yc1); min(xc1) max(xc1)];
    rendered1 = RenderMList([xc1, yc1], 'gaussianWidth', width, 'ROI', ROI, 'imageScale', a1);
    rendered2 = RenderMList([xc2, yc2], 'gaussianWidth', width, 'ROI', ROI, 'imageScale', a2);

    % Normalize and resize Channel 2 to match Channel 1
    rendered2 = imresize(rendered2, a1/a2);

    % Background subtraction
    rendered1 = rendered1 - imopen(rendered1, strel('disk', 25));
    rendered2 = rendered2 - imopen(rendered2, strel('disk', 25));

    % Binarization
    bw1 = imbinarize(rendered1, graythresh(rendered1) * b);
    bw2 = imbinarize(rendered2, graythresh(rendered2) * b);

    % Remove small objects (noise)
    bw1 = bwareaopen(bw1, 10);
    bw2 = bwareaopen(bw2, 10);

    % Find the largest connected component in each channel
    CC1 = bwconncomp(bw1);
    CC2 = bwconncomp(bw2);

    if CC1.NumObjects > 0
        sizes1 = cellfun(@numel, CC1.PixelIdxList);
        [~, idx1] = max(sizes1);
        bw1_largest = zeros(size(bw1));
        bw1_largest(CC1.PixelIdxList{idx1}) = 1;
    else
        bw1_largest = bw1;
    end

    if CC2.NumObjects > 0
        sizes2 = cellfun(@numel, CC2.PixelIdxList);
        [~, idx2] = max(sizes2);
        bw2_largest = zeros(size(bw2));
        bw2_largest(CC2.PixelIdxList{idx2}) = 1;
    else
        bw2_largest = bw2;
    end

    % Compute centroids
    stats1 = regionprops(bw1_largest, 'Centroid');
    stats2 = regionprops(bw2_largest, 'Centroid');

    if isempty(stats1) || isempty(stats2)
        fprintf('Skipping due to missing centroid: %s\n', filename(ii).name);
        continue;
    end

    centroid1 = stats1.Centroid;
    centroid2 = stats2.Centroid;

    % Calculate Euclidean Distance Between Ch1 and Ch2 Centroids
    dist_ch1_ch2 = sqrt((centroid1(1) - centroid2(1))^2 + (centroid1(2) - centroid2(2))^2);

    % ===== Merge Based on Overlap Condition =====
    overlap = bw1_largest & bw2_largest;
    if any(overlap(:))
        mask = bw1_largest | bw2_largest; % Combined mask
    else
        mask = bw2_largest; % Only Ch2
    end

    % ===== Manual Mask Selection Step =====
    figure;
    imshow(mask);
    title('Select the valid mask region. Press Enter when done.');
    user_mask = roipoly();
    mask = mask & user_mask;

    % Calculate the size of Ch1 and Ch2 **within the selected mask**
    size_ch1_mask = sum(bw1_largest(:) & mask(:));
    size_ch2_mask = sum(bw2_largest(:) & mask(:));

    % Compute Overlapped Area (Absolute Pixel Count)
    mask_ch2 = bw2_largest & mask;
    actual_overlap = mask_ch2 & bw1_largest;
    overlap_area = sum(actual_overlap(:));

    % ===== Detect Sub-dots (High Sensitivity) within the Final Mask =====
    filtered1 = imgaussfilt(rendered1, 1) - imgaussfilt(rendered1, 3);
    filtered2 = imgaussfilt(rendered2, 1) - imgaussfilt(rendered2, 3);

    filtered1_masked = filtered1 .* mask;
    filtered2_masked = filtered2 .* mask;

    bw_dots1 = imbinarize(filtered1_masked, graythresh(filtered1_masked) * 1.2);
    bw_dots2 = imbinarize(filtered2_masked, graythresh(filtered2_masked) * 1.2);

    bw_dots1 = bwareaopen(bw_dots1, 3);
    bw_dots2 = bwareaopen(bw_dots2, 3);

    CC_dots1 = bwconncomp(bw_dots1);
    CC_dots2 = bwconncomp(bw_dots2);

    num_dots_ch1 = CC_dots1.NumObjects;
    num_dots_ch2 = CC_dots2.NumObjects;

    % ===== Generate Four-Panel Output Image =====
    originalImage = cat(3, imadjust(rendered2), imadjust(rendered1), zeros(size(mask)));
    largestRegions = cat(3, bw2_largest, bw1_largest, zeros(size(mask)));
    maskImage = repmat(mask, 1, 1, 3);  % <<< UPDATED: Show mask in white
    dotImage = cat(3, bw_dots2, bw_dots1, zeros(size(mask)));

    combinedImage = [originalImage, largestRegions, maskImage, dotImage];

    figure;
    imshow(combinedImage);

    % Save TIFF
    PCCfile = fullfile(saving_path, sprintf('Mask_Combined_%03d.tif', ii));
    imwrite(combinedImage, PCCfile);

    % Store results
    results = [results; {filename(ii).name, size_ch1_mask, size_ch2_mask, ...
                         corr2(rendered1, rendered2), size_ch2_mask / size_ch1_mask, ...
                         dist_ch1_ch2, overlap_area, num_dots_ch1, num_dots_ch2}];

    fprintf('Processed file: %s\n', filename(ii).name);
end

% Convert results to table and save
results_table = cell2table(results, 'VariableNames', headers);
writetable(results_table, fullfile(saving_path, 'MaskAnalysisResults.csv'));

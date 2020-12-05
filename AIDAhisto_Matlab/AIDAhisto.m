function [peaks, peaks_coord, settings]=AIDAhisto(inputPath,varargin)
%% AIDAhisto: Atlas-based imaging data analysis tool for mouse brain histology
%
% Syntax:  peaks = AIDAhisto(input_Image, 'WIDTH', width, 'CHANNEL', int ,...
%                           'MIN_DIST', double, 'BAR_FILTER', boolean...
%                            'DARK_PEKAS',boolean,'SAVE_DATA'str,...
%                            'ROI_PATH',str, 'SETTINGS', structure)
% Inputs:
%   inputPath  - path to input grayscale image
%                or color image (see parameter 'CHANNEL')                 
%   
% Parameters: 
%   width      - cell width (double, estimation)
%   CHANNEL    - color channel RGB (int) 1 - 3 
%   REF_CHANNEL- background channel that will be used to mark artifacts
%   MIN_DIST   - min distance between cells (double)
%   BAR_FILTER - set to 1 for detecting non circular cells e.g GFAP or IBA1
%   DARK_PEAKS - set to 1 for detecting dark cells in bright background (int)
%   SAVE_DATA  - set to 0 for avoid saving images (.png) and peak coordinates (.txt)
%   THRES_W    - weighting factor for isodata threshold - default 10
%   SETTINGS   - A structure with multiple parameters as its fields
%
% Outputs:
%   peaks       - 2D binary images with peaks
%   peaks_coord - coordinates of peaks
%   settings    - a struct with the settings used to run AIDAhist
%
% NOTE: Input arguments are processed in the order they are submitted, thus
% in the following example, the thres_w will be set to 5:
%
%   >> settings.tres_w = 10;
%   >> AIDAHist(input_image, width, 'SETTINGS', settings, 'THRES_W', 5);
%
%
% Author: Niklas Pallast
% Neuroimaging & Neuroengineering
% Department of Neurology
% University Hospital Cologne
% Kerpener Str. 62
% 50937 Cologne, Germany
% Mar 2019; Last revision: 17-June-2019
%
%
% Adapted by: Johannes de Jong
% UC Berkeley
% Last revision: September 10th, 2019
%
%------------- BEGIN CODE --------------


%% read parameters
tic;

% Set default parameters
settings.width = 5;
settings.ch = 1;
settings.ref_ch = 0;
settings.darkPeaks = 0;
settings.filterType = 0;
settings.thresw = 10;

% Set meta parameters
save_data = 0;
plot_data = 0;

i = 1 ;
while length(varargin)>i
    parameter = varargin{i};
    val = varargin{i+1};
    i = i+2;
    switch upper(parameter)
        case 'WIDTH'
            settings.width = val;
        case 'CHANNEL'
            settings.ch = val;
        case 'REF_CHANNEL'
            settings.ref_ch = val;
        case 'MIN_DIST'
            settings.min_dist = val;
        case 'DARK_PEAKS'
            settings.darkPeaks = val;
        case 'BAR_FILTER'
            settings.filterType = val;
        case 'SAVE_DATA'
            save_data = val;      
        case 'THRES_W'
            settings.thresw = val;
        case 'PLOT_DATA'
            plot_data = val;
        case 'SETTINGS'
            settings = val;
        otherwise
            error(['Unknown parameter name ' parameter])
    end
end



% Error handeling 1, does the file actually exist
if ~exist(inputPath,'file')
    error('Input path does not exist')
end


% Error handeling 2, check if we can actually apply the filters to cells of
% this size. (Convolutions break down below 3 pixels)
if settings.width<3
    error("The Cell width has to be bigger than 3 pixels")     
end


% Collect the output file string in case we want to save
paths = HAN_get_paths(inputPath);


% Load the image, if it's a multy-chanel (color) image, load only the
% selected channel
source = imread(paths.source);
if size(source,3)>1
    input_Image = source(:,:,settings.ch);
end
fprintf('Using channel %i for analysis.\n',settings.ch);

% Subtract a reference if the user asked that
if settings.ref_ch ~=0
    p = polyfit(double(source(:,:,settings.ref_ch)), double(input_Image), 1);
    input_Image = input_Image - (source(:,:,settings.ref_ch).*p(1) + p(2));
    fprintf('Fitted and subtracted channel %i to remove artifacts.\n',settings.ref_ch);
end


% Figure out the minimum distance between the cells
if ~exist('mind_dist','var')
    min_dist = round(settings.width/2);
end


% If the cells are darker then the background, just invert the image.
% However, darkPeaks seems to work really good on bright cells as well. In
% fact, the AIDAhist example code runs darkPeaks on bright cells.
if settings.darkPeaks == 1
    input_Image = imcomplement(input_Image);
end


% Error handeling 3, check if the image used for analysis is correct
if sum(input_Image(:))==0
    error(['There is no image information in channel ' num2str(settings.ch)])
end


%% poof Image size - the calculation is proceeded on a downscale image
x_s = size(input_Image,1);
y_s = size(input_Image,2);
if (x_s*y_s)>100000000
    input_Image =imresize(input_Image,0.2,'nearest');
end


%% apply 2D convolution with 1 of 2 filter sets
% the filters are automatically scaled to the expected cell width
disp("Finding Kernel");
if settings.filterType == 1
    
    kernel = makeBarFilters(settings.width);
    
    imageVec = zeros([size(input_Image'),size(kernel,3)]);
    disp("Convolution")
    for i=1:size(kernel,3)
        imageVec(:,:,i)=conv2(input_Image',kernel(:,:,i),'same');
    end
    image = (mean(imageVec,3));
else
    kernel = makeSfilters(settings.width);
    imageVec = zeros([size(input_Image'),size(kernel,3)]);
    disp("Convolution")
    for i=1:size(kernel,3)
        imageVec(:,:,i)=conv2(input_Image',kernel(:,:,i),'same');
    end
    image = imcomplement(mean(imageVec,3));
end



% Apply the scaled treshold to the image
threshold = settings.thresw*isodataAlgorithm(input_Image);
disp(['Applying a scaled threshold of: ' num2str(threshold)]);
image(image<threshold) = threshold;
image = image - threshold;


% Find the maximums to every pixel above the threshold (presumably the
% center of the cell). By diluting the pixel with the highest value will
% sort of 'drown out' the pixels around it.
disp("Finding Maximums");
neighborhood = filldisc(floor(min_dist));
J = imdilate(image,neighborhood);
peaks= image==J;
peaks(image==0)=0;


% Find the coordiantes of the maximums, note that the pixels with the
% highest values has the same value as the diluted pixels around it in J,
% the other pixels are 'drowned out' by it.
peaks=peaks';
[r,c] = find(peaks==1);
peaks_coord=[c,r]; % X Y from now on instead of r c (which would be Y X).
fprintf("Found %i cells\n", size(peaks_coord,1));


%% save image and peak coordinates if SAVE_DATA=1
se = strel('disk',2);
if save_data==1
    disp("Save data");
    
    
    % Bringing it back to the original size if we worked with huge images.
    if (x_s*y_s)>100000000
        peaks =imresize(peaks,[x_s,y_s]);
        peaks = imerode(peaks,se);
        input_Image = imresize(input_Image,[x_s,y_s]);
    end
    
    % save the cell data to a txt file
    fileID = fopen(paths.cells,'w');
    fprintf(fileID, 'AIDAhisto - Number of detected cells: %i\ncell postions (xy):\n\n',length(peaks_coord));
    for i=1:length(peaks_coord)
        fprintf(fileID,'%i %i\n',peaks_coord(i,1),peaks_coord(i,2));
    end
    fclose(fileID);
    
    % save the settings to a .mat file
    save(paths.AIDAHis_settings, 'settings')
end


%% Plot the data if that is requested by the user
if plot_data==1
    figure;
    imshow(imoverlay(uint8(255*mat2gray(input_Image)),imdilate(peaks,se),'r'),[])
end

end


%% get neighborhood
function y = filldisc(radius)
    %FILLDISC  Returns an array with a disc full of ones.
    %    Y = FILLDISC(R)  Returns a square array of size 2*ROUND(R)+1.
    %    Elements less than R from the central elements are set to 1, rest to
    %    0. This is done fairly efficiently, taking advantage of symmetries.

    % Copyright David Young 2010

    rsq = radius * radius;
    R = round(radius);
    c = R + 1;
    t = 2*R + 1;
    y = zeros(t);

    % points on St George's cross
    y(:, c) = 1;
    y(c, :) = 1;

    for rstart = 1:R
        rend = round(sqrt(rsq - rstart * rstart));
        if rend < rstart
            break
        end

        % points on St Andrew's cross
        y(c+rstart, c+rstart) = 1;
        y(c+rstart, c-rstart) = 1;
        y(c-rstart, c+rstart) = 1;
        y(c-rstart, c-rstart) = 1;

        % fill in octants
        y(c+rstart:c+rend, c+rstart) = 1;
        y(c+rstart:c+rend, c-rstart) = 1;
        y(c-rend:c-rstart, c+rstart) = 1;
        y(c-rend:c-rstart, c-rstart) = 1;
        y(c+rstart, c+rstart:c+rend) = 1;
        y(c-rstart, c+rstart:c+rend) = 1;
        y(c+rstart, c-rend:c-rstart) = 1;
        y(c-rstart, c-rend:c-rstart) = 1;
    end
end


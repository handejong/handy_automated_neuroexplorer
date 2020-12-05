function [output_ref, training_data] = HAN_automatic_reference(source_image, varargin)
%HAN_AUTOMATIC_REFERENCE Automatically scales and transforms a slide from
%the reference atlas (50um or 10um) to the input source_image.
%   Source image should be an image file. Currently HAN_automatic_reference
%   does not select which atlas slice would be the best fit for you. You
%   have to put in the proposed slice using the arguments 'start slice'. So
%   like this:
%
%   >>> HAN_automatic_reference('3.jpg', 'start slice', 201, '50um atlas');
%
%   This will fit the 201th slice of the 50um atlas to your input slice.
%   Make sure to use 'save data' and/or 'plot data' to save and/or plot
%   the results.
%
%   HAN_automatic_reference is part of the Handy Automated Neuroexplorere
%   (HAN). HAN is made by Johannes de Jong at UC Berkeley.
%   j.w.dejong@berkeley.edu.


% Deal with input arguments
atlas_path = 'Allen Reference Atlas/Slides 10um/slide_';
atlas_type = '10um atlas';
save_data = false;
plot_data = false;
channel = 1;
for i=1:length(varargin)
    switch varargin{i}
        
        case 'skipp'
            continue;
            
        case '50um atlas' % Select the 50um atlas
            atlas_path = 'Allen Reference Atlas/Slides 50um/slide_';
            atlas_type = '50um atlas';
            
        case '10um atlas' % Select the 10um atlas
            % This is the default, so doesn't change anyting for now
            atlas_path = 'Allen Reference Atlas/Slides 10um/slide_';
            atlas_type = '10um atlas';
            
        case 'start slice' % set the atlas start slice
            c_atlas_slice = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        case 'channel'
            channel = varargin{i+1};
            varagin{i+1} = 'skipp';
            
        case 'save data' 
            save_data = true;
            
        case 'plot data'
            plot_data = true;
            
        otherwise
            disp('Unknown input argument:')
            varargin{i};
    end
end


% Get the data
paths = HAN_get_paths(source_image);
source = imread(source_image);
if size(source,3)>1
    source = source(:,:,channel);
end
ref = imread([atlas_path num2str(c_atlas_slice) '.tiff']);

% Grab the colormap, add a line for the 'other' group
load('Allen_colormap.mat','cmap')
cmap = [cmap; [0.9 0.9 0.9]];

% Find the threshold between background and tissue
% Note that the corners are usually empty, while the center contains tissue
corners = zeros(25,25,4);
corners(:,:,1) = source(1:25,1:25);
corners(:,:,2) = source(end-24:end,1:25);
corners(:,:,3) = source(end-24:end, end-24:end);
corners(:,:,4) = source(1:25,end-24:end);
half_X = round(size(source,2)/2);
half_Y = round(size(source,1)/2);
center = source(half_Y-25:half_Y+25, half_X-25: half_X+25);
threshold = 0.5 * (mean(center(:)) + mean(corners(:))); 

% Make the shadows
source_shadow = zeros(size(source));
source_shadow(source>=threshold) = 1;
ref_shadow = zeros(size(ref));
ref_shadow(ref~=1) = 1;

% Scale the atlas
scale_factor = sum(source_shadow(:))/sum(ref_shadow(:));
ref_shadow = imresize(ref_shadow, sqrt(scale_factor), 'nearest');
output_ref = imresize(ref, sqrt(scale_factor), 'nearest');

% Perform the allignment
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 3;
optimizer.MaximumIterations = 300;
tform = imregtform(ref_shadow, source_shadow, 'affine', optimizer, metric);

% Transform
outputView = imref2d(size(source));
output_ref = imwarp(output_ref, tform, 'OutputView', outputView);

% Display
if plot_data
    figure
    source_image = repmat(source,[1 1 3]);
    image_projection = image(source_image,'Tag','source');
    image_projection.CData = image_projection.CData*100;
    axis equal; axis off; hold on;
    output_reference_projection = imagesc(output_ref,...
        'Visible','on',...
        'AlphaData',ones(size(output_ref)).*0.2,...
        'Tag','ref');
    axis equal
    colormap(cmap)
end

% Store the trainign_data
training_data.source = source_image;
training_data.atlas = atlas_type;
training_data.ref_nr = c_atlas_slice;
training_data.atlas_landmarks = [];
training_data.source_landmarks = [];
training_data.tform = tform;
training_data.optimizer = optimizer;
training_data.scale_factor = scale_factor;


% Save the data if requested
if save_data
    % Corrected refference
    imwrite(output_ref, paths.ref);
    save(paths.training_data, 'training_data');
    
    % Update the atlas data
    HAN_update_atlas_data(source_image, training_data);
end




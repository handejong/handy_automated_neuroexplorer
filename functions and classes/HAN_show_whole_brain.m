function flag = HAN_show_whole_brain(whole_brain_file, varargin)
%HAN_SHOW_WHOLE_BRAIN Displays the whole brain data file in 3D
%   After running 

flag = 1; % not completed

% Deal with input arguments
plot_data = false;
save_data = false;
level = 1;
min_cells = 1;
for i=1:length(varargin)
    switch varargin{i}
        
        case 'skipp'
            continue;
            
        case 'level'
            level = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        case 'min cells'
            min_cells = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        otherwise
            disp('Unknown input argument:')
            varargin{i};
    end
end


% Get the data
if ischar(whole_brain_file)
    file_ID = fopen(whole_brain_file);
    data = textscan(file_ID, '%f %f %f %f %f %f %f');
    fclose(file_ID);
    X = data{1};
    Y = data{2};
    Z = data{3};
    region = data{6};
else
    X = whole_brain_file(:,1);
    Y = whole_brain_file(:,2);
    Z = whole_brain_file(:,3);
    region = whole_brain_file(:,6);
end


% Exclude cells that are outside the brain
indexer = region==1;
disp(['Not including ' num2str(sum(indexer)) ' cells that are outside the brain.']);
region = region(~indexer);
X = X(~indexer);
Y = Y(~indexer);
Z = Z(~indexer);


% Update the region to the correct level and label all cells in regions
% with to few cells as 'other'.
[region, unique_regions, region_names] = HAN_update_level(region, level, min_cells);


% Grab the colormap, add a line for the 'other' group
load('Allen_colormap.mat','cmap')
cmap = [cmap; [0.9 0.9 0.9]];


% Plot all the regions seperate
main_figure = figure('Position',[10, 10, 1000, 800]);
set(gcf, 'Visible','on')
axis ij
% background_image = imread('sagital_50um.tiff');
% image(background_image);
% colormap(cmap);
hold on
axis equal
hold on
for i = 1:length(unique_regions)
    index = region == unique_regions(i);
    
    % Find the displayname
    displayName = region_names{i};

    % Plot the data
    scatter3(X(index), Y(index), Z(index),'.',...
        'DisplayName', displayName,...
        'SizeData', 1,...
        'CData',cmap(unique_regions(i),:));
    drawnow();
    
    % Update the user how much we have done
    disp([num2str(round(100*i/length(unique_regions))) '% done'])
end

% Some formatting on the plot
set(gca,'Color','k')
hAxis = gca;
hAxis.Position(1) = 0.01;

% Some formatting on the legend
legend;
hLegend = findobj(gcf, 'Type', 'Legend');
hLegend.Color = [1 1 1];
hLegend.Position(1) = 0.7;
hLegend.Position(2) = max([hLegend.Position(2), 0.1]);

% Well, apparently we made it through the end
flag = 0;

end


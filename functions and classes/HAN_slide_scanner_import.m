function data = HAN_slide_scanner_import( filename, varargin)
%SLIDE_SCANNER_IMPORT Imports a set of sections from one slide.
%   Filename should be a .czi file. Usually .czi files store the slides at
%   different resolutions. The output results is a stuct with n cells.
%   Where n is the number of channels in the picture. Every cell will
%   contain a 3D matrix with all sections in the 3rd dimension. So the
%   matrix will be (y resolution, x resolution, # of sections).
%
%   This is an overview optional arguments
%           'tiff export'       Will export pixtures to .tiff files
%           'cropper'           Will run a cropper that will allow the user
%                               to crop figure (so save space).
%           'start number'      will name the ouput files (if requested)
%                               starting with the number that is in the 
%                               next argument
%           'resolution'        The resolution (from 1=high to 4 = low)
%                               that whould be used.
%
%   Example:     data=slide_scanner_import('filename.czi','tiff export');
%   Example2:    data=slide_scanner_import('filename.czi','tiff export','cropper','start_number',28);
%
%   note: To open a file browser dialog, type a random filename
%
%   Slide_scanner_import is made by Han de Jong (j.w.dejong@berkeley.edu)
%   Microscopy data is imported using the Bioformat Importer API, which is
%   copyrighed by the Open Microscopy Environment and freely available
%   here: www.openmicroscopy.org.


% Baseline variables
settings.tiff_export = false;         % Will not export tif files
settings.cropper = false;             % Will not run the figure cropper
settings.same_size = false;           % Will not make .tiff files same size
settings.start_number = 1;            % Start number for ordering the slices
settings.resolution = 0;              % Default will use a GUI to decide

% Dealing with input arguments
skip_var=-1;
for i=1:length(varargin)
    argument=varargin{i};
    
    if i~=skip_var
        switch argument
            case 'tiff export'
                disp('Exporting figures to tiff...')
                settings.tiff_export = true;
            case 'cropper'
                disp('Running the cropper...')
                settings.cropper = true;
            case 'start number'
                settings.start_number = varargin{i+1};
                skip_var = i+1;
            case 'resolution'
                settings.resolution = varargin{i+1};
                skip_var = i+1;
            otherwise
                disp('Argument not recognized, please see overview below:')
                disp(' ')
                help slide_scanner_import
        end
    end
end


% importing the file and storing the selected picture in a cell struct
% named raw data. Figuring out the resolution structure
raw_data = bfopen(filename);
counter = 1;
i = size(raw_data{1,1}{1,1});
i2 = size(raw_data{2,1}{1,1});

% Deal with different color channels
nr_channels = size(raw_data{1,1},1);

% Count the different resolutions
while i2<i
    counter=counter+1;
    i=i2;
    i2=size(raw_data{counter+1,1}{1,1});
end
resolution=[1, counter];


% let the user select a resolution and save only those images unless the
% resolution if preset.
if settings.resolution == 0
    res_figure=figure('Name','Click on the prefered resolution.');
    res_figure.Position(3)=1200;
    for i=1:resolution(2)
        subplot(1,4,i)
        average=mean2(raw_data{i,1}{1,1}); %this is for auto contrast
        bd_text=['set(gcf,''UserData'',''',num2str(i),'''); uiresume(gcbf)'];
        res_image{i}=imagesc(raw_data{i,1}{1,1},[1, 5*average]);
        res_image{i}.ButtonDownFcn=bd_text;
        axis equal tight
        colormap gray
    end
    disp('Click on the prefered resolution.')
    disp('Note: not shown at real contrast.')
    disp('Note2: feel free to use figure controls (zoom etc.) but these need to be ''closed'' (clicked again) before selecting one resolution.')
    uiwait(res_figure)
    resolution(1) = str2num(res_figure.UserData);
    close(res_figure)
else
    resolution(1) = settings.resolution;
end

np_sections=(length(raw_data)-2)/resolution(2);
counter=1;
for i=resolution(1):resolution(2):length(raw_data)-2
    % merge channels
    temp = size(raw_data{i,1}{1,1});
    if nr_channels == 1
        data{counter} = zeros(temp);
    else
        data{counter} = zeros([temp 3]); % It has to be three channels
    end
    for j=1:nr_channels
        data{counter}(:,:,j) = raw_data{i,1}{j,1};
    end
    counter=counter+1;
end


% Does the user want to run the cropper?
if settings.cropper
    for i=1:np_sections
        data{i}=cropper(data{i});
    end 
end
    

% Exports Tiff's if that's what the user wants
if settings.tiff_export
    
    % Make a folder where we will store the output
    [~, project_name] = fileparts(filename);
    try
        mkdir(project_name)
    catch
        project_name=inputdlg('Project name:');
        project_name=project_name{1};
        mkdir(project_name)
    end
    
   %Save every image
   for i=1:np_sections
        filename=[project_name '/' num2str(settings.start_number) '.tiff'];
        imwrite(uint16(data{i}),filename)
        settings.start_number=settings.start_number+1;
   end
    
    disp('Data saved to .tiff files')
end
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cropped_figure]=cropper(input_figure)
% Will allow the user to crop the image

main_figure=figure();
mouser.empty=true;
main_figure.UserData=mouser;
top_colormap=mean2(input_figure)*5; 
image=imagesc(input_figure(:,:,1),[1, top_colormap]);
image.Tag='image';
colormap gray
axis equal tight
hold on
axis manual

%Callbacks
set(image,'ButtonDownFcN',@mouse_down);

%Wait for the user to press a button (presumably after cropping)
done_button=uicontrol(main_figure,'String','done','Callback','uiresume(gcbf)');
uiwait(main_figure)
mouser=main_figure.UserData;
close(main_figure)

dims=mouser.selection;
cropped_figure=input_figure(dims(2):dims(2)+dims(4),dims(1):dims(1)+dims(3), :);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mouse_down(scr, ev)
% Controls mouse clicks
mouser=scr.Parent.Parent.UserData; % thiw works for all clicked objects


if ev.Button==1 %left mouse click
    % Delete previous drawing
    switch scr.Tag % what was clicked?
        case 'image' % clicked on image, draw a new rectangle  
            if ~mouser.empty % there is not rectangle currently
                delete(mouser.drawing)
                delete(mouser.ancors.north)
                delete(mouser.ancors.west)
                delete(mouser.ancors.south)
                delete(mouser.ancors.east)
            end

            % Coordinates of the click
            punter=get(scr.Parent,'CurrentPoint');

            mouser.release=false;
            mouser.start=[punter(1,1), punter(1,2)];
            mouser.drawing=rectangle('Position',[punter(1,1), punter(1,2), 1, 1],'EdgeColor',[0.8 0.4 0]);
            % I mixed up east and west! But at least I'm consistent in the
            % whole script.
            mouser.ancors.north=plot(punter(1,1),punter(1,2),'Marker','s','MarkerSize',20,'MarkerEdgeColor','g','ButtonDownFcn',@mouse_down,'Tag','north');
            mouser.ancors.west=plot(punter(1,1),punter(1,2),'Marker','s','MarkerSize',20,'MarkerEdgeColor','g','ButtonDownFcn',@mouse_down,'Tag','west');
            mouser.ancors.south=plot(punter(1,1),punter(1,2),'Marker','s','MarkerSize',20,'MarkerEdgeColor','g','ButtonDownFcn',@mouse_down,'Tag','south');
            mouser.ancors.east=plot(punter(1,1),punter(1,2),'Marker','s','MarkerSize',20,'MarkerEdgeColor','g','ButtonDownFcn',@mouse_down,'Tag','east');
            mouser.figure=scr;
            mouser.empty=false;
            mouser.action='standard';          
        case 'north' %clicked on one of the markers
            mouser.action='north';
        case 'west'
            mouser.action='west';
        case 'south'
            mouser.action='south';
        case 'east'
            mouser.action='east';
        otherwise
            disp('no idea what the mouse is doing')
    end
    scr.Parent.Parent.UserData=mouser;
    scr.Parent.Parent.WindowButtonMotionFcn=@mouse_motion;
    scr.Parent.Parent.WindowButtonUpFcn=@mouse_release;
end
end

function mouse_motion(src, ev)
% Deals with mouse motion

mouser=src.UserData;
punter=get(mouser.figure.Parent,'CurrentPoint'); %this is to get the coordinates on the figure axes

switch mouser.action
    case 'standard'
        block=[mouser.start(1), mouser.start(2), punter(1,1)-mouser.start(1), punter(1,2)-mouser.start(2)];
    case 'north'
        block=mouser.drawing.Position;
        block(4)=block(4)-(punter(1,2)-block(2));
        block(2)=punter(1,2);
    case 'west' % still flipped with east
        block=mouser.drawing.Position;
        block(3)=punter(1,1)-mouser.start(1);
    case 'south'
        block=mouser.drawing.Position;
        block(4)=punter(1,2)-mouser.start(2);
    case 'east' % still flipped with west
        block=mouser.drawing.Position;
        block(3)=block(3)-(punter(1,1)-block(1));
        block(1)=punter(1,1);
    otherwise
        disp('This is an ignored and really strange error.')
end

% Deal with rectangles drawn up or to the left.
if block(3)<=0
    block(1)=block(1)+block(3);
    block(3)=block(3)*-1;   
end
if block(4)<=0
    block(2)=block(2)+block(4);
    block(4)=block(4)*-1;
end
mouser.drawing.Position=block;

% Again, note east and west are flipped
set(mouser.ancors.north,'XData',block(1)+block(3)/2,'YData',block(2))
set(mouser.ancors.west,'XData',block(1)+block(3),'YData',block(2)+block(4)/2)
set(mouser.ancors.south,'XData',block(1)+block(3)/2,'YData',block(2)+block(4))
set(mouser.ancors.east,'XData',block(1),'YData',block(2)+block(4)/2)
            
end

function mouse_release(scr, ev)
% Deals with mouse motion end
scr.WindowButtonMotionFcn=[];
scr.WindowButtonUpFcn=[];

mouser=scr.UserData;
mouser.release=true;
mouser.selection=round(mouser.drawing.Position); %[X Y W H]

scr.UserData=mouser;

end


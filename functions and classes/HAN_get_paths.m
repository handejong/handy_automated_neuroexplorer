function [paths] = HAN_get_paths(source_file)
%HAN_GET_PATHS Returns the paths and filenames to each file associated with
%the input file.
%
%   During processing of source_files the neuroexplorer creates several
%   files at different stages of the process. These are:
%   filename.jpeg (or tiff, etc.) The original microsopy image.
%   filename_ref.tiff  The transformed reference atlas.
%   filename_training_data.mat Data related to the atlas transformaiton
%   filename_AIDAhis_settings.mat The AIDAhis settings
%   filename_AIDAhis.txt Identified cells in the input image
%   filename_AIDAhis_eddited.txt The same as #3, but eddited
%   filename_cells_mapped.txt The cells mapped to the atlas
%
%   HAN_get_paths finds the paths to these associated files and stores them
%   in a handy struct.
%
%   HAN_GET_PATHS is part of the Handy Automated Neuroexplorer (HAN). HAN
%   is made by Johannes de Jong at UC Berkeley. j.w.dejong@berkeley.edu

% The source file
paths.source = source_file;

% Brake up the file parts
[path, filename] = fileparts(source_file);

% The reference
paths.ref = [path '/' filename '_ref.tiff'];

% The training data
paths.training_data = [path '/' filename '_training_data.mat'];

% The atlas data
paths.atlas_data = [path '/' 'atlas_data.mat'];

% The convolution settings
paths.AIDAHis_settings = [path '/' filename '_AIDAHis_settings.mat'];

% The cells
paths.cells = [path '/' filename '_AIDAHis.txt'];

% The eddited cells
paths.cells_edited = [path '/' filename '_AIDAHis_eddited.txt'];

% The cells mapped to the reference
paths.cells_mapped = [path '/' filename '_cells_mapped.txt'];



end


function HAN_plot_cells(cell_file, reference)
    %HAN_PLOT_CELLS will plot the found cells on top of the refference image.
    
    % Get the refference image
    ref_image = imread(reference);
    
    % Get the cells
    file_ID = fopen(cell_file);
        cells_temp = textscan(file_ID, '%u %u', 'HeaderLines', 3);
    fclose(file_ID);
    cells(:,1) = cells_temp{1};
    cells(:,2) = cells_temp{2};
    
    % Make a figure
    figure
    imagesc(ref_image);
    caxis([1 1000])
    hold on
    scatter(cells(:,1), cells(:,2),'.','r');

end
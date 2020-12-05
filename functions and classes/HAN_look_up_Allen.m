function [region_name, region_atlas_id, no_find] = HAN_look_up_Allen(input_ID, level)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

region_name = {};
region_ID = nan;
no_find = [];

legend = load('Allen Reference Atlas/Allen_structure_tree.mat','structure_tree');
legend = legend.structure_tree;

for i=1:length(input_ID)
    % Find the input ID in the legend
    index = input_ID(i);
    try
        c_level = legend(index).depth;
    catch
        if index==1328
            region_name{i} = 'other';
        else
            disp(['Can''t find ' num2str(input_ID(i))])
            region_name{i} = 'unkown';
            region_atlas_id(i) = index;
            no_find(end+1) = input_ID(i);
        end
        continue
    end
    index = input_ID(i);
    c_level = legend(index).depth;
    % Go up the tree to the right level
    while c_level > level
        index = find([legend.id] == legend(index).parent_structure_id);
        c_level = legend(index).depth;
    end

    % Put out the name and the id at this level
    region_name{i} = legend(index).name;
    region_atlas_id(i) = legend(index).atlas_id;
end
        

end


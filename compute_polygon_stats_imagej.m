function compute_polygon_stats_imagej

% Resolution of image
resolution = 542/500 * 1000; % 542 pixels per 500 nm to pixels per micron

[ui_filename,ui_pathname] = uigetfile({'*.zip'},'Select the roi zip file');

file_to_read = [ui_pathname, ui_filename];

[~,filename_double_extension,~] = fileparts(file_to_read);
[~,filename,~] = fileparts(filename_double_extension);

% Add external libraries to path
addpath(genpath('External_Libraries'))

% Get the ROIs
cvsROIs = ReadImageJROI(file_to_read);

% Fix coordinates
for i = 1:length(cvsROIs)
    cvsROIs{i}.mnCoordinates = (cvsROIs{i}.mnCoordinates + 0.5)/resolution;
end

% Compute centroids of the polygons
polygon_centroids = zeros(length(cvsROIs), 2);
for i = 1:length(cvsROIs)
    polygon_centroids(i, :) = polygonCentroid(cvsROIs{i}.mnCoordinates);
%     if isnan(polygon_centroids(i, 1))
%         polygon_centroids(i, :) = xml_struct(i).points(1,:);
%     end
end

% Iterate through the list of all the polygons and make a new list of
% matched polygons
unpaired_polygon = []; % List to contain any unmatched polygons
unpaired_polygon_counter = 0; % Unmatched polygon count
is_paired = zeros(length(cvsROIs), 1); % Whether we have paired the current polygon with another yet and its position
polygon_struct = repmat(struct(), floor(length(cvsROIs)/2), 1); % Structure array to contain the matched polygons
completed_pairs_counter = 0; % How many polygon pairs we have matched
for i = 1:length(cvsROIs)
    % If we have already paired this polygon, skip it
    if is_paired(i, 1)
        continue
    end
    
    % Iterate through all the polygons to look for pairs
    matching_idx = -1;
    for j = 1:length(cvsROIs)
        % If we are comparing the polygon to itself, skip
        if i == j
            continue
        end
        
        % Check if the current polygon i is inside any of the other
        % polygons
        point_in_polygon = inpolygon(cvsROIs{i}.mnCoordinates(:, 1), cvsROIs{i}.mnCoordinates(:, 2), cvsROIs{j}.mnCoordinates(:, 1), cvsROIs{j}.mnCoordinates(:, 2));
        if any(point_in_polygon)
            % If it is output that it is matching and set that i is the
            % inner polygon
            matching_idx = j;
            is_current_inner = true;
            break
        end
        % Check if any polygons are inside this polygon i
        point_in_polygon = inpolygon(cvsROIs{j}.mnCoordinates(:, 1), cvsROIs{j}.mnCoordinates(:, 2), cvsROIs{i}.mnCoordinates(:, 1), cvsROIs{i}.mnCoordinates(:, 2));
        if any(point_in_polygon)
            % If there are, set that we have a match and note that this is
            % the outer polygon
            matching_idx = j;
            is_current_inner = false;
            break
        end
    end
    
    % Check that the polygon we are trying to match is already paired or if
    % we found no pair
    if (is_paired(matching_idx, 1) || matching_idx == -1)
        unpaired_polygon(unpaired_polygon_counter+1).polygon = cvsROIs{i}.mnCoordinates; %#ok<AGROW>
        unpaired_polygon(unpaired_polygon_counter+1).centroid = polygon_centroids(i, :); %#ok<AGROW>
        unpaired_polygon(unpaired_polygon_counter+1).possible_pair = matching_idx; %#ok<AGROW>
        unpaired_polygon_counter = unpaired_polygon_counter+1;
        continue
    end
    
    if is_current_inner
        % If we are the inner polygon, record data in the matched polygon
        % array as such
        polygon_struct(completed_pairs_counter+1).inner_polygon = cvsROIs{i}.mnCoordinates;
        polygon_struct(completed_pairs_counter+1).outer_polygon = cvsROIs{matching_idx}.mnCoordinates;
        polygon_struct(completed_pairs_counter+1).inner_centroid = polygon_centroids(i, :);
        polygon_struct(completed_pairs_counter+1).outer_centroid = polygon_centroids(matching_idx, :);
    else
        % If we are the outer polygon, record data in the matched polygon
        % array as such
        polygon_struct(completed_pairs_counter+1).inner_polygon = cvsROIs{matching_idx}.mnCoordinates;
        polygon_struct(completed_pairs_counter+1).outer_polygon = cvsROIs{i}.mnCoordinates;
        polygon_struct(completed_pairs_counter+1).inner_centroid = polygon_centroids(matching_idx, :);
        polygon_struct(completed_pairs_counter+1).outer_centroid = polygon_centroids(i, :);
    end
    
    % Set that current polygon has been matched
    is_paired(i, 1) = true;
    is_paired(matching_idx, 1) = true;
    
    % Increase the number of polygons matched
    completed_pairs_counter = completed_pairs_counter+1;
end

% Possibly shrink the array in case we haven't matched the expected number
polygon_struct = polygon_struct(1:completed_pairs_counter);

raw_image = imread([ui_pathname, filename_double_extension, '.tif']);

% Plot the matched polygons
fig_num = figure;
imshow(raw_image)
color_order = get(gca,'ColorOrder');
colour_counter = 1;
hold on
for i = 1:length(polygon_struct)
    % Plot inner and outer using the same colour
    plot(polygon_struct(i).inner_polygon([1:end,1],1)*resolution, polygon_struct(i).inner_polygon([1:end,1],2)*resolution, polygon_struct(i).outer_polygon([1:end,1],1)*resolution, polygon_struct(i).outer_polygon([1:end,1],2)*resolution, 'Color', color_order(mod(colour_counter,7)+1, :))
    % Add a text label at the inner polygon's centroid
    text(polygon_struct(i).inner_centroid(1)*resolution, polygon_struct(i).inner_centroid(2)*resolution, num2str(i), 'Color', color_order(mod(colour_counter,7)+1, :))
    colour_counter = colour_counter + 1;
end
hold off

% Set axis to be equal and set axis labels
axis equal
xlabel('x (micron)')
ylabel('y (micron)')

% Plot any unmatched polygons
original_polygon_struct_length = length(polygon_struct);
if unpaired_polygon_counter > 0
    % Warn that we have unmatched polygons
    fprintf(sprintf('Warning! Has %d unpaired polygons.\\n', unpaired_polygon_counter))
    hold on
    for i = 1:length(unpaired_polygon)
        plot(unpaired_polygon(i).polygon([1:end,1],1)*resolution, unpaired_polygon(i).polygon([1:end,1],2)*resolution, 'Color', color_order(mod(colour_counter,7)+1, :))
        % Put text specifying it is unmatched
        text(unpaired_polygon(i).centroid(1)*resolution, unpaired_polygon(i).centroid(2)*resolution, sprintf('Unpaired %d', original_polygon_struct_length + i), 'Color', color_order(mod(colour_counter,7)+1, :))
        colour_counter = colour_counter + 1;
        
        % Add this to the end of the polygon_struct
        polygon_struct(end+1).inner_polygon = unpaired_polygon(i).polygon; %#ok<AGROW>
        polygon_struct(end).inner_centroid = unpaired_polygon(i).centroid;
        polygon_struct(end).outer_polygon = unpaired_polygon(i).polygon;
        polygon_struct(end).outer_centroid = unpaired_polygon(i).centroid;
    end
    hold off
end

% Save the figure
savefig(fig_num, [filename, ' - Figure.fig'])

% Compute the polygon statistcs
for i = 1:length(polygon_struct)
    % Compute the areas
    polygon_struct(i).inner_area = abs(polygonArea(polygon_struct(i).inner_polygon));
    polygon_struct(i).outer_area = abs(polygonArea(polygon_struct(i).outer_polygon));
    
    % Find the matching second order moment ellipse
    polygon_struct(i).inner_ellipse = polygonEquivalentEllipse(polygon_struct(i).inner_polygon);
    polygon_struct(i).outer_ellipse = polygonEquivalentEllipse(polygon_struct(i).outer_polygon);
    
    % Record the major and minor axis of the inner ellipse
    polygon_struct(i).inner_major_axis = polygon_struct(i).inner_ellipse(3)*2;
    polygon_struct(i).inner_minor_axis = polygon_struct(i).inner_ellipse(4)*2;
    
    % Record the angle of the inner ellipse
    temp_angle = polygon_struct(i).inner_ellipse(5);
    % Convert from -90 to 90 to 0 to 180 angle
    if temp_angle < 0
        temp_angle = 180 + temp_angle;
    end
    polygon_struct(i).inner_ellipse_angle = temp_angle;
    
    % Record the major and minor axis of the outer ellipse
    polygon_struct(i).outer_major_axis = polygon_struct(i).outer_ellipse(3)*2;
    polygon_struct(i).outer_minor_axis = polygon_struct(i).outer_ellipse(4)*2;
    
    % Record the angle of the outer ellipse
    temp_angle = polygon_struct(i).outer_ellipse(5);
    % Convert from -90 to 90 to 0 to 180 angle
    if temp_angle < 0
        temp_angle = 180 + temp_angle;
    end
    polygon_struct(i).outer_ellipse_angle = temp_angle;
    
    % Find the ferret distance via brute force (maximum pairwise distance)
    dx = polygon_struct(i).inner_polygon(:,1) - polygon_struct(i).inner_polygon(:,1).';
    dy = polygon_struct(i).inner_polygon(:,2) - polygon_struct(i).inner_polygon(:,2).';
    pairwise_dist = hypot(dx,dy);
    [max_dist,idx] = max(pairwise_dist(:));
    [k1,k2] = ind2sub(size(pairwise_dist),idx);
    polygon_struct(i).inner_feret_diameter = max_dist;
    % Find the ferret angle
    temp_angle = rad2deg(angle2Points(polygon_struct(i).inner_polygon(k2,:), polygon_struct(i).inner_polygon(k1,:)));
    temp_angle = temp_angle - (temp_angle > 180)*180; % Convert to 0 to 180
    polygon_struct(i).inner_feret_angle = temp_angle;
    
    % Find the ferret distance via brute force (maximum pairwise distance)
    dx = polygon_struct(i).outer_polygon(:,1) - polygon_struct(i).outer_polygon(:,1).';
    dy = polygon_struct(i).outer_polygon(:,2) - polygon_struct(i).outer_polygon(:,2).';
    pairwise_dist = hypot(dx,dy);
    [max_dist,idx] = max(pairwise_dist(:));
    [k1,k2] = ind2sub(size(pairwise_dist),idx);
    polygon_struct(i).outer_feret_diameter = max_dist;
    % Find the ferret angle
    temp_angle = rad2deg(angle2Points(polygon_struct(i).outer_polygon(k2,:), polygon_struct(i).outer_polygon(k1,:)));
    temp_angle = temp_angle - (temp_angle > 180)*180; % Convert to 0 to 180
    polygon_struct(i).outer_feret_angle = temp_angle;
    
    try %#ok<TRYNC>
        % Find the minimum feret distance and angle for inner polygon
        [polygon_struct(i).inner_minimum_feret_diameter, temp_angle]= minimumCaliperDiameter(polygon_struct(i).inner_polygon);
        polygon_struct(i).inner_minimum_feret_angle = rad2deg(temp_angle);
    end
    % Find the minimum feret distance and angle for outer polygon
    try  %#ok<TRYNC>
        [polygon_struct(i).outer_minimum_feret_diameter, temp_angle] = minimumCaliperDiameter(polygon_struct(i).outer_polygon);
        polygon_struct(i).outer_minimum_feret_angle = rad2deg(temp_angle);
    end
end

% hold on
% for i = 1:length(polygon_struct)
%     % Plot inner and outer using the same colour
%     drawEllipse(polygon_struct(i).inner_ellipse, 'Color', color_order(mod(i,7)+1, :))
%     drawEllipse(polygon_struct(i).outer_ellipse, 'Color', color_order(mod(i,7)+1, :))
% end
% hold off



% Convert to table
polygon_table = struct2table(polygon_struct, 'AsArray', true);
% Write out a table omitting the polygons and centroids and ellipses but
% including a first row with index numbers
if exist([filename, ' - Stats.xlsx'], 'file')
    delete([filename, ' - Stats.xlsx']);
end
writetable([array2table((1:size(polygon_table,1)).'), polygon_table(:, [5:6, 9:end])], [filename, ' - Stats.xlsx'])

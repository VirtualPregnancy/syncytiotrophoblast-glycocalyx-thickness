function compute_wall_distance_shadow_coverage

% Resolution of image
resolution = 542/500 * 1000; % 542 pixels per 500 nm to pixels per micron

% Tolerance 
coverage_tolerance = 1.25E-6;

% Change this to change where the temporary figures are saved
temp_folder = 'Temp\';

% Set to false to not plot one by one and just compute everything and
% output the xlsx
is_iterate_plot = true;

% Numerical parameters, ideally should check method convergence by doubling
% these parameters and ensuring results do not change
num_line_divisions = 2; % How subdivisions per polygon line
n_even = 100000; % Number of points for histogram, median etc

% Add external libraries to path
addpath(genpath('External_Libraries'))




% Make the tempoary directory
if ~exist(temp_folder, 'dir')
    % Directory does not exist so create it
    mkdir(temp_folder)
end

% Get the figure file path
[ui_filename,ui_pathname] = uigetfile({'*.fig'},'Select the fig file');

file_to_read = [ui_pathname, ui_filename];

[~,filename_double_extension,~] = fileparts(file_to_read);
[~,filename,~] = fileparts(filename_double_extension);

% Open the figure
input_fig = openfig(file_to_read);

% Find the polygons
h = findobj(input_fig, 'Type', 'line');

% Convert the polygon data to structures, also scale by the resolution
xml_struct = struct;
for i = 1:length(h)
    xml_struct(i).points(:, 1) = h(i).XData/resolution;
    xml_struct(i).points(:, 2) = h(i).YData/resolution;
end

% Also get the text so we can get the polygon numbers
h_text = findobj(input_fig, 'Type', 'text');
h_text_position = vertcat(h_text.Position)/resolution;
h_text_position = h_text_position(:, 1:2);


% Compute centroids of the polygons so that we can match them with each
% other
polygon_centroids = zeros(length(xml_struct), 2);
for i = 1:length(xml_struct)
    polygon_centroids(i, :) = polygonCentroid(xml_struct(i).points);
    %     if isnan(polygon_centroids(i, 1))
    %         polygon_centroids(i, :) = xml_struct(i).points(1,:);
    %     end
end

% Work out the distance between all pairs of centroids and find the one
% with minimum distance, this should match inner and outer polygons
dx = polygon_centroids(:,1) - polygon_centroids(:,1).';
dy = polygon_centroids(:,2) - polygon_centroids(:,2).';
dx(1:1+size(dx,1):end) = Inf;
dy(1:1+size(dy,1):end) = Inf;
pairwise_dist = hypot(dx,dy);
[~, min_dist_idx] = min(pairwise_dist);

% Iterate through the list of all the polygons and make a new list of
% matched polygons
unpaired_polygon = struct; % List to contain any unmatched polygons
unpaired_polygon_counter = 0; % Unmatched polygon count
is_paired = zeros(length(xml_struct), 1); % Whether we have paired the current polygon with another yet and its position
polygon_struct = repmat(struct(), floor(length(xml_struct)/2), 1); % Structure array to contain the matched polygons
completed_pairs_counter = 0; % How many polygon pairs we have matched
failed_match_label = []; % Array for labels that we fail to match the centroid exactly and go to closest
for i = 1:length(xml_struct)
    % If we have already paired this polygon, skip it
    if is_paired(i, 1)
        continue
    end
    
    % Check that the polygon we are trying to match with says that we are
    % its pair, otherwise we are an unmatched polygon
    if (i ~= min_dist_idx(min_dist_idx(i)))
        unpaired_polygon(unpaired_polygon_counter+1).polygon = xml_struct(i).points;
        unpaired_polygon(unpaired_polygon_counter+1).centroid = polygon_centroids(i, :);
        unpaired_polygon(unpaired_polygon_counter+1).possible_pair = min_dist_idx(i);
        unpaired_polygon_counter = unpaired_polygon_counter+1;
        continue
    end
    
    % Check if the current polygon is the inner or outer polygon
    is_current_inner = isPointInPolygon(xml_struct(i).points(1,:), xml_struct(min_dist_idx(i)).points);
    if is_current_inner
        % If we are the inner polygon, record data in the matched polygon
        % array as such
        polygon_struct(completed_pairs_counter+1).inner_polygon = xml_struct(i).points;
        polygon_struct(completed_pairs_counter+1).outer_polygon = xml_struct(min_dist_idx(i)).points;
        polygon_struct(completed_pairs_counter+1).inner_centroid = polygon_centroids(i, :);
        polygon_struct(completed_pairs_counter+1).outer_centroid = polygon_centroids(min_dist_idx(i), :);
    else
        % If we are the outer polygon, record data in the matched polygon
        % array as such
        polygon_struct(completed_pairs_counter+1).inner_polygon = xml_struct(min_dist_idx(i)).points;
        polygon_struct(completed_pairs_counter+1).outer_polygon = xml_struct(i).points;
        polygon_struct(completed_pairs_counter+1).inner_centroid = polygon_centroids(min_dist_idx(i), :);
        polygon_struct(completed_pairs_counter+1).outer_centroid = polygon_centroids(i, :);
    end
    
    label_dist = ((h_text_position(:, 1) - polygon_struct(completed_pairs_counter+1).inner_centroid(1)).^2 + ...
        (h_text_position(:, 2) - polygon_struct(completed_pairs_counter+1).inner_centroid(2)).^2);
    
    % Also match the text label to the inner polygon
    text_label_idx = find(label_dist < 1E-14);
    
    % If we fail to match exactly, print and warning and go to the closest
    if isempty(text_label_idx)
        [dist_check, text_label_idx] = min(label_dist);
        fprintf(sprintf('Warning! Label %d failed to match exactly, matching closest instead.\\n', str2double(h_text(text_label_idx).String)))
        failed_match_label = [failed_match_label, str2double(h_text(text_label_idx).String)];
    end
    
    % Convert the text label to a number
    polygon_struct(completed_pairs_counter+1).text_label = h_text(text_label_idx).String;
    polygon_struct(completed_pairs_counter+1).num = str2double(h_text(text_label_idx).String);
    
    % Set that current polygon has been matched
    is_paired(i, 1) = true;
    is_paired(min_dist_idx(i), 1) = true;
    
    % Increase the number of polygons matched
    completed_pairs_counter = completed_pairs_counter+1;
end

if length([polygon_struct.num]) ~= length(unique([polygon_struct.num]))
    error('Error: A single label was used twice so we have a duplicate number')
end

% Sort the polygons by number
[~, sort_idx] = sort([polygon_struct.num]);
polygon_struct = polygon_struct(sort_idx);



% Create a summary figure and copy the input figure into the first subplot
if is_iterate_plot
    summary_fig = figure;
    cur_ax = subplot(2, 2, 1);
    original_ax = findobj(input_fig, 'Type', 'Axes');
    copied_ax = copyobj(original_ax, summary_fig);
    set(copied_ax, 'Position', get(cur_ax,'position'));
    delete(cur_ax);
    
    % Make a distance check figure
    distance_check_fig = figure;
end

% Compute the wall distances by dividing the outer polygon into pieces and
% finding the closest distance between those points and the inner polygon
% to do this we first get user input of which polygon(s) to plot / compute
total_tic = tic;

% Prompt for what we want to do initially
prompt_start = sprintf('Enter polygon number (up to %d) or compute (a)ll polygons and output. Blank for polygon 1, (Q)uit: ', size(polygon_struct, 1));
while true
    % Get user input
    user_input = input(prompt_start, 's');
    % Also try to convert the input into a numerical value
    [num_input, num_convert_status] = str2num(user_input); %#ok<ST2NM>
    % Check what the input is
    if strcmpi(user_input, 'A')
        % If it is compute all polygons, see that we want to compute all
        % polygons (polygon to plot = -1 indicates this)
        polygon_num_to_compute = -1;
        break
    elseif isempty(user_input)
        % If the input is blank, plot the first polygon
        polygon_num_to_compute = 1;
        break
    elseif num_convert_status
        % Otherwise if the input is a numerical value
        if (num_input > 0) && (num_input <= size(polygon_struct, 1))
            % Check if it is in the value range and plot it if it is
            polygon_num_to_compute = num_input;
            break
        else
            % Otherwise print that it is outside the range
            fprintf('Input exceeds range of 1 to %d.\n', size(polygon_struct, 1))
        end
    elseif strcmpi(user_input, 'Q')
        return
    else
        % If input is not any of the choices, print that it is invalid
        fprintf('Invalid input, please try again.\n')
    end
end

% Prompt for what to do next the in the loop, this is used at the end of
% the loop
prompt_loop = sprintf('(N)ext, (P)revious, Compute (A)ll, enter a polygon number up to %d, (Q)uit: ', size(polygon_struct, 1));
% While we are still plotting / getting use input
while true
    % Work out the polygon to plot / compute
    if polygon_num_to_compute == -1
        % If we are trying to plot polygon -1, that means compute all polygons
        j_to_compute = 1:size(polygon_struct, 1);
    else
        % Otherwise plot the polygon in question
        j_to_compute = polygon_num_to_compute;
    end
    
    % Compute the wall distances of polygon j
    for j = j_to_compute
        single_tic = tic;
        % Compute the areas
        polygon_struct(j).inner_area = abs(polygonArea(polygon_struct(j).inner_polygon));
        polygon_struct(j).outer_area = abs(polygonArea(polygon_struct(j).outer_polygon));
        
        % We need to clean up any self intersections / broken bits of the
        % outer polygon to ensure that the visibility / shadow is correct
        % (i.e. we don't have somewhere completely shadowed because of an
        % incorrect join)
        % Start with the original outer polygon
        cleaned_outer_polygon = polygon_struct(j).outer_polygon;
        % Remove the last repeated index if it is repeated
        if (sum((cleaned_outer_polygon(end, :) - cleaned_outer_polygon(1, :)).^2) < 1E-15)
            cleaned_outer_polygon(end, :) = [];
        end
        
        % We are going to go through the points of the polygon one by one
        % and check the angle and delete any large angles between adjacent
        % lines
        points_checked = 0;
        % While we still have points left to check
        while points_checked <= size(cleaned_outer_polygon, 1) + 1
            % Work out the index of this point and the points before and
            % after. Handle the wrap around and the end
            current_point1_index = mod(points_checked - 1, size(cleaned_outer_polygon, 1))+1;
            current_point2_index = mod(points_checked, size(cleaned_outer_polygon, 1))+1;
            current_point3_index = mod(points_checked + 1, size(cleaned_outer_polygon, 1))+1;
            % Work out the two edges these points make up and compute the
            % angle between them via dot product
            edge1 = cleaned_outer_polygon(current_point2_index, :) - cleaned_outer_polygon(current_point1_index, :);
            edge2 = cleaned_outer_polygon(current_point3_index, :) - cleaned_outer_polygon(current_point2_index, :);
            angle = acosd(dot(edge1,edge2)/(norm(edge1)*norm(edge2)));
            
            % If the angle between is is very large
            if angle > 120
                % Delete the point in questions
                cleaned_outer_polygon(current_point2_index, :) = [];
                % Move back one point to check the angle of the previous
                % line which would change due to our deletion
                points_checked = points_checked - 1;
            else
                % Otherwise, if the angle is small, check the next point
                points_checked = points_checked + 1;
            end
        end
        
        % Once we have deleted all points with large angles, repeat the
        % last index angle for a closed polygon
        cleaned_outer_polygon = cleaned_outer_polygon([1:end, 1], :);
        
        % Divide the outer polygon into parts which we will use to compute the
        % distance
        test_points_outer = zeros((size(cleaned_outer_polygon, 1)-1)*num_line_divisions, 2);
        for i = 1:(size(cleaned_outer_polygon, 1)-1)
            % Divide each line into num_line_divisions subdivisions
            test_points_outer_temp = ...
                [linspace(cleaned_outer_polygon(i,1), cleaned_outer_polygon(i+1,1), num_line_divisions).', ...
                linspace(cleaned_outer_polygon(i,2), cleaned_outer_polygon(i+1,2), num_line_divisions).'];
            test_points_outer(((i-1)*num_line_divisions+1):(i*num_line_divisions), :) = test_points_outer_temp;
        end
        
        % Compute the distance between each subdivision as this will not be
        % the same between lines (the different lines in the polygon will be
        % different length). This gives distance along the perimeter
        % This would be the x axis of our wall distance against perimeter plot
        test_points_outer_dist_x = [0; sqrt(sum(diff(test_points_outer).^2,2))];
        test_points_outer_dist_x = cumsum(test_points_outer_dist_x);
        
        % Compute the close distance between the subdivision and the inner
        % polygon. This is computing closest distance from points on the outer
        % polygon to the inner polygon
        %     [test_points_outer_dist_y, test_points_outer_dist_pos] = distancePointPolyline(test_points_outer, polygon_struct(j).inner_polygon);
        
        % However need to deal with the case that this closest point might pass
        % through the outer polygon. So loop through our test points and work
        % out the visible surface of the inner polygon from that point
        % We will work out the closest distances of the points one by one
        test_points_outer_dist_y = zeros(size(test_points_outer, 1), 1);
        test_points_outer_dist_pos = test_points_outer_dist_y;
        % Also some points may be completely shadowed, if that occurs,
        % exclude those points and remember them
        failed_test_points = [];
        failed_test_point_indices = [];
        % Work out the visibility and wall distance of each point in turn
        for point_i = 1:size(test_points_outer, 1)
            % Current point we are working on to compute wall distance
            cur_test_point_outer = test_points_outer(point_i, :);
            
            % Create a bunch of line segments connecting our current point
            % we want to compute the wall distance of to the possible
            % locations on the inner polygon which have the closest
            % distance. Make these into one big polyline for intersection
            % calculation using polyxpoly
            x_line_of_sight = [repmat(cur_test_point_outer(1), size(polygon_struct(j).inner_polygon, 1), 1), polygon_struct(j).inner_polygon(:, 1), nan(size(polygon_struct(j).inner_polygon, 1), 1)].';
            y_line_of_sight = [repmat(cur_test_point_outer(2), size(polygon_struct(j).inner_polygon, 1), 1), polygon_struct(j).inner_polygon(:, 2), nan(size(polygon_struct(j).inner_polygon, 1), 1)].';
            % Compute any intersections between our candidate distances and
            % the outer polygon. No need check inner polygon because if we
            % hit the inner polygon we would have the closest distance and
            % going through it won't give us a closer distance
            [intersection_x, intersection_y, intersection_indices] = polyxpoly(cleaned_outer_polygon(:, 1), cleaned_outer_polygon(:, 2), x_line_of_sight(:), y_line_of_sight(:));
            
            % Remove all the self intersections that are happening at our
            % current point on the outer polygon
            origin_intersection_indices = abs((intersection_x - cur_test_point_outer(1)).^2 + (intersection_y - cur_test_point_outer(2)).^2) > 1E-15;
            intersection_x_cleaned = intersection_x(origin_intersection_indices);
            intersection_y_cleaned = intersection_y(origin_intersection_indices);
            intersection_indices_cleaned = intersection_indices(origin_intersection_indices, :);
            
            % Then any remaining intersections are points on the inner
            % polygon that are not visible from the current position
            shadowed_point_indices = intersection_indices_cleaned(:, 2);
            
            % Therefore, remove those as possibilities of closest distance
            % for our current point we are compute
            inner_polygon_visible = polygon_struct(j).inner_polygon;
            inner_polygon_visible(ceil(unique(shadowed_point_indices)/3), :) = NaN;
            
            % Then compute the closest distance for the remaining visible
            % points
            [cur_test_points_outer_dist_y, cur_test_points_outer_dist_pos] = distancePointPolyline(cur_test_point_outer, inner_polygon_visible);
            test_points_outer_dist_y(point_i) = cur_test_points_outer_dist_y;
            test_points_outer_dist_pos(point_i) = cur_test_points_outer_dist_pos;
            
            % Finally, if we find any points that completely shadowed,
            % remember those points so that we can exclude them and print
            % a warning
            if isnan(cur_test_points_outer_dist_y)
                warning('Failed to find wall distance for a point on the outer wall')
                failed_test_point_indices = [failed_test_point_indices; point_i];
                failed_test_points = [failed_test_points; cur_test_point_outer];
            end
        end
        
        % Delete any points that were competely shadowed for our wall
        % distance function, this will have the effect of the distance at
        % that location being interpolated from the wall distances before
        % and after
        test_points_outer_dist_y(failed_test_point_indices) = [];
        test_points_outer_dist_x(failed_test_point_indices) = [];
        test_points_outer_dist_pos(failed_test_point_indices) = [];
        test_points_outer(failed_test_point_indices, :) = [];
        
        
        
        % Next, we will interploate over the wall distances against perimeter
        % function so that we have equal spacing in the perimeter. This will
        % let us treat the wall distances as discrete and draw boxplots etc
        
        % First, remove duplicate points from our perimeters and wall distances
        temp_cleaned = [test_points_outer_dist_x, test_points_outer_dist_y];
        temp_cleaned = unique(temp_cleaned,'rows');
        test_points_outer_dist_x_cleaned = temp_cleaned(:, 1);
        test_points_outer_dist_y_cleaned = temp_cleaned(:, 2);
        
        % Then interpolate to give us even distances
        test_points_outer_dist_x_even = linspace(test_points_outer_dist_x(1), test_points_outer_dist_x(end), n_even);
        test_points_outer_dist_y_even = interp1(test_points_outer_dist_x_cleaned, test_points_outer_dist_y_cleaned, test_points_outer_dist_x_even);

        % Work out the coverage percentage
        covered_locations = test_points_outer_dist_y >= coverage_tolerance;
        crossing_locations = diff(covered_locations);
        crossing_indices = find(crossing_locations);
        perimeter_expanded = test_points_outer_dist_x;
        thickness_expanded = test_points_outer_dist_y;
        covered_locations_expanded = covered_locations;
        for i = flip(crossing_indices(:).')
            crossing_perimeter = interp1(test_points_outer_dist_y(i:i+1), test_points_outer_dist_x(i:i+1), coverage_tolerance);
            perimeter_expanded = [perimeter_expanded(1:i); crossing_perimeter; perimeter_expanded(i+1:end)];
            thickness_expanded = [thickness_expanded(1:i); coverage_tolerance; thickness_expanded(i+1:end)];
            covered_locations_expanded = [covered_locations_expanded(1:i); true; covered_locations_expanded(i+1:end)];
        end
        
        % figure
        % plot(test_points_outer_dist_x, test_points_outer_dist_y, '-x', perimeter_expanded, thickness_expanded, '+', [0, test_points_outer_dist_x(end)], [coverage_tolerance, coverage_tolerance])

        perimeter_distances = diff(perimeter_expanded);
        perimeter_distances_to_count = thickness_expanded(1:end-1) >= coverage_tolerance & thickness_expanded(2:end) >= coverage_tolerance;
        coverage_distance = sum(perimeter_distances(perimeter_distances_to_count));
        
        % Compute the mean, sd and median
        polygon_struct(j).mean = mean(test_points_outer_dist_y_even);
        polygon_struct(j).sd = std(test_points_outer_dist_y_even);
        polygon_struct(j).median = median(test_points_outer_dist_y_even);

        % Coverage percentage
        polygon_struct(j).coverage_percentage = coverage_distance/perimeter_expanded(end)*100;
        
        % Compute the quantiles
        polygon_struct(j).quantiles = quantile(test_points_outer_dist_y_even, [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]);
        polygon_struct(j).min = min(test_points_outer_dist_y_even);
        polygon_struct(j).q_5 = polygon_struct(j).quantiles(1);
        polygon_struct(j).q_10 = polygon_struct(j).quantiles(2);
        polygon_struct(j).q_25 = polygon_struct(j).quantiles(3);
        polygon_struct(j).q_50 = polygon_struct(j).quantiles(4);
        polygon_struct(j).q_75 = polygon_struct(j).quantiles(5);
        polygon_struct(j).q_90 = polygon_struct(j).quantiles(6);
        polygon_struct(j).q_95 = polygon_struct(j).quantiles(7);
        polygon_struct(j).max = max(test_points_outer_dist_y_even);
        
        
        % Plot results one by one
        if is_iterate_plot
            
            % Move the axes of the copied input plot to the right position
            figure(summary_fig)
            cur_ax = subplot(2, 2, 1);
            poly_x_min = min(polygon_struct(j).outer_polygon(:, 1))*resolution;
            poly_x_max = max(polygon_struct(j).outer_polygon(:, 1))*resolution;
            poly_y_min = min(polygon_struct(j).outer_polygon(:, 2))*resolution;
            poly_y_max = max(polygon_struct(j).outer_polygon(:, 2))*resolution;
            poly_x_diff = poly_x_max - poly_x_min;
            poly_y_diff = poly_y_max - poly_y_min;
            ax_limits = [poly_x_min - 0.1*poly_x_diff, poly_x_max + 0.1*poly_x_diff, ...
                poly_y_min - 0.1*poly_y_diff, poly_y_max + 0.1*poly_y_diff];
            axis(cur_ax, ax_limits)
            
            % Plot the wall distance against perimeter plot
            cur_ax = subplot(2, 2, 2);
            plot(cur_ax, test_points_outer_dist_x, test_points_outer_dist_y, [0, test_points_outer_dist_x(end)], [coverage_tolerance, coverage_tolerance])
            xlabel('Perimeter (micron)')
            ylabel('Wall distance (micron)')
            sgtitle(sprintf("Polygon %d", polygon_struct(j).num))
            
            % Draw a boxplot of the wall distances
            cur_ax = subplot(2, 2, 3);
            boxplot(cur_ax, test_points_outer_dist_y_even, ...
                repmat(sprintf("P%d", polygon_struct(j).num), ...
                size(test_points_outer_dist_y_even, 2), 1), 'Orientation', 'horizontal')
            xlabel('Wall distance (micron)')
            
            % Draw a histogram of wall distances
            cur_ax = subplot(2, 2, 4);
            cur_hist = histogram(cur_ax, test_points_outer_dist_y_even);
            xlabel('Wall distance (micron)')
            
            savefig(summary_fig, [temp_folder, filename, sprintf(' - Polygon %d', polygon_struct(j).num), ' - Statistics Summary.fig'])
            
            % Create a distance check plot. Draw lines showing all the wall
            % distances we found to check the distances make sense.
            figure(distance_check_fig)
            % Inner polygon
            plot(polygon_struct(j).inner_polygon(:,1), polygon_struct(j).inner_polygon(:,2), 'k');
            cur_ax = gca;
            cur_ax.YDir = copied_ax.YDir;
            axis equal
            hold on
            % Cleaned outer polygon
            plot(cleaned_outer_polygon(:,1), cleaned_outer_polygon(:,2), 'r')
            % Original otuer polygon
            plot(polygon_struct(j).outer_polygon(:,1), polygon_struct(j).outer_polygon(:,2), 'k')
            % Lines showing the actual wall distances
            test_points_outer_dist_points = polylinePoint(polygon_struct(j).inner_polygon, test_points_outer_dist_pos);
            plot([test_points_outer_dist_points(:, 1), test_points_outer(:,1)].', [test_points_outer_dist_points(:, 2), test_points_outer(:,2)].')
            % Any excluded points
            if ~isempty(failed_test_points)
                plot(failed_test_points(:, 1), failed_test_points(:, 2), '*m')
            end
            hold off
            title(sprintf("Polygon %d", polygon_struct(j).num))
            xlabel('x (micron)')
            ylabel('y micron)')
            savefig(distance_check_fig, [temp_folder, filename, sprintf(' - Polygon %d', polygon_struct(j).num), ' - Wall Distance Check.fig'])
            
            % Convert the current results to a table and display them
            cur_result = struct2table(polygon_struct(j), 'AsArray', true);
            cur_result = cur_result(:, [6, 7:12, 14:22])
            
            % Wait for a button press before plotting the next one
            %         waitforbuttonpress
            drawnow
        end
        toc(single_tic)
    end
    
    % If we have just computed all polygon, then end our input loop
    if polygon_num_to_compute == -1
        break
    end
    
    % Otherwise get user input of what to do next
    while true
        % Get input
        user_input = input(prompt_loop, 's');
        % Try to convert input to numeric value
        [num_input, num_convert_status] = str2num(user_input);
        if strcmpi(user_input, 'N')
            % If input is next, try to compute next polygon
            if polygon_num_to_compute == size(polygon_struct, 1)
                % If we are at last polygon, print and get another user
                % input
                fprintf('Reached last polygon, no more to go next to.\n')
            else
                % Otherwise, set that we are polygon next polygon and
                % stop getting input
                polygon_num_to_compute = polygon_num_to_compute+1;
                break
            end
        elseif strcmpi(user_input, 'P')
            % Compute previous polygon
            if polygon_num_to_compute == 1
                % If we are at first polygon, print this and get another
                % input
                fprintf('At first polygon, no more to go previous to.\n')
            else
                % Otherwise, set previous polygon to be computed then
                % stop getting input
                polygon_num_to_compute = polygon_num_to_compute-1;
                break
            end
        elseif strcmpi(user_input, 'A')
            % If we are asked to compute all, set -1 to indiciate this fact
            polygon_num_to_compute = -1;
            break
        elseif num_convert_status
            % If inputted numerical value, check it is in valid range
            if (num_input > 0) && (num_input <= size(polygon_struct, 1))
                % If it is, set we are computing that and stop getting
                % input
                polygon_num_to_compute = num_input;
                break
            else
                % Otherwise print it is outside range and get another input
                fprintf('Input exceeds range of 1 to %d.\n', size(polygon_struct, 1))
            end
        elseif strcmpi(user_input, 'Q')
            return
        else
            % If input invalid, print that and get another input
            fprintf('Invalid input, please try again.\n')
        end
    end
end
toc(total_tic)



% Convert to table
polygon_table = struct2table(polygon_struct);
% Write out a table omitting the polygons and centroids and ellipses but
% including a first row with index numbers
if exist([filename, ' - Wall Distances.xlsx'], 'file')
    delete([filename, ' - Wall Distances.xlsx']);
end
writetable(polygon_table(:, [6, 7:12, 14:22]), [filename, ' - Wall Distances.xlsx'])


for i = 1:length(failed_match_label)
    fprintf(sprintf('Warning! Label %d failed to match exactly, matching closest instead.\\n', failed_match_label(i)))
end








% After we have calculated everything, get input again to check the figures
prompt_check = sprintf('(N)ext, (P)revious, enter a polygon number up to %d, (Q)uit: ', size(polygon_struct, 1));
polygon_num_to_compute = j;
while true
    % Get input
    user_input = input(prompt_check, 's');
    % Try to convert input to numeric value
    [num_input, num_convert_status] = str2num(user_input);
    if strcmpi(user_input, 'N')
        % If input is next, try to compute next polygon
        if polygon_num_to_compute == size(polygon_struct, 1)
            % If we are at last polygon, print and get another user
            % input
            fprintf('Reached last polygon, no more to go next to.\n')
            continue
        else
            % Otherwise, set that we are polygon next polygon and
            % stop getting input
            polygon_num_to_compute = polygon_num_to_compute+1;
        end
    elseif strcmpi(user_input, 'P')
        % Compute previous polygon
        if polygon_num_to_compute == 1
            % If we are at first polygon, print this and get another
            % input
            fprintf('At first polygon, no more to go previous to.\n')
            continue
        else
            % Otherwise, set previous polygon to be computed then
            % stop getting input
            polygon_num_to_compute = polygon_num_to_compute-1;
        end
    elseif num_convert_status
        % If inputted numerical value, check it is in valid range
        if (num_input > 0) && (num_input <= size(polygon_struct, 1))
            % If it is, set we are computing that and stop getting
            % input
            polygon_num_to_compute = num_input;
        else
            % Otherwise print it is outside range and get another input
            fprintf('Input exceeds range of 1 to %d.\n', size(polygon_struct, 1))
            continue
        end
    elseif strcmpi(user_input, 'Q')
        break
    else
        % If input invalid, print that and get another input
        fprintf('Invalid input, please try again.\n')
        continue
    end
    
    % Reopen the summary window, open it hidden and the copy it to the
    % correct figure
    fig_read_temp = openfig([temp_folder, filename, sprintf(' - Polygon %d', polygon_num_to_compute), ' - Statistics Summary.fig'], 'invisible');
    clf(summary_fig)
    copyobj(fig_read_temp.Children, summary_fig);
    close(fig_read_temp);
    
    % Reopen the wall distance window, open it hidden and the copy it to 
    % the correct figure
    fig_read_temp = openfig([temp_folder, filename, sprintf(' - Polygon %d', polygon_num_to_compute), ' - Wall Distance Check.fig'], 'invisible');
    clf(distance_check_fig)
    copyobj(fig_read_temp.Children, distance_check_fig);
    close(fig_read_temp);
    
    drawnow
end

if is_iterate_plot
    close(input_fig, summary_fig, distance_check_fig)
else
    close(input_fig)
end

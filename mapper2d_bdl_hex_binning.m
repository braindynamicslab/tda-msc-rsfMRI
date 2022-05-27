%WGGGG
%
% aim - to code mapper 2d in matlab
% author - saggar@stanford.edu (7.11.2019)
% 
%
% output: adja, num_vertices, level_of_vertex, pts_in_vertex,
% points_in_level, and vertices_in_level
%
% inspired by TDAMapper R 
% and Gurjeet Singh's original paper - Singh et al. 2007
% 
% Version History
% - [7.11.19] Originally wrote it to check across 4 directions in the binning space
% for finding overlap across bins. Thus, constrained the nodes to have a max
% degree 4. [TDAMapper R style]
% - [7.15.19] Modified to include the other 4 diagonal directions for
% finding overlappping bins. Thus, constrained the nodes to a max degree of 8. 
% - [7.16.19] Todo: After talking with Samir C. we will now try to expand
% the binning strategy from 2D bins to several other styles.
% - [7.17.19] Update [non-metric version] - instead of looking over the neighbors for overlap, we will 
% instead just connect levels irrespective of neighorhood status - as long
% as they share a data point
% - [7.22.19] Updated [non-metric-version] - added prunning of
% nodes/vertices - thus vertices with same set of data points are merged
% into one vertex.
% - [11.25.19] Adding hex-binning to the Mapper code
%
function [adja, adja_pruned, pts_in_vertex, pts_in_vertex_pruned] = mapper2d_bdl_hex_binning(distMat, filter_values, num_intervals, percent_overlap, num_bins_clustering, nsides)

 %initialize variables
 vertex_index = 0;
 if ~exist('nsides','var')
    nsides = 6; %for hex, and 4 for regular.
 end

 % indexed from 1 to the number of vertices
 level_of_vertex = [];
 pts_in_vertex = [];
 
 % indexed from 1 to the number of levels
 pts_in_level = {};
 vts_in_level = {}; 
 
 filter_min_1 = floor(min(filter_values(:,1)));
 filter_max_1 = ceil(max(filter_values(:,1)));
 filter_min_2 = floor(min(filter_values(:,2)));
 filter_max_2 = ceil(max(filter_values(:,2)));

 interval_length_1 = (filter_max_1 - filter_min_1) / (num_intervals(1) - (num_intervals(1) - 1) * percent_overlap/100 );
 interval_length_2 = (filter_max_2 - filter_min_2) / (num_intervals(2) - (num_intervals(2) - 1) * percent_overlap/100 );

 step_size_1 = interval_length_1 * (1 - percent_overlap/100);
 step_size_2 = interval_length_2 * (1 - percent_overlap/100);

 num_levels = num_intervals(1) * num_intervals(2);

 level_indices_1 = repmat(1:num_intervals(1), [1,num_intervals(2)]);
 level_indices_2 = reshape(repmat(1:num_intervals(2), num_intervals(1),1),[num_levels,1])';
 %figure; 
 for level = 1:1:num_levels
    level_1 = level_indices_1(level); 
    level_2 = level_indices_2(level);
    
    min_val_level_1 = filter_min_1 + (level_1-1)*step_size_1;
    min_val_level_2 = filter_min_2 + (level_2-1)*step_size_2;
    
    max_val_level_1 = min_val_level_1 + interval_length_1;
    max_val_level_2 = min_val_level_2 + interval_length_2;
    
    % step 1: find center of min & max and create hexagon ther
    % step 2: find pts that are within that hexagon
    % step 3: we will start intersecting
    center_level_1 = (min_val_level_1+max_val_level_1)./2;
    center_level_2 = (min_val_level_2+max_val_level_2)./2;
    hex_poly = nsidedpoly(nsides, 'Center', [center_level_1, center_level_2], 'SideLength',  (interval_length_1+interval_length_2)./2);  
    pts_level = isinterior(hex_poly, filter_values);
    %plot(hex_poly); hold on;
    %pts_level = (filter_values(:,1) >= min_val_level_1) & (filter_values(:,2) >= min_val_level_2) & (filter_values(:,1) <= max_val_level_1) & (filter_values(:,2) <= max_val_level_2);
    
    num_pts_level = sum(pts_level);
    pts_in_level{level} = find(pts_level);
    
    if num_pts_level == 0
       %fprintf(1,'Level set is empty\n');
       vts_in_level{level} = -1;
       continue;        
    elseif num_pts_level == 1
       %fprintf(1,'Level set has only 1 pt\n');
       num_vts_level = 1;
       cluster_indices_within_level = [1];
       vts_in_level{level} = vertex_index + (1:num_vts_level);
    elseif num_pts_level > 1
       %fprintf(1,'Level set has only %d pts\n', num_pts_level);
       level_distMat = distMat(pts_level, pts_level);
       level_max_dist = max(level_distMat(:));
       [Z, cutoff] = find_cluster_cutoff(level_distMat, num_bins_clustering);
       if cutoff <0
           fprintf(2,'cutoff <0\n');
           cutoff = 1;
           
       end
       cluster_indices_within_level = cluster(Z, 'cutoff',cutoff, 'criterion','distance');
       num_vts_level = max(cluster_indices_within_level);
       vts_in_level{level} = vertex_index + (1:num_vts_level);
    end
    
    for j = 1:1:num_vts_level
       vertex_index = vertex_index + 1;
       level_of_vertex(vertex_index) = level;
       pts_in_vertex{vertex_index} = pts_in_level{level}(cluster_indices_within_level==j);
        
    end
     
 end
 %scatter(filter_values(:,1), filter_values(:,2), 15,'k');
 % pruning nodes/vertices that have same number of points in them
 pts_in_vertex_pruned = pts_in_vertex;
 for i = 1:1:vertex_index
     for j = i:1:vertex_index
         if i == j
             continue;
         end
         k1 = pts_in_vertex{i};
         k2 = pts_in_vertex{j};
         if isequal(sort(k1),sort(k2))
             %fprintf(1,'found one vertex for pruning i=%d, j=%d\n', i,j);
             pts_in_vertex_pruned{j} = [];
         end
     end
 end
 tmp = {};
 k = 1;
 for i = 1:1:vertex_index
    if ~isempty(pts_in_vertex_pruned{i})
        tmp{k} = pts_in_vertex_pruned{i};
        k = k+1;
    end
 end
 pts_in_vertex_pruned = tmp;
 fprintf(1,'Pruning done\n');
 
 vertex_index_pruned = length(pts_in_vertex_pruned);
 adja_pruned =  zeros(vertex_index_pruned, vertex_index_pruned);
 for i = 1:1:vertex_index_pruned
     for j = i:1:vertex_index_pruned
         if i == j
             continue;
         end
         k1 = pts_in_vertex_pruned{i};
         k2 = pts_in_vertex_pruned{j};
         adja_pruned(i,j) = (length(intersect(k1,k2))>0);
         adja_pruned(j,i) = adja_pruned(i,j);
     end
 end
 fprintf(1,'Prunned Mapper Done\n');
 
 
 adja = zeros(vertex_index, vertex_index);
%  for i = 1:1:vertex_index
%      for j = i:1:vertex_index
%          if i == j
%              continue;
%          end
%          k1 = pts_in_vertex{i};
%          k2 = pts_in_vertex{j};
%          adja(i,j) = (length(intersect(k1,k2))>0);
%          adja(j,i) = adja(i,j);
%      end
%  end
%  fprintf(1,'Mapper Done\n');
 
end

function [Z, cutoff] = find_cluster_cutoff(distMat, num_bins_clustering)
   Z = linkage(distMat(tril(true(length(distMat)),-1))');
   if length(Z(:,3)) == 1 % only two points hence one cluster
       fprintf(1,'Only two points found for clustering\n');
       cutoff = Inf;
       return;
   end
   
   lens = [Z(:,3)' max(max(distMat))];
   [numBins, bc] = hist(lens, num_bins_clustering);
   z = find(numBins==0);
   if (sum(z) == 0)
       cutoff = Inf;
       return;
   else
       cutoff = bc(z(1)); % pickup the smallest index of z
   end


end

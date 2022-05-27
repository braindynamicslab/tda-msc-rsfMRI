% WGGGG
%
% aim - to code neighborhood based filter function in Matlab
% author - saggar@stanford.edu (7.11.2019)
% 
%
% output: knnGraph (as a table, binary and weighted versions)
%
% inspired by the paper on penalized KNN graph in Baya & Granitto 2011 BMC
% Bioinformatics
% 
% Version History
% - [7.11.19] Wrote it to perform knn graph construction with outlier
% removal and penalized edges to construct a single connected component
%
%
function [knnGraphTbl, knnGraph_dense_bin, knnGraph_dense_wtd, knnGraph_dense_bin_conn, knnGraph_dense_wtd_conn]= createPKNNG_bdl(distMat, num_k)
   % create neighborhood graph
   knnGraphTbl = zeros(size(distMat,1), num_k);
   knnGraph_bin = zeros(size(distMat));
   knnGraph_wtd = zeros(size(distMat));
   for n = 1:1:size(distMat,1)
      [tmp idx] = sort(distMat(n,:),'ascend');
      knnGraphTbl(n,:) = idx(2:2+num_k-1);
      for l = 1:1:num_k
         knnGraph_bin(n, knnGraphTbl(n,l)) = 1; 
         knnGraph_wtd(n, knnGraphTbl(n,l)) = distMat(n,knnGraphTbl(n,l)); 
          
      end
   end
   
   % remove outliers in knn graph to find densely connected subgraphs
   % first remove non-reciprocal connections
   knnGraph_dense_bin = knnGraph_bin;
   knnGraph_dense_wtd = knnGraph_wtd;
   for n = 1:1:size(distMat,1)
      tmp = knnGraphTbl(n,:);
      for l = 1:1:num_k         
         if ~ismember(n,knnGraphTbl(tmp(l),:))
             knnGraph_dense_bin(n,tmp(l)) = 0;
             knnGraph_dense_bin(tmp(l),n) = 0;

             knnGraph_dense_wtd(n,tmp(l)) = 0;
             knnGraph_dense_wtd(tmp(l),n) = 0;
             
         end
         
      end
       
   end
   
   % connect disconnected components of the graph
   g = graph(knnGraph_dense_wtd);
   knnGraph_dense_wtd_conn = knnGraph_dense_wtd;
   knnGraph_dense_bin_conn = knnGraph_dense_bin;
   bins = conncomp(g);
   nComp = max(bins);
   if nComp > 1 
      for c = 1:1:nComp          
          for d = c+1:1:nComp
              nodes_c = find(bins==c);
              nodes_d = find(bins==d);
              max_edge_c = max(max(knnGraph_dense_wtd(nodes_c, nodes_c)));
              max_edge_d = max(max(knnGraph_dense_wtd(nodes_d, nodes_d)));
              if length(nodes_d) == 1
                 [val, best_nodes_c] = min(distMat(nodes_d, nodes_c));
                 knnGraph_dense_wtd_conn(nodes_d, best_nodes_c) = val * exp(val/max_edge_c);
                 knnGraph_dense_wtd_conn(best_nodes_c, nodes_d) = knnGraph_dense_wtd_conn(nodes_d, best_nodes_c);
                 knnGraph_dense_bin_conn(nodes_d, best_nodes_c) = 1;
                 knnGraph_dense_bin_conn(best_nodes_c, nodes_d) = 1;
              elseif length(nodes_c) == 1
                 [val, best_nodes_d] = min(distMat(nodes_c, nodes_d));
                 knnGraph_dense_wtd_conn(nodes_c, best_nodes_d) = val * exp(val/max_edge_d);
                 knnGraph_dense_wtd_conn(best_nodes_d, nodes_c) = knnGraph_dense_wtd_conn(nodes_c, best_nodes_d);
                 knnGraph_dense_bin_conn(nodes_c, best_nodes_d) = 1;
                 knnGraph_dense_bin_conn(best_nodes_d, nodes_c) = 1;                  
              else % find the best pair of nodes between nodes_c and nodes_d that can be connected
                 tmp = distMat(nodes_c, nodes_d);
                 val = min(tmp(:));
                 [nodes_c_min, nodes_d_min] = find(tmp==val);
                 nodes_c_min = nodes_c(nodes_c_min);
                 nodes_d_min = nodes_d(nodes_d_min);
                 knnGraph_dense_wtd_conn(nodes_c_min, nodes_d_min) = val * exp(val/max(max_edge_c,max_edge_d));
                 knnGraph_dense_wtd_conn(nodes_d_min, nodes_c_min) = knnGraph_dense_wtd_conn(nodes_c_min, nodes_d_min);
                 knnGraph_dense_bin_conn(nodes_c_min, nodes_d_min) = 1;
                 knnGraph_dense_bin_conn(nodes_d_min, nodes_c_min) = 1;
              end
                  
          end
          
      end
   end
   
   
end
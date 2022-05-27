%WGGGG
%
%
% Aims - 
%
% [1] analyzes Mapper output and create figures from the paper below
%
% Please cite:
% Saggar, M., Shine, J.M., Liegeois, R., Dosenbach, N.U.F., Fair, D. 2021. 
% Precision dynamical mapping using topological data analysis reveals a 
% unique hub-like transition state at rest. BioRxiv
%
% date - 5.31.2021
% author - saggar@stanford.edu
%

%% generating figure 2, evidence of hub nodes in the real as compared to null data
clear; close all; clc;

data_folder = 'output';
session = 'odd';
metric = 'euclidean';
numNull = 25; 
doi = [21, 36]; % degree of interest (hump observed in the degree dist plot)
p=[];
prop_doi = [];
prop_doi_ar = [];
prop_doi_pr = [];  
N_re = [];
N_ar = [];
N_pr = [];
output_name = 'mapperout';
for null_n = 1:1:numNull
    null_n
    for s = 1:1:10
        mapperout = load(sprintf('%s/sub-MSC%02d/sub-MSC%02d_%s_runs_mat_metric_%s_%s.mat', data_folder, s, s, session, metric, output_name));           
        try
            mapperout_ar = load(sprintf('%s/sub-MSC%02d/sub-MSC%02d_AR_null%03d_%s_runs_mat_metric_%s_%s.mat', data_folder, s, s, null_n, session, metric, output_name));
        catch
            continue;
        end
        try 
            mapperout_pr = load(sprintf('%s/sub-MSC%02d/sub-MSC%02d_PR_null%03d_%s_runs_mat_metric_%s_%s.mat', data_folder, s, s, null_n, session, metric, output_name));
        catch
            continue;
        end

        % calculating degree-distribution for the mapper generated graph        
        deg_ar = degrees_und(mapperout_ar.nodeBynode);
        deg_pr = degrees_und(mapperout_pr.nodeBynode);                
        deg_re = degrees_und(mapperout.nodeBynode);
                
        [N_re(null_n, s,:), ~] = histcounts(deg_re, 'Normalization', 'probability', 'BinLimits', [0 50], 'BinWidth',1);
        [N_ar(null_n, s,:), ~] = histcounts(deg_ar, 'Normalization', 'probability', 'BinLimits', [0 50], 'BinWidth',1);
        [N_pr(null_n, s,:), ~] = histcounts(deg_pr, 'Normalization', 'probability', 'BinLimits', [0 50], 'BinWidth',1);        
        
        prop_doi(null_n, s) = sum((deg_re>doi(1)))./length(deg_re);
        prop_doi_ar(null_n, s) = sum((deg_ar>doi(1)))./length(deg_ar);
        prop_doi_pr(null_n, s) = sum((deg_pr>doi(1)))./length(deg_pr);       
    end
end

figure; 
plot(squeeze(mean(N_re,1)));hold on
plot(squeeze(mean(N_ar,1)));
plot(squeeze(mean(N_pr,1)));
legend('Real','AR','PR');

%% generating figure 3 & 5, annotating graphs using resting state networks and topographic gradient
% In the current version we create metaInfo data that can be visualized outside
% Matlab using our python implementation at https://braindynamicslab.github.io/dyneusr/
% Matlab version for similar viz. is coming soon...

clear; close all; clc;
cmap = [
254	4	2
254	206	153
255	251	0
30	253	4
253	176	245
68	185	217
0	0	0
91	22	137
18	255	235
184	127	191
255	255	255]./255;
nclust = 10;
num_k = 30; 
res_val = 30; 
gain_val = 70; 

data_folder = 'output';
session = 'odd';
metric = 'euclidean';
net_names_allelse = {'DMN','VIS','FRP','DAN','MOT','VAN','SAL','COP','SM','AUD'};
net_names_s4 =      {'DMN','VIS','FRP','DAN','MOT',            'COP','SM','AUD'}; % subj4 has van and sal unavailable
net_names_s10 =     {'DMN','VIS','FRP','DAN','MOT','VAN',      'COP','SM','AUD'}; % subj10 has sal unavailable

json_thrval_pos = 0.5;
json_saveOrNot = 1;
winnerTakeAll = 0; 

doi = [21, 36]; % degree of interest (hump observed in the degree dist plot)

corrProp_sub = [];

for s = 1%:1:10
    sbj_name = sprintf('sub-MSC%02d-Session-%s', s, session)
    mapperout = load(sprintf('%s/sub-MSC%02d/sub-MSC%02d_%s_runs_mat_metric_%s_mapperout_may31_2021.mat', data_folder, s, s, session, metric));           

    % degree annotation 
    deg = degrees_und(mapperout.nodeBynode);    

    % doi_nodes
    doi_nodes = (deg>doi(1));

    % plotting graphs using d3.js    
    if s==4
        net_names = net_names_s4;
        % pushing the metaInfo as networks 7 and 8 are missing in Sub04
        metaInfoNetData = zeros(size(mapperout.networkDataVertex,1),10);
        metaInfoNetData(:,1:5) = mapperout.networkDataVertex(:,1:5);
        metaInfoNetData(:,8:10) = mapperout.networkDataVertex(:,6:8);        
    elseif s==10
        net_names = net_names_s10;
        % pushing the metaInfo as network 8 is missing in Sub10
        metaInfoNetData = zeros(size(mapperout.networkDataVertex,1),10);
        metaInfoNetData(:,1:6) = mapperout.networkDataVertex(:,1:6);
        metaInfoNetData(:,8:10) = mapperout.networkDataVertex(:,7:9);        

    else
        net_names = net_names_allelse;
        metaInfoNetData = mapperout.networkDataVertex;
    end

    [entropy1] = calcEntropy(mapperout.nodeTpMat, mapperout.nodeBynode, metaInfoNetData);
    [hubness95, hubness99] = computeHubness(mapperout.nodeBynode, doi, deg);       
    [dom, corrProp] = calcDominance(mapperout.nodeTpMat, mapperout.nodeBynode, metaInfoNetData, json_thrval_pos);
 
    % write json for proportional rsn based piecharts
    metaInfo.value = metaInfoNetData;
    metaInfo.deg = 100*zscoreScaling(deg);
    metaInfo.ses = mapperout.sess_id;
    metaInfo.ent = 100*zscoreScaling(entropy1);
    metaInfo.doi = doi_nodes+1; %adding 1 for plotting purposes
    metaInfo.hubs99 = hubness99 + 1; %adding 1 for plotting purposes
    metaInfo.dom = dom; 
    
    
    corrProp_sub(s,:) = corrProp(:);
    metaInfo_s(s) = metaInfo;        
end

%% attempting to D3.js like force layouts in Matlab. 
% please note - we didn't generate pie-charts in Matlab, hence the hub nodes
% don't look uniformly distributed.
% RSN based notation
figure; plot(graph(mapperout.nodeBynode),'Layout', 'force','UseGravity','on','MarkerSize',2+5*log10(sum(mapperout.nodeTpMat,2)),'nodecdata',metaInfo.dom);
colormap(cmap);
% SD based notation - Topographic gradient
figure; plot(graph(mapperout.nodeBynode),'Layout', 'force','UseGravity','on','MarkerSize',2+5*log10(sum(mapperout.nodeTpMat,2)),'nodecdata',metaInfo.ent);
caxis([0 60])
%% generating figure 4: subject specificity
clear; close all; clc;
set(groot, 'DefaultTextInterpreter', 'none')

sub_spec = zeros(22,121);
load('output/metaInfo_s_msc_odd.mat')
sub_spec(1:2:20,:) = corrProp_sub;
sub_spec(21,:) = mean(corrProp_sub);
metaInfo_s_odd = metaInfo_s;

load('output/metaInfo_s_msc_even.mat')
sub_spec(2:2:20,:) = corrProp_sub;
sub_spec(22,:) = mean(corrProp_sub);
metaInfo_s_even = metaInfo_s;     

[r,p]=corrplot(sub_spec');
imagesc(r);

%% generating figure 6: traversal in time domain
clear; close all; clc;
set(groot, 'DefaultTextInterpreter', 'none')

net_names_allelse = {'DMN','VIS','FRP','DAN','MOT','VAN','SAL','COP','SM','AUD'};
net_names_s4 =      {'DMN','VIS','FRP','DAN','MOT',            'COP','SM','AUD'};
net_names_s10 =     {'DMN','VIS','FRP','DAN','MOT','VAN',      'COP','SM','AUD'};

data_folder = 'output';
session = 'odd';
metric = 'euclidean';
json_thrval_pos = 0.5;
load(sprintf('output/metaInfo_s_msc_%s',session),'metaInfo_s','corrProp_sub')

doi = [21, 36]; % degree of interest
cmap = [
254	4	2
254	206	153
255	251	0
30	253	4
253	176	245
68	185	217
0	0	0
91	22	137
18	255	235
184	127	191
255	255	255]./255;

ciu_bin_s = {};
ciu_s = {};
tr_color_s = {};

for s = 1:1:10
    sbj_name = sprintf('sub-MSC%02d-Session-%s', s, session)
    mapperout = load(sprintf('%s/sub-MSC%02d/sub-MSC%02d_%s_runs_mat_metric_%s_mapperout_may31_2021.mat', data_folder, s, s, session, metric));           
 
    % getting fd frames for proper Markov chain estimation
    fd_files = dir(sprintf('data/sub-MSC%02d/*TMASK.txt',s));   
    
    if strcmp(session,'odd') == 1
        ses_idx = 1:2:10
    elseif strcmp(session,'even') == 1
        ses_idx = 2:2:10
    end
    
    tr_label = [];
    sess_label = [];
    k_itr = 1;
    fd_across_ses = [];
    for ses = ses_idx
        fd = load([fd_files(ses).folder '/' fd_files(ses).name]);
        
        fd_across_ses = [fd_across_ses; fd];
        
       % fd(end) = -1;            
        if fd(end) == 1 %last frame of a session is good, then make it 0 so that we can tag end of session
            fd(end) = -1;            
        end
        tr_label = [tr_label; fd];
        sess_label = [sess_label; ones(length(fd),1)*k_itr];
        k_itr = k_itr + 1;
    end
    
    % plotting graphs using d3.js    
    if s==4
        net_names = net_names_s4;
    elseif s==10
        net_names = net_names_s10;        
    else
        net_names = net_names_allelse;
    end    
    
    nodes_hub = find(metaInfo_s(s).hubs99==2); % 2 hubs, 1 non hubs
    nodes_nothub = find(metaInfo_s(s).hubs99==1);
    dom = metaInfo_s(s).dom;
    
    hub_trs_wtd = sum(mapperout.nodeTpMat(nodes_hub,:))./length(nodes_hub);    
    nothub_trs_wtd = sum(mapperout.nodeTpMat(nodes_nothub,:))./length(nodes_nothub);
    rsn_trs_wtd = [];
    net_names{end+1} = 'NOTA';
    for net = 1:1:length(net_names)
        nodes_dom = find(dom==net);
        nodes_dom = setdiff(nodes_dom, nodes_hub); %excluding hub nodes
        rsn_trs_wtd(net,:) = sum(mapperout.nodeTpMat(nodes_dom,:))./length(nodes_dom);
        
    end
      
    % colorful - across all rsn and hub nodes
    [i,ciu]=max([hub_trs_wtd; rsn_trs_wtd]);
    tmp = tr_label;
    tmp2 = tr_label;
    tmp(tr_label~=0) = ciu;
    tmp2(tr_label~=0) = ciu;
    tmp(tr_label==-1) = [];
    [tm2, mc2] = runMC2(tmp, 0, 'ciu'); % removes bad/stiching frames while counting for transitions
    ciu_s{s} = tmp;
    
    mc_s(s) = mc2;    
    tmp3 = tmp2';
    tmp3(tmp3==0) = []; % remove badframes.
    tr_color_s{s}= tmp3;    
end

clc; close all;
figure; set(gcf,'color','w','Position',[  1         744        2560         448]);

for s=1:1:10
    ax = subtightplot(10,1,s), imagesc(tr_color_s{s}(1:600));
    if s==4
        map = getColorMap([{'HUB'},net_names_s4,{'NOTA'}]);
    elseif s==10
        map = getColorMap([{'HUB'},net_names_s10,{'NOTA'}]);
    else
        map = getColorMap([{'HUB'},net_names_allelse,{'NOTA'}]);
    end
    colormap(ax,map);
    set(gca, 'FontSize',20)
end

%% generating figure 6E

clear; clc; close all;
tmp = load('output/mc_even_hub99.mat');
ev = reshape(tmp.mc_all, [10,144]);
tmp = load('output/mc_odd_hub99.mat');
od = reshape(tmp.mc_all, [10,144]);

mc_nodiag = [];

for s = 1:1:10
    tmp1 = squeeze(od(s,:,:));
    mc_nodiag = [mc_nodiag;  tmp1(:)'];

    tmp1 = squeeze(ev(s,:,:));
    mc_nodiag = [mc_nodiag;  tmp1(:)'];

end

hypo_mat = [1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
            1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
            0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
            0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
            0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0
            0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0
            0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0
            0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0
            0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0
            0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0
            0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0
            0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0
            0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0
            0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1
            0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1
            ];

[r,p]=corrplot(mc_nodiag','type','spearman','testR','on');
close all;
r = suppress_diag(r);
hypo_mat2 = suppress_diag(hypo_mat);
anova1(fisherz([r(hypo_mat2>0);r(hypo_mat2==0)]), [ones(20,1); 2*ones(380,1)])
figure; nhist({r(hypo_mat2>0) r(hypo_mat2==0)})


%% utility functions
function new_map = getColorMap(net_names_needed)
    
    net_names = {'HUB','DMN','VIS','FRP','DAN','MOT','VAN','SAL','COP','SM','AUD','NOTA'};
    map =   [
            201,201,201;
            254,3,2;
            255,250,0;
            255,206,153;
            30,253,3;
            254,176,245;
            71,185,216;
            0,0,0;
            92,22,137;
            16,255,235;
            195,122,193;
            255,255,255;
            ];
    map = map./255;
    new_map = [];
    %new_map = [new_map; map(1,:)];
    for n = 1:1:length(net_names_needed)
        new_map = [new_map; map(find(ismember(net_names,net_names_needed{n})),:)];
    end
    

end
function [tm, mc] = runMC2(ciu, plotFlag, plotTitle)
    % modifying the code to account for frame dropouts and session
    % boundaries
    
    nciu = max(ciu);%sum(unique(ciu)>0);
    tm = zeros(nciu, nciu) + 1e-10; % a small value is added for subjects which don't have all networks.
    for i = 1:1:length(ciu)-1
        j = i + 1;
        if ciu(i) <= 0 || ciu(j) <=0
            %fprintf(1,'skipping due to 0\n');
            continue;
        else
            tm(ciu(i),ciu(j)) = tm(ciu(i),ciu(j)) + 1;
        end

    end
    
    mc = dtmc(tm);
    
    if plotFlag == 1
        figure; set(gcf, 'Position',[596         496        1722         408]);
        subplot(1,3,1),graphplot(mc,'ColorNodes',true,'ColorEdges',true,'LabelEdges',true);
        subplot(1,3,2),eigplot(mc);
        X = redistribute(mc,25);
        subplot(1,3,3),distplot(mc,X,'Type','Histogram','FrameRate',0.01);
        suptitle(replace(plotTitle,'_','-')); 
    end

end
function [hubness95, hubness99] = computeHubness(mat, doi, deg)
    g = graph(mat);
    closeness_matlab = centrality(g, 'closeness');   

    hubness99 = zeros(size(mat,1),1);
    hubness99(intersect(find(deg>doi(1)), find(closeness_matlab>prctile(closeness_matlab,99)))) = 1;

    hubness95 = zeros(size(mat,1),1);
    hubness95(intersect(find(deg>doi(1)), find(closeness_matlab>prctile(closeness_matlab,95)))) = 1;
    
end
function scaled = zscoreScaling(value)
    scaled = zscore(value);
    scaled = scaled./max(scaled);
end
function [dom1, corrProp] = calcDominance(nodeTpMat, nodeBynode, netData, thr)
    numMetaGroups = size(netData,2)+1; %last one is for neither of the 11 networks
    numNodes = size(nodeBynode, 1);
    prop = zeros(numNodes, numMetaGroups); 
    for n = 1:1:numNodes
        trs = find(nodeTpMat(n,:));

        for tr = trs
            if any(netData(tr,:)>thr)
                prop(n, netData(tr,:)>thr) = prop(n, netData(tr,:)>thr) + 1;
            else
                prop(n, numMetaGroups) = prop(n, numMetaGroups) + 1;
            end
        end
        
       
    end
    
    [~, dom1] = max(prop, [], 2);
    corrProp = corr(prop);
end
function [entropy1] = calcEntropy(nodeTpMat, nodeBynode, netData)

    numNodes = size(nodeBynode, 1);
    entropy1 = [];
    for n = 1:1:numNodes
        trs = find(nodeTpMat(n,:));
        if length(trs) > 1
            mean_nets = mean(netData(trs,:));        
            entropy1(n) = std(mean_nets);
        else
            entropy1(n) = -1;            
        end
    end
end
function parsave(myfile, sbj_name, runType, metricType, nodeTpMat, nodeBynode, tpMat, parcelData,  networkDataVertex, ses_idx, sess_id, filter)   
    matObj = matfile(myfile, 'Writable', true);
    matObj.nodeTpMat = nodeTpMat;
    matObj.nodeBynode = nodeBynode;
    matObj.tpMat = tpMat;
    matObj.parcelData = parcelData;
    matObj.networkDataVertex = networkDataVertex;
    matObj.sbj_name = sbj_name;
    matObj.ses_idx = ses_idx;
    matObj.sess_id = sess_id;
    matObj.metricType = metricType;
    matObj.runType = runType;
    matObj.filter = filter;
end
function [nodeTpMat, nodeBynode, tpMat, filter] = runBDLMapper_wrapper(parcelData, metricType)

    % run mapper
    data_z = parcelData;
    nclust = 10;
    num_bin_clusters = nclust;
    num_k = 30; 
    res_val = 30; 
    gain_val = 70; 
    
    [nodeTpMat, nodeBynode, tpMat, filter] = runBDLMapper(data_z, metricType, res_val, gain_val, num_k, num_bin_clusters);
            
end

function [nodeTpMat, nodeBynode, tpMat, filter] = runBDLMapper(data, metricType, res_val, gain_val, num_k, num_bin_clusters)

    X = data;


    resolution = [res_val res_val];
    gain = gain_val;
    
    fprintf(1,'Estimating distance matrix\n');
    tic
    distMat = estimateDistance(X, metricType);
    toc
        
    fprintf(1,'Estimating knn graph\n');
    tic
    % create knn graph, estimate geodesic distances, embed using cmdscale and apply mapper
    [knnGraphTbl, knnGraph_dense_bin, knnGraph_dense_wtd, knnGraph_dense_bin_conn, knnGraph_dense_wtd_conn]= createPKNNG_bdl(distMat, num_k);

    knn_g_wtd = graph(knnGraph_dense_bin_conn);

    % estimate geodesic distances
    dist_geo_wtd = round(distances(knn_g_wtd,'Method','positive'));
    toc
    
    fprintf(1,'Estimating embedding\n');
    tic
    
    % embed using cmdscale
    [y,e] = cmdscale(dist_geo_wtd);

    filter = [y(:,1), y(:,2)];
    toc
    
    fprintf(1,'Running mapper\n');
    tic
    
    
    
    [adja, adja_pruned, pts_in_vertex, pts_in_vertex_pruned] = mapper2d_bdl_hex_binning(distMat, filter, resolution, gain, num_bin_clusters, 3); 
    
    toc
    
    fprintf(1,'Creating final output\n');
    tic
    
    % creating matrices for d3 visualization
    numNodes = length(pts_in_vertex_pruned);
    numTp = size(X,1);
    nodeTpMat = zeros(numNodes, numTp);
    for node = 1:1:numNodes
        tmp = pts_in_vertex_pruned{node};
        nodeTpMat(node, tmp) = 1;
    end

    nodeBynode = adja_pruned;
    tpMat = getMatTp_wtd(nodeBynode, nodeTpMat);
    fprintf(1,'Done\n');
    toc

end
function distMat = estimateDistance(X, metricType)
    distMat = squareform(pdist(X, metricType));
end
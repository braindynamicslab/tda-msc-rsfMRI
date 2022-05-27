%WGGGG
%
%
% aims - 
% [1] Runs Mapper on odd and even runs of MSC data and
% saves the mapperout file that can be used later for analysis and viz.
% [2] Runs corresponding null models for comparison.
% 
% Please cite:
% Saggar, M., Shine, J.M., Liegeois, R., Dosenbach, N.U.F., Fair, D. 2021. 
% Precision dynamical mapping using topological data analysis reveals a 
% unique hub-like transition state at rest. BioRxiv
%
% date - 5.31.2021
% author - saggar@stanford.edu
%
%% load data
cd ~/data
clear; clc; close all;
metricType = 'euclidean';
runType = 'odd';
output_name = 'mapperout';

for s = 1:1:10 
    sbj_name = sprintf('sub-MSC%02d',s);
    cd(sbj_name);
    
    ses_files = dir('vc*.dtseries.nii');
    ses_files = {ses_files.name};
    fd_files = dir('*TMASK.txt');
    fd_files = {fd_files.name};

    parcel_ids = cifti_read(sprintf('%s_parcels.dtseries.nii', sbj_name));
    parcel_ids = parcel_ids.cdata;
    vertex_net = cifti_read(sprintf('%s_parcel_networks.dscalar.nii', sbj_name));
    vertex_net = vertex_net.cdata; 

    vertexNetworksIds = unique(vertex_net);
    parcelNetworkIds = unique(parcel_ids);
    
    net_ids_in_parcels = {};
    for net = 1:1:length(vertexNetworksIds)
        tmp = unique(parcel_ids(find(vertex_net == vertexNetworksIds(net))));
        [c,ia,ib] = intersect(parcelNetworkIds, tmp);
        net_ids_in_parcels{net} = ia;
    end
    parcel_order=[];
    for n=1:1:length(net_ids_in_parcels) 
        parcel_order= [parcel_order; net_ids_in_parcels{n}]; 
    end
    networkDataVertex = [];
    parcelData = [];
    parcelData_noTmask = [];

    X_odd = [];
    
    sess_id = [];
    if strcmp(runType,'odd') == 1
        ses_idx = 1:2:10
    elseif strcmp(runType,'even') == 1
        ses_idx = 2:2:10
    end
    
    fd_all_ses = [];
    for ses = ses_idx
        ses_files{ses}
        st = cifti_read(ses_files{ses});
        data = st.cdata;
        fd  = load(fd_files{ses});
        fd_all_ses = [fd_all_ses; fd];

        data_notMask = data;
        data = data(:,fd>0);
        
        X_odd = [X_odd; zscore(data')];

        sess_id = [sess_id; ones(size(data,2),1)*ses];

        % create Network based mean time series for vertexwise data
        tmp = [];
        tmp2 = [];
        for p = 1:1:length(vertexNetworksIds)
           tmp = find(vertex_net == vertexNetworksIds(p));
           tmp2(p,:) = mean(data(tmp,:),1);
        end
        tmp2 = zscore(tmp2');
        tmp2(:,1) = [];
        if s==4 % accounting for two misssing networks in s4
            tmp2 = [tmp2(:,1), tmp2(:,2), tmp2(:,3), tmp2(:,5), tmp2(:,6), tmp2(:,9), mean(tmp2(:,[10,11]),2), tmp2(:,12)];
        elseif s==10 % accounting for one misssing network in s10
            tmp2 = [tmp2(:,1), tmp2(:,2), tmp2(:,3), tmp2(:,5), tmp2(:,6), tmp2(:,7), tmp2(:,9), mean(tmp2(:,[10,11]),2), tmp2(:,12)];
        else
            tmp2 = [tmp2(:,1), tmp2(:,2), tmp2(:,3), tmp2(:,5), tmp2(:,6), tmp2(:,7), tmp2(:,8), tmp2(:,9), mean(tmp2(:,[10,11]),2), tmp2(:,12)];
        end
        
        networkDataVertex = [networkDataVertex; tmp2];

        % create parcel based mean timeseries
        tmp = [];
        tmp2 = [];
        tmp3 = [];
        for p = 1:1:length(parcelNetworkIds)
           tmp = find(parcel_ids == parcelNetworkIds(p));
           tmp2(p,:) = mean(data(tmp,:),1);
           tmp3(p,:) = mean(data_notMask(tmp,:),1);
        end
        tmp2 = zscore(tmp2');
        tmp3 = zscore(tmp3');
        parcelData = [parcelData; tmp2];
        parcelData_noTmask = [parcelData_noTmask; tmp3];
        
        
    end
    
    % real data 
    [nodeTpMat, nodeBynode, tpMat, filter] = runBDLMapper_wrapper(parcelData, metricType);
    myfile = sprintf('%s_%s_runs_mat_metric_%s_%s.mat', sbj_name, runType, metricType, output_name);
    parsave(myfile, sbj_name, runType, metricType, nodeTpMat, nodeBynode, tpMat, parcelData,  parcel_order, networkDataVertex, ses_idx, sess_id, filter);
    
    % null data
    num_surr = 25;
    surr_pr = [];
    ar_model_order = 1;
    % using CBIG code to generate null datasets.
    surr_ar = CBIG_RL2017_get_AR_surrogate(parcelData, num_surr, ar_model_order, 'gaussian', length(fd_all_ses), parcelData_noTmask);
    surr_pr = CBIG_RL2017_get_PR_surrogate(parcelData, num_surr);
    
    for it = 1:1:num_surr
        it
        parcelData_null = squeeze(surr_pr(:,:,it));
        sbj_name_null = sprintf('%s_PR_null%03d',sbj_name,it);
        
        
        %creating metaInfo for null data
        parcelDataNetworks_null = [];
        tmp = [];
        tmp2 = [];
        for net = 1:1:length(vertexNetworksIds)
            tmp2(:,net) = mean(parcelData_null(:,net_ids_in_parcels{net}),2);
        end
        tmp2 = zscore(tmp2);
        tmp2(:,1) = [];
        if s==4
            tmp2 = [tmp2(:,1), tmp2(:,2), tmp2(:,3), tmp2(:,5), tmp2(:,6), tmp2(:,9), mean(tmp2(:,[10,11]),2), tmp2(:,12)];
        elseif s==10
            tmp2 = [tmp2(:,1), tmp2(:,2), tmp2(:,3), tmp2(:,5), tmp2(:,6), tmp2(:,7), tmp2(:,9), mean(tmp2(:,[10,11]),2), tmp2(:,12)];
        else
            tmp2 = [tmp2(:,1), tmp2(:,2), tmp2(:,3), tmp2(:,5), tmp2(:,6), tmp2(:,7), tmp2(:,8), tmp2(:,9), mean(tmp2(:,[10,11]),2), tmp2(:,12)];
        end        
        parcelDataNetworks_null = [parcelDataNetworks_null; tmp2];
        

        try
            [nodeTpMat_null, nodeBynode_null, tpMat_null, filter_null] = runBDLMapper_wrapper(parcelData_null, metricType);        
        catch
            fprintf(2,'Error due to CMD SCale in PR it%03d', it);
            continue;
        end
        
        myfile = sprintf('%s_%s_runs_mat_metric_%s_%s', sbj_name_null, runType, metricType, output_name);
        parsave(myfile, sbj_name_null, runType, metricType, nodeTpMat_null, nodeBynode_null, tpMat_null, parcelData_null, parcel_order, parcelDataNetworks_null, ses_idx, sess_id, filter_null)
        
        close all;
    end
    

    for it = 1:1:num_surr
        it
        parcelData_null_ar = squeeze(surr_ar(:,:,it));

        % Jul 22, 2021 - adding temporal mask to null to see if that
        % changes degree distribution
        parcelData_null_ar = parcelData_null_ar(fd_all_ses>0,:);
        
        sbj_name_null = sprintf('%s_AR_null%03d',sbj_name,it);
        
        %creating metaInfo for null data
        parcelDataNetworks_null_ar = [];
        tmp = [];
        tmp2 = [];
        for net = 1:1:length(vertexNetworksIds)
            tmp2(:,net) = mean(parcelData_null_ar(:,net_ids_in_parcels{net}),2);
        end
        tmp2 = zscore(tmp2);
        tmp2(:,1) = [];
        if s==4
            tmp2 = [tmp2(:,1), tmp2(:,2), tmp2(:,3), tmp2(:,5), tmp2(:,6), tmp2(:,9), mean(tmp2(:,[10,11]),2), tmp2(:,12)];
        elseif s==10
            tmp2 = [tmp2(:,1), tmp2(:,2), tmp2(:,3), tmp2(:,5), tmp2(:,6), tmp2(:,7), tmp2(:,9), mean(tmp2(:,[10,11]),2), tmp2(:,12)];
        else
            tmp2 = [tmp2(:,1), tmp2(:,2), tmp2(:,3), tmp2(:,5), tmp2(:,6), tmp2(:,7), tmp2(:,8), tmp2(:,9), mean(tmp2(:,[10,11]),2), tmp2(:,12)];
        end        
        parcelDataNetworks_null_ar = [parcelDataNetworks_null_ar; tmp2];
        
        try
            [nodeTpMat_null_ar, nodeBynode_null_ar, tpMat_null_ar, filter_null_ar] = runBDLMapper_wrapper(parcelData_null_ar, metricType);        
        catch
            fprintf(2,'Error due to CMD SCale in AR it%03d', it);
            continue;
        end
        
        myfile = sprintf('%s_%s_runs_mat_metric_%s_%s', sbj_name_null, runType, metricType, output_name);
        parsave(myfile, sbj_name_null, runType, metricType, nodeTpMat_null_ar, nodeBynode_null_ar, tpMat_null_ar, parcelData_null_ar,  parcel_order, parcelDataNetworks_null_ar, ses_idx, sess_id, filter_null_ar)
        
        close all;
    end
    toc
end


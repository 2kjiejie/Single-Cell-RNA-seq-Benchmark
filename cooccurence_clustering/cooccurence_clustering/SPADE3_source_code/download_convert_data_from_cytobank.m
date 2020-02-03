function [is_successful] = download_convert_data_from_cytobank(directoryname)

is_successful = 0;
try
    run(fullfile(directoryname,'cytobank_data_source.m'));
catch
    fid = fopen(fullfile(directoryname,'cytobank_data_source.m'));
    d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
    d(d==char(13))=[];
    d = regexp(d,char(10),'split');
    for i=1:length(d)
        if isempty(setdiff(d{i},' ')) || d{i}(1)=='%'
            continue;
        end
        tmp = d{i};
        tmp(tmp==' ') = [];
        tmp(tmp=='''') = [];
        tmp(tmp==';') = [];
        tmp = regexp(tmp, '=', 'split');
        if isequal(tmp{1},'CYTOBANK')
            CYTOBANK = tmp{2};
        end
        if isequal(tmp{1},'username')
            username = tmp{2};
        end
        if isequal(tmp{1},'password')
            password = tmp{2};
        end
        if isequal(tmp{1},'EXPT_ID')
            EXPT_ID = tmp{2};
        end
        if isequal(tmp{1},'SPADE_ID')
            SPADE_ID = tmp{2};
        end
    end
end


% authenticate
tic; fprintf('Connecting to cytobank (~5 seconds) ...')
url = ['https://',CYTOBANK,'.cytobank.org/cytobank/api/v1/authenticate'];
uri = matlab.net.URI(url);
request = matlab.net.http.RequestMessage;
request.Header = matlab.net.http.HeaderField('accept', 'application/json', 'charset', 'UTF-8', 'content-type', 'application/json');
request.Method = 'post';
request.Body = matlab.net.http.MessageBody;
request.Body.Payload = ['{ "username": "', username, '", "password": "', password, '" }'];
response = sendRequest(uri,request);
if ~isequal(response.StatusCode,200)
    is_successful = 0;
    return
end
response.Body.Data.user;
authToken = response.Body.Data.user.authToken;
toc



if exist('SPADE_ID')  
    % download SPADE analysis results, full package
    tic; fprintf('Downloading SPADE results from Cytobank to local (number of files * 40~240 seconds per file) ...')
    url = ['https://',CYTOBANK,'.cytobank.org/cytobank/api/v1/experiments/',EXPT_ID,'/advanced_analyses/spade/',SPADE_ID,'/download?item=full_data'];
    uri = matlab.net.URI(url);
    request = matlab.net.http.RequestMessage;
    request.Header = matlab.net.http.HeaderField('authorization', authToken);
    request.Method = 'get';
    response = sendRequest(uri,request);
    header = jsondecode(jsonencode(response.Header));
    header_names = {header.Name};
    header_values = {header.Value};
    tmp = header_values{find(ismember(header_names, 'Content-Disposition'))};
    filename = tmp(strfind(tmp,'filename')+10:end-1);
    fileID = fopen(fullfile(directoryname,filename),'w');
    fwrite(fileID,response.Body.Data,'uint8');
    fclose(fileID);
    toc

    % unzip the downloaded data 
    tic; fprintf('Unzip downloaded data ...')
    mkdir(fullfile(directoryname,'downloaded_from_cytobank'))
    unzip(fullfile(directoryname,filename),fullfile(directoryname,'downloaded_from_cytobank'))
    delete(fullfile(directoryname,filename))
    toc

    % convert downloaded SPADE results into matlab format
    tic; fprintf('Convert downloaded data to matlab format ...');

    source_directory = fullfile(directoryname,'downloaded_from_cytobank');
    target_directory = directoryname;
    parameter_filename = 'SPADE_parameters.mat';
    pooled_downsampled_filename = 'SPADE_pooled_downsampled_data.mat';
    cluster_mst_upsample_filename='SPADE_cluster_mst_upsample_result.mat';


    table_content = read_xls_v2(fullfile(source_directory,'clusters.table'), char(0), 0);
    number_of_desired_clusters = length(table_content)-1; 
    used_channels = regexp(table_content{1},'"','split'); used_channels = used_channels(2:2:end);

    source_filenames = getfilenames(source_directory);
    source_filenames = source_filenames(isInListEnd(source_filenames, 'cluster.fcs'));
    for i=1:length(source_filenames)
        [data, marker_names, channel_names, scaled_data, compensated_data, fcshdr] = readfcs_v2(fullfile(source_directory, source_filenames{i}));
        [C,IA,IB] = intersect(used_channels,channel_names);
        used_markers = marker_names(sort(IB));
        local_density = data(ismember(marker_names,'density'),:);
        all_assign{i} = data(ismember(marker_names,'cluster'),:);
        writefcs_v2(fullfile(target_directory, source_filenames{i}(1:strfind(source_filenames{i},'.density.fcs.cluster.fcs')-1)), data(1:end-2,:), marker_names(1:end-2,:),channel_names(1:end-2,:));
    end


    all_fcs_filenames = getfilenames(directoryname);
    all_fcs_filenames = all_fcs_filenames(isInListEnd(all_fcs_filenames,'.fcs'));
    all_markers = cell(0);
    all_overlapping_markers = cell(0);
    for i=1:length(all_fcs_filenames)
        filename = fullfile(directoryname,all_fcs_filenames{i});
        [fcshdr] = fca_readfcs_only_header(filename);
        display(['Get marker names from fcs file:'])
        display(['   ',all_fcs_filenames{i}]);
        marker_names = cell(0);
        for i=1:length(fcshdr.par), 
            marker_names{i,1} = fcshdr.par(i).name2;  
            if isequal(unique(fcshdr.par(i).name2),' ')
                marker_names{i,1} = fcshdr.par(i).name;  
            end
        end
        [C,I] = setdiff(marker_names, all_markers);
        if ~isempty(I)
            all_markers = [all_markers;marker_names(sort(I))];
        end
        if length(all_overlapping_markers)==0
            all_overlapping_markers = marker_names;
        else
            [C,IA,IB] = intersect(all_overlapping_markers, marker_names);
            all_overlapping_markers = all_overlapping_markers(sort(IA));
        end
    end


    apply_compensation = 0;
    arcsinh_cofactor = 5;
    clustering_algorithm = 'kmeans';
    density_estimation_optimization_factor = 1.5;
    file_annot = all_fcs_filenames;
    file_used_to_build_SPADE_tree = all_fcs_filenames;
    kernel_width_factor = 5;
    max_allowable_events = 5000000;
    outlier_density = 1;
    target_cell_number = 20000;
    target_density = 10;
    target_density_mode = 1; % 1 means using target density percentile, 2 means choose a TD such that a fixed number of cells survive downsampling
    transformation_option = 1; % 0 means no transformation, 1 means arcsinh, 2 means arcsinh followed by 0-mean 1-var
    used_markers = all_markers(ismember(all_markers, used_markers));
    save(fullfile(directoryname, parameter_filename),... 
        'all_fcs_filenames', 'all_markers','all_overlapping_markers', 'apply_compensation', ...
        'arcsinh_cofactor', 'clustering_algorithm', 'density_estimation_optimization_factor', ...
        'file_annot', 'file_used_to_build_SPADE_tree', 'kernel_width_factor', 'max_allowable_events', ...
        'number_of_desired_clusters', 'outlier_density', 'target_cell_number', 'target_density', ...
        'target_density_mode', 'transformation_option', 'used_markers');


    % generate local density data files
    for i=1:length(all_fcs_filenames)
        [data, marker_names, channel_names] = readfcs_v2(fullfile(directoryname, all_fcs_filenames{i}));
        local_density = data(ismember(marker_names,'density'),:);
        data = flow_arcsinh(data(1:end,:), arcsinh_cofactor);
        marker_names = marker_names(1:end,:);
        kernel_width = 1;
        save(fullfile(directoryname, [all_fcs_filenames{i}(1:end-3),'mat']), 'data', 'kernel_width', 'local_density', 'marker_names', 'used_markers');
    end

    % generate the pooled downsampled data
    downsampled_source_filenames = getfilenames(source_directory);
    downsampled_source_filenames = downsampled_source_filenames(isInListEnd(downsampled_source_filenames, 'downsample.fcs'));
    all_data=[];
    tube_channel=[];
    all_local_density=[];
    idx = [];
    all_marker_names = all_markers;
    used_marker_names = used_markers;
    for i=1:length(downsampled_source_filenames)
        [data, marker_names, channel_names] = readfcs_v2(fullfile(source_directory, source_filenames{i}));
        [data2, marker_names2, channel_names2] = readfcs_v2(fullfile(source_directory, downsampled_source_filenames{i}));
        keep_ind = ismember(data(ismember(marker_names,'EventNum'),:), data2(ismember(marker_names2,'EventNum'),:));
        idx = [idx,  data(ismember(marker_names,'cluster'),keep_ind)];
        local_density = data(ismember(marker_names,'density'),keep_ind);
        data = flow_arcsinh(data(:, keep_ind), arcsinh_cofactor);

        if isequal(marker_names,all_marker_names)
            all_data = [all_data, data];
        else
            data_tmp = zeros(length(all_marker_names),size(data,2))+NaN;
            [C,IA,IB] = intersect(marker_names,all_marker_names);
            data_tmp(IB,:) = data(IA,:);
            all_data = [all_data, data_tmp];
        end
        all_local_density = [all_local_density,local_density];
        tube_channel = [tube_channel,repmat(i,1,size(data,2))];
    end
    all_data = [all_data;tube_channel];
    all_marker_names{end+1} = 'FileInd';

    data = all_data; 
    marker_names = all_marker_names;
    local_density = all_local_density;
    save(fullfile(directoryname, pooled_downsampled_filename), 'all_data', 'all_local_density', 'data', 'local_density', 'marker_names', 'used_markers');

    % generate the mst results
    fid = fopen(fullfile(source_directory,'mst.gml')); d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
    source_ind = strfind(d,'source');
    target_ind = strfind(d,'target');
    mst_tree = sparse(number_of_desired_clusters, number_of_desired_clusters);
    for k=1:length(source_ind)
        a = str2num(d(source_ind(k)+6 : source_ind(k)+min(strfind(d(source_ind(k):end),char(10)))))+1;
        b = str2num(d(target_ind(k)+6 : target_ind(k)+min(strfind(d(target_ind(k):end),char(10)))))+1;
        mst_tree(a,b) = 1;
        mst_tree(b,a) = 1;
    end

    table_content = read_xls_v2(fullfile(source_directory,'layout.table'),' ',0);
    node_positions = cell2num(table_content)';
    node_positions = node_positions - repmat((max(node_positions,[],2)+min(node_positions,[],2))/2,1,size(node_positions,2));
    node_positions = node_positions/max(abs(node_positions(:)))*40;

    node_size = zeros(1,size(node_positions,2));
    for i=1:length(node_size), node_size(i) = sum(local_density(idx==i)); end
    node_size = flow_arcsinh(node_size, median(node_size)/2);
    node_size = ceil(node_size/max(node_size)*10);
    node_size(node_size<5)=5;
    node_size = node_size * 1.2;

    tree_annotations = [];
    tree_bubble_contour = [];
    save(fullfile(directoryname, cluster_mst_upsample_filename), 'all_assign','all_fcs_filenames','file_annot', 'data', 'local_density', 'marker_names', 'used_markers','idx','mst_tree','node_positions','node_size','tree_annotations','tree_bubble_contour');


    marker_node_average=cell(0); counter = 1;
    % get the average from the pooled data
    for j=1:length(marker_names)
        marker_node_average{counter,1} = 'POOLED';
        marker_node_average{counter,2} = marker_names{j};
        [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(data(j,:), idx);
        group_avg(group_idx_values==0)=[];
        counts(group_idx_values==0)=[];
        group_idx_values(group_idx_values==0)=[];
        tmp = zeros(1,max(idx))+NaN;
        tmp(group_idx_values) = group_avg;
        marker_node_average{counter,3} = tmp;
        counter = counter + 1;
    end
    marker_node_average{counter,1} = 'POOLED';
    marker_node_average{counter,2} = 'CellFreq';
    [dummy, tmp] = SPADE_compute_one_marker_group_mean(ones(1,length(idx)),idx);    
    marker_node_average{counter,3} = tmp(:)';
    counter = counter + 1;
    % get the average from individual files
    for i=1:length(all_fcs_filenames)
        fprintf('\b\b\b\b\b%5d',i+1);
        load(fullfile(directoryname, [all_fcs_filenames{i}(1:end-3),'mat']),'data','marker_names');
        for j=1:length(marker_names)
            marker_node_average{counter,1} = file_annot{i};
            marker_node_average{counter,2} = marker_names{j};
            [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(data(j,:), all_assign{i}); % the following few lines are for the purpose that: some file may not have any cell belong to one particular node, and therefore, the "group_avg" does not have information for every node
            group_avg(group_idx_values==0)=[];
            counts(group_idx_values==0)=[];
            group_idx_values(group_idx_values==0)=[];
            tmp = zeros(1,max(idx))+NaN;
            tmp(group_idx_values) = group_avg;
            marker_node_average{counter,3} = tmp;
            counter = counter + 1;
        end    
        marker_node_average{counter,1} = file_annot{i};
        marker_node_average{counter,2} = 'CellFreq';
        [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(ones(1,length(all_assign{i})),all_assign{i}); 
        marker_node_average{counter,3} = zeros(1,max(idx));
        marker_node_average{counter,3}(group_idx_values)=counts;
        counter = counter + 1;
    end

    save(fullfile(directoryname, cluster_mst_upsample_filename), 'marker_node_average', '-append');

    toc

    is_successful = 1;

else  % download data files only
    
    
    % get the list of fcs files
    tic; fprintf('Downloading FCS files from Cytobank to local ... \n')
    url = ['https://',CYTOBANK,'.cytobank.org/cytobank/api/v1/experiments/',EXPT_ID,'/fcs_files'];
    uri = matlab.net.URI(url);
    request = matlab.net.http.RequestMessage;
    request.Header = matlab.net.http.HeaderField('authorization', authToken);
    request.Method = 'get';
    response = sendRequest(uri,request);
    fcsFiles_info = response.Body.Data.fcsFiles;
    
    for i=1:length(fcsFiles_info)
        tic
        filename = fcsFiles_info(i).filename;
        file_ID = fcsFiles_info(i).id;
        fprintf(['Downloading FCS file ', num2str(i), '/', num2str(length(fcsFiles_info)), ' ', filename,' (~30 seconds) ... '])
        
        url = ['https://',CYTOBANK,'.cytobank.org/cytobank/api/v1/experiments/',EXPT_ID,'/fcs_files/',num2str(file_ID),'/download'];
        uri = matlab.net.URI(url);
        request = matlab.net.http.RequestMessage;
        request.Header = matlab.net.http.HeaderField('authorization', authToken);
        request.Method = 'get';
        response = sendRequest(uri,request);
        file_write_ID = fopen(fullfile(directoryname,filename),'w');
        fwrite(file_write_ID,response.Body.Data,'uint8');
        fclose(file_write_ID);
        toc
    end
    
    is_successful = 1;
    
end



% logout
tic; fprintf('Disconnect from cytobank (~3 seconds)...')
url = ['https://', CYTOBANK, '.cytobank.org/cytobank/api/v1/logout'];
uri = matlab.net.URI(url);
request = matlab.net.http.RequestMessage;
request.Header = matlab.net.http.HeaderField('authorization', authToken);
request.Method = 'post';
request.Body = matlab.net.http.MessageBody;
request.Body.Payload = ' '; 
response = sendRequest(uri,request);
response.Body.Data.message;
toc






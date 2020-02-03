function [is_successful] = upload_annotation_to_cytobank(directoryname, tree_annotations)

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


if ~exist('SPADE_ID')  
    fprintf('SPADE_ID not specified, cannot upload bubbles to cytobank !! \n\n');
    return
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



% copy SPADE results to create a new SPADE
tic; fprintf('Generate another copy of the SPADE analysis in cytobank ...')
url = ['https://',CYTOBANK,'.cytobank.org/cytobank/api/v1/experiments/',EXPT_ID,'/advanced_analyses/spade/',SPADE_ID,'/copy_results'];
uri = matlab.net.URI(url);
request = matlab.net.http.RequestMessage;
request.Header = matlab.net.http.HeaderField('authorization', authToken);
request.Method = 'post';
request.Body = matlab.net.http.MessageBody;
request.Body.Payload = ' '; 
response = sendRequest(uri,request);
new_SPADE_ID = num2str(response.Body.Data.spade.id)
toc


% upload bubbles
tic; fprintf('Upload bubbles to cytobank ...')
clear spade
for i=1:length(tree_annotations)
    clear tmp
    tmp.name = ['annotation - ',num2str(i)];
    tmp.nodes = tree_annotations{i};
    spade.bubbles(i) = tmp;
end    
clear payload
payload.spade = spade;

url = ['https://',CYTOBANK,'.cytobank.org/cytobank/api/v1/experiments/',EXPT_ID,'/advanced_analyses/spade/',new_SPADE_ID,'/set_bubbles'];
uri = matlab.net.URI(url);
request = matlab.net.http.RequestMessage;
request.Header = matlab.net.http.HeaderField('authorization', authToken, 'content-type', 'application/json');
request.Method = 'post';
request.Body = matlab.net.http.MessageBody;
request.Body.Payload = jsonencode(payload);
response = sendRequest(uri,request);
toc



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





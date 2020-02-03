function [c,Q,bestpartition,bestpartitionhierarchy] = LouvainfromBin( filename, numiters )
%the path to louvain in each cyt
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
curr_path = [curr_path filesep];
ps = [curr_path 'Louvain' filesep];

filename = strrep( filename, '.bin', '' );
if exist([filename '.bin'])~=2
    error('G.bin file does not exist yet!')
end
% begin
%fprintf(1, 'MATLAB: calling convert:\n');
command = [ps 'convert -i ' filename '.bin -o ' filename '_graph.bin -w ' filename '_graph.weights' ];
%fprintf(1,'%s\n', command );
system( command );

% run community detection
fprintf(1,'Running community detection, ITERATION %4d', 0);
for iter = 1:numiters    
    fprintf(1,'\b\b\b\b%4d', iter);
    
    command = [ps 'community ' filename '_graph.bin -l -1 -v -w ' filename '_graph.weights > ' filename '.tree'];
    %fprintf( 1, '%s\n', command );
    [~,r] = system( command );
    %fprintf(1, '\n');
    try
        % find each iteration's modularity
        q = find_modularity( r );
        % fprintf( 1, 'MATLAB: modularity scores:\n' );
    catch
        cleanup(filename)
        error('Unable to find modularity score in the stderr: %s.\nCheck that the correct path to the Louvain code is specified.', r);
    end

    
    % find number of lvevls
    command = [ps 'hierarchy ' filename '.tree' ];
    %fprintf(1, '%s\n', command );
    [~,r] = system( command );
    %fprintf(1, '\n' );
    
    r = strtok(r, 10);
    r = regexprep( r, 'Number of levels: ', '' );
    num_levels = str2double( r )-1;
    
    %fprintf( 1, 'MATLAB: max level is %d\n', num_levels );
    
    % import levels
    for level = 1:num_levels
        %fprintf( 1, 'MATLAB: importing level %d\n', level );
        command = [ps 'hierarchy ' filename '.tree -l ' num2str( level ) ' > ' filename '.tmp' ];
        %fprintf(1, '%s\n', command );
        system( command );
        hierarchy_output = load( [filename '.tmp'] );
        c{iter,level} = hierarchy_output(:,2) + 1;
        Q{iter,level} = q(level);
    end
    
end
fprintf(1,'\n');

% find best partition
maxmod = -1;
for i = 1:numel(Q)
    if Q{i} > maxmod
        maxmod = Q{i};
        [I,J] = ind2sub( size(Q), i );
    end
end
bestpartition = c{I,J};
bestpartitionhierarchy = c(I,:);
% delete temporary files
cleanup(filename)

%-------------------------
function Q = find_modularity( r )
% Q = find_modularity( r )
% convert the text output into modularity score of each iteration
signature = '  modularity increased from %f to %f';
idx = 0;
while( ~isempty( r ) )
    % read a line and match it to the signature
    [token, r] = strtok( r, char( 10 ) );
    a = sscanf( token, signature );
    
    if( ~isempty( a ) )
        % signature matched copy new modularity
        idx = idx + 1;
        Q( idx ) = a( 2 );
    end
end
%-------------------------
function cleanup(filename)
files{1} = dir([filename '*.tmp']);
files{2} = dir([filename '*.tree']);
files{3} = dir([filename '*.weights']);
files{4} = dir([filename '*.bin']);
for i = 1:4
    for j = 1:length(files{i})
        delete( files{i}(j).name )
    end
end


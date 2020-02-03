function mappedX = fast_tsne(X, initial_dims, perplexity, theta)
%FAST_TSNE Runs the (landmark) C++ implementation of t-SNE
%
%   mappedX = fast_tsne(X, initial_dims, perplexity, theta)
%
% Runs the C++ implementation of Barnes-Hut-SNE. The high-dimensional 
% datapoints are specified in the NxD matrix X. The dimensionality of the 
% datapoints is reduced to initial_dims dimensions using PCA (default = 50)
% before t-SNE is performed. Next, t-SNE reduces the points to two 
% dimensions. The perplexity of the input similarities may be specified
% through the perplexity variable (default = 30). The variable theta sets
% the trade-off parameter between speed and accuracy: theta = 0 corresponds
% to standard, slow t-SNE, while theta = 1 makes very crude approximations.
% Appropriate values for theta are between 0.1 and 0.7 (default = 0.5).
% The function returns the two-dimensional data points in mappedX.
%
% NOTE: The function is designed to run on large (N > 5000) data sets. It
% may give poor performance on very small data sets (it is better to use a
% standard t-SNE implementation on such data).


% Copyright (c) 2013, Laurens van der Maaten (Delft University of Technology)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the Delft University of Technology.
% 4. Neither the name of the Delft University of Technology nor the names of 
%    its contributors may be used to endorse or promote products derived from 
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
% EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
% IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
% OF SUCH DAMAGE.


    if ~exist('initial_dims', 'var') || isempty(initial_dims)
        initial_dims = 50;
    end
    if ~exist('perplexity', 'var')
        perplexity = 30;
    end
    if ~exist('theta', 'var')
        theta = 0.5;
    end
    
    work_dir = pwd;
    
    [curr_path, ~, ~] = fileparts(mfilename('fullpath'));
    curr_path = [curr_path filesep];
    
    cd(curr_path);
    try
        
        %checking for old tsne results for the same data
        fileCheck= exist ('tsneResults.mat', 'file');
        
        % create or load cache
        if (fileCheck==0) 
            mapMat = containers.Map(); % hash map for matrix name and tsne output
        else
            file= load('tsneResults.mat'); %loading the old tsne results
            mapMat=file.mapMat;
        end

        % getting the  matrix hash
        hashMat= DataHash(X+perplexity+theta); 

        % check for has in cache -> no need to run tsne again
        if (isKey(mapMat,hashMat)) 
            
            value=values(mapMat,{hashMat});
            mappedX=value{1};
            
            cd(work_dir); % cd to original directory
            
%             fprintf('\nfast_tsne: cached results found.');
            
            % returning cached result
            return;  
        end
        
        % Perform the initial dimensionality reduction using PCA
        X = double(X);
        X = bsxfun(@minus, X, mean(X, 1));
        covX = X' * X;
        [M, lambda] = eig(covX);
        [~, ind] = sort(diag(lambda), 'descend');
        if initial_dims > size(M, 2)
            initial_dims = size(M, 2);
        end
        M = M(:,ind(1:initial_dims));
        X = X * M;
        clear covX M lambda

        % Run the fast diffusion SNE implementation
        write_data(X, theta, perplexity);
        bh_tsne='';
        %Comment the two lines following this comment and the 
        %corresponding if...else statements if you dont wanna use the 
        %64 bit versions in case they are giving you trouble.
        %In both windows and linux the 64 bit version is abt 2 seconds
        %faster than the corresponding 32 bit version for 6000 datapoints

        arch_str=computer('arch');
        arch=arch_str(length(arch_str)-1:length(arch_str));
        if ismac==1
            bh_tsne='bh_tsne_mac64';
        elseif isunix==1
            if str2double(arch)==64
                bh_tsne='bh_tsne_linux64';
            else
                bh_tsne='bh_tsne_linux32';
            end
        elseif ispc==1
            if str2double(arch)==64
                bh_tsne='bh_tsne_win64';
            else        
                bh_tsne='bh_tsne_win32';
            end
        end
        tic, system([curr_path bh_tsne]); 
        toc
        [mappedX, landmarks, costs] = read_data;   
        landmarks = landmarks + 1;              % correct for Matlab indexing
        delete('data.dat');
        delete('result.dat');

        % while the hash map is too big removing the first element
        % TODO, implement LRU cachce
        if length(mapMat)>5 
            lstKeys= keys(mapMat); 
            remove(mapMat,lstKeys(1));
        end

        % Save results to cache file
        mapMat(hashMat)=mappedX; 
        save('tsneResults.mat','mapMat');
        
        cd(work_dir);
    catch ME
        cd(work_dir);
        rethrow(ME);
    end
end


% Writes the datafile for the fast t-SNE implementation
function write_data(X, theta, perplexity)
    [n, d] = size(X);
    h = fopen('data.dat', 'wb');
	fwrite(h, n, 'integer*4');
	fwrite(h, d, 'integer*4');
    fwrite(h, theta, 'double');
    fwrite(h, perplexity, 'double');
	fwrite(h, X', 'double');
	fclose(h);
end


% Reads the result file from the fast t-SNE implementation
function [X, landmarks, costs] = read_data
    h = fopen('result.dat', 'rb');
	n = fread(h, 1, 'integer*4');
	d = fread(h, 1, 'integer*4');
	X = fread(h, n * d, 'double');
    landmarks = fread(h, n, 'integer*4');
    landmarks = landmarks + 1;
    costs = fread(h, n, 'double');      % this vector contains only zeros
    X = reshape(X, [d n])';
	fclose(h);
end

function Graph2Binary( G, filename )

	% Write to binary file
    G = tril(G);
    fprintf(1,'Writing graph to temporary .bin file ... %3d%% ', 0);
    textname = [filename '.bin'];
    out = fopen( textname, 'w+' );
    f = find( G > 0 );
    [i,j] = ind2sub( size(G), f );
    srctarget = uint32([i-1 j-1]);
    weights = G(f);
    tick = round(length(i)/10);
    for idx = 1:length(f)
        fwrite(out,srctarget(idx,:),'uint32');
        fwrite(out,full(weights(idx)),'double');
        if ~mod(idx,tick) || idx==length(srctarget)
            fprintf(1,'\b\b\b\b\b%3d%% ',round(idx/length(srctarget)*100));
        end
    end
    fprintf(1,'\n');
    fclose(out);
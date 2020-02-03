function ref_genesets = read_prepare_ref_genesets
    filename = 'downloads/msigdb.v6.1.symbols.gmt';
    fid = fopen(filename);
    d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
    d(d==char(13))=[];
    ind_line_end = find(d==char(10));
    all_genes = cell(0);
    tic
    fprintf('Parsing each of %d genesets ... %5d ', length(ind_line_end), 0)
    for i=1:length(ind_line_end)
        fprintf('\b\b\b\b\b\b%5d ', i)
        if i==1
            this_line = d(1:ind_line_end(1)-1);
        else
            this_line = d(ind_line_end(i-1)+1:ind_line_end(i)-1);
        end
        tmp = regexp(this_line,char(9),'split');
        geneset_name(i,1) = tmp(1);
        geneset_ref(i,1) = tmp(2);
        geneset_members(i,1) = {tmp(3:end)};
        all_genes = union(all_genes, tmp(3:end));
    end
    toc
    
    MSigDB_gene_symbol = unique(all_genes);
    MSigDB_geneset_name = geneset_name;
    MSigDB_geneset_ref = geneset_ref;
    MSigDB_geneset_membership = sparse(length(MSigDB_gene_symbol), length(MSigDB_geneset_name));
    fprintf('Saving each of %d genesets to sparse matrix format ... %5d ', length(ind_line_end), 0)
    for i=1:length(ind_line_end)
        fprintf('\b\b\b\b\b\b%5d ', i)
        MSigDB_geneset_membership(ismember(MSigDB_gene_symbol, geneset_members{i}),i)=1;
    end
    toc
    
    
    filenames = getfilenames('downloads');
    filenames(isInListFront(filenames,'msigdb.'))=[];
    MSigDB_geneset_category = cell(length(MSigDB_geneset_name),1);
    tic
    fprintf('Attaching category information ... %5d ', 0)
    for file_ind=1:length(filenames)
        fprintf('\b\b\b\b\b\b%5d ', file_ind)
        category_name = filenames{file_ind}(1:strfind(filenames{file_ind},'.v')-1);
        fid = fopen(fullfile('downloads',filenames{file_ind}));
        d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
        d(d==char(13))=[];
        ind_line_end = find(d==char(10));
        geneset_name = cell(0);
        geneset_ref = cell(0);
        for i=1:length(ind_line_end)
            if i==1
                this_line = d(1:ind_line_end(1)-1);
            else
                this_line = d(ind_line_end(i-1)+1:ind_line_end(i)-1);
            end
            tmp = regexp(this_line,char(9),'split');
            geneset_name(i,1) = tmp(1);
            geneset_ref(i,1) = tmp(2);
        end
        [C,IA,IB] = intersect(MSigDB_geneset_name, geneset_name);
        if length(C)~=length(geneset_name)
            error('problem');
        end
        MSigDB_geneset_category(IA) = {category_name};
    end
    toc
    
    tic
    fprintf('Saving processed MSigDB genesets ... ', length(ind_line_end), 0)
    save('MSigDB_processed_data.mat', 'MSigDB_geneset_membership', 'MSigDB_gene_symbol', 'MSigDB_geneset_name', 'MSigDB_geneset_ref', 'MSigDB_geneset_category');
    toc
    
end

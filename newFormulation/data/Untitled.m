% Treating these data

fullProtein = model_adj.genes;
data1 = readcell('GSM1054768_pre_gene_exp.txt');
data2 = readcell('GSM1054769_3_min_gene_exp.txt');
data3 = readcell('GSM1054770_10_min_gene_exp.txt');
data4 = readcell('GSM1054771_30_min_gene_exp.txt');
data5 = readcell('GSM1054772_60_min_gene_exp.txt');

data_new = cell(0,0);

for i = 1:length(fullProtein)
    protein = fullProtein{i};
    
    if ~contains(protein,'Cre')
        data_new{i,1} = 0;
        data_new{i,2} = 0;
        data_new{i,3} = 0;
        data_new{i,4} = 0;
        data_new{i,5} = 0;
        continue;
    end
    
    protein(end-4:end) = [];
    idx = find(strcmp(data1(:,1),protein));
    
    if ~isempty(idx)
        data_new{i,1} = data1{idx,10};
        data_new{i,2} = data2{idx,10};
        data_new{i,3} = data3{idx,10};
        data_new{i,4} = data4{idx,10};
        data_new{i,5} = data5{idx,10};
    else
        data_new{i,1} = 0;
        data_new{i,2} = 0;
        data_new{i,3} = 0;
        data_new{i,4} = 0;
        data_new{i,5} = 0;
    end
end

writecell(data_new,'Transcriptome_new');

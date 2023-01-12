% post treatment of quantified data

if ~exist('fst','var')
    fst = fastaread('nucleotide_nuc.txt');
    fst = struct2table(fst);
end

fullGenes = fst.Header;
rm = [];
for i = 1:length(fullGenes)
    st = split(fullGenes{i},'] [');
    idx = find(contains(st,'db_xref=Phytozome:'));
    
    if ~isempty(idx)
        fullGenes{i} = erase(st{idx},'db_xref=Phytozome:');
    else
        idx = find(contains(st,'locus_tag='));
        fullGenes{i} = erase(st{idx},'locus_tag=');
    end
end
clear idx i rm st;
fullGenes = unique(fullGenes);

MP = zeros(length(fullGenes),height(fst));
for i = 1:length(fullGenes)
    MP(i,find(contains(fst.Header,fullGenes{i}))) = 1;
end

hder = {};
for i = 1:17
    if i < 10
        dt_raw = readcell(['quant/SRR117440',num2str(i),'/quant.sf'],'FileType','text');
        hder{1,end+1} = ['SRR117440',num2str(i)];
    else
        dt_raw = readcell(['quant/SRR11744',num2str(i),'/quant.sf'],'FileType','text');
        hder{1,end+1} = ['SRR11744',num2str(i)];
    end

    if i == 1
        dt = dt_raw(2:end,4);
    else
        dt = [dt,dt_raw(2:end,4)];
    end
end

dt_new = MP * cell2mat(dt);
dt_new = ['Gene',hder;fullGenes,num2cell(dt_new)];
writecell(dt_new,'output_TAG_transcriptome.txt');

% compare newly quant result with processed data from the study
dt_old = readtable('/Users/16hy16/Documents/MATLAB/projects/Algal/case_study/GSE55253_Expression_GEO_2014Feb.txt');;

for i = 1:height(dt_old)
    dt_old.ID{i} = dt_old.ID{i}(1:end-5);
end

commonGenes = intersect(fullGenes,dt_old.ID);
cprTable = zeros(length(commonGenes),2);

for i = 1:length(commonGenes)
    cprTable(i,1) = dt_new{find(strcmp(dt_new(:,1),commonGenes{i})),2};
    cprTable(i,2) = dt_old.UGTWs6Na0aHS1(find(strcmp(dt_old.ID,commonGenes{i})));
%     cprTable(i,1) = dt_new{find(strcmp(dt_new(:,1),commonGenes{i})),3};
%     cprTable(i,2) = dt_old.UGTWs6na0bHS1(find(strcmp(dt_old.ID,commonGenes{i})));
end

figure;
scatter(cprTable(:,1),cprTable(:,2));
xlabel('New TPM');
ylabel('Old TPM');
set(gca,'YScale','log');
set(gca,'XScale','log');

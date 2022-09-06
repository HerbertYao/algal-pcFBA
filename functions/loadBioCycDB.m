function biocyc = loadBioCycDB(pth)
% NEED TO BE REWORK ON
% read raw data compounds, enzrxns, genes, protcplxs, rxns, and
% metabolicRxns sbml

% standarized format:
%   .enzrxns: col 1: unique-id (ENZRXN66-19089)
%             col 2: description in word (proline dehydrogenase)
%             col 3: gene product / protein (HS01958-MONOMER)
%             col 4: rxn unique-id (RXN0-7005)
%             col 5: formula ((S)-1-pyrroline-5-carboxylate + 3 H+ + 2 e-  <-->  L-proline)
%   .genes: col 1: unique-id (HS04248)
%           col 2: gene name (PKD2)
%           col 3: GENE ID THATS USED IN M-MODEL (5311)
%           col 4: gene product / protein (ENSG00000118762-MONOMER)
%   .rxns: col 1: unique-id (RXN-14181)
%          col 2: enzrxns unique-id, separated by comma (ENZRXNIO2-8554,ENZRXNIO2-7365)
%          col 3: reactants (Beta-D-Galactosides,WATER)
%          col 4: products (Non-Galactosylated-Galactose-Acceptors,D-galactopyranose)
%          col 5: reactants stoich coef ([-1,-1])
%          col 6: products stoich coef ([1,1])
%          col 7: formula from metabolic-reactions.xml
%   .protcplx: col 1: unique-id (CPLX4LZ-110)
%              col 2: subunit stoichiometry (3*ENSG00000118762-MONOMER)
%              col 3: name (Complex III)

addpath(pth);

biocyc.compounds_raw = readcell([pth,'/compounds.dat'],'Delimiter',' - ');
biocyc.compounds_raw(find(startsWith(biocyc.compounds_raw(:,1),'#')),:) = [];
biocyc.compounds_raw = biocyc.compounds_raw(:,1:2);

biocyc.enzrxns_raw = readcell([pth,'/enzrxns.dat'],'Delimiter',' - ');
biocyc.enzrxns_raw(find(startsWith(biocyc.enzrxns_raw(:,1),'#')),:) = [];
biocyc.enzrxns_raw = biocyc.enzrxns_raw(:,1:2);

biocyc.genes_raw = readcell([pth,'/genes.dat'],'Delimiter',' - ');
biocyc.genes_raw(find(startsWith(biocyc.genes_raw(:,1),'#')),:) = [];
biocyc.genes_raw = biocyc.genes_raw(:,1:2);

biocyc.protcplxs_raw = readcell([pth,'/protcplxs.col'],'FileType','text','Delimiter','\t');
biocyc.protcplxs_raw(find(startsWith(biocyc.protcplxs_raw(:,1),'#')),:) = [];

biocyc.rxns_raw = readcell([pth,'/reactions.dat'],'Delimiter',' - ');
biocyc.rxns_raw(find(startsWith(biocyc.rxns_raw(:,1),'#')),:) = [];
biocyc.rxns_raw = biocyc.rxns_raw(:,1:2);

biocyc.metabolicRxns = sbmlimport('metabolic-reactions.xml');
% parse out metabolic Reaction names
biocyc.metabolicRxnNames = cell(length(biocyc.metabolicRxns.Reactions),1);
for i = 1:length(biocyc.metabolicRxnNames)
    biocyc.metabolicRxnNames{i} = biocyc.metabolicRxns.Reactions(i).Name;
end

% Process .dat file info
biocyc.compounds = datParse(biocyc.compounds_raw,{'UNIQUE-ID','COMMON-NAME','INCHI','INCHI-KEY'});
biocyc.genes = datParse(biocyc.genes_raw,{'UNIQUE-ID','COMMON-NAME','DBLINKS','PRODUCT'},...
    {'','','NCBI-GENE',''});
biocyc.rxns = datParse(biocyc.rxns_raw,{'UNIQUE-ID','ENZYMATIC-REACTION','LEFT','RIGHT'});
biocyc.enzrxns = datParse(biocyc.enzrxns_raw,{'UNIQUE-ID','COMMON-NAME','ENZYME','REACTION'});

% process the third column of hmcyc.genes to pure ncbi id
for i = 1:length(biocyc.genes)
    if isempty(biocyc.genes{i,3})
        continue;
    end

    str = split(biocyc.genes{i,3},'"');
    biocyc.genes{i,3} = '';

    for j = 1:((length(str)-1)/2)
        biocyc.genes{i,3} = [biocyc.genes{i,3},str{2*j}];
        if j ~= ((length(str)-1)/2)
            biocyc.genes{i,3} = [biocyc.genes{i,3},','];
        end
    end
end

% add reaction stoichiometry to hmcyc.rxns
% first put everything to default: -1 for reagents, +1 for products
biocyc.rxns(:,5:6) = cell(length(biocyc.rxns),2);
for i = 1:length(biocyc.rxns)
    if ~isempty(biocyc.rxns{i,3})
        biocyc.rxns{i,5} = -1*ones(1,length(split(biocyc.rxns{i,3},',')));
    end
    if ~isempty(biocyc.rxns{i,4})
        biocyc.rxns{i,6} = ones(1,length(split(biocyc.rxns{i,4},',')));
    end
end

% search for coef other than 1 in rxns_raw
rxnIdx = 1;
for i = 1:length(biocyc.rxns_raw)
    if strcmp(biocyc.rxns_raw{i,1},'//')
        rxnIdx = rxnIdx + 1;

    elseif strcmp(biocyc.rxns_raw{i,1},'^COEFFICIENT')
        if strcmp(biocyc.rxns_raw{i-1,1},'LEFT')
            idx = find(strcmp(split(biocyc.rxns{rxnIdx,3},','),biocyc.rxns_raw{i-1,2}));
            try
                biocyc.rxns{rxnIdx,5}(idx) = -biocyc.rxns_raw{i,2};
            catch
                warning('Unknown stoich assignment error: line %d, rxn %d\n',i,rxnIdx);
            end
        elseif strcmp(biocyc.rxns_raw{i-1,1},'RIGHT')
            idx = find(strcmp(split(biocyc.rxns{rxnIdx,4},','),biocyc.rxns_raw{i-1,2}));
            try
                biocyc.rxns{rxnIdx,6}(idx) = biocyc.rxns_raw{i,2};
            catch
                warning('Unknown stoich assignment error: line %d, rxn %d\n',i,rxnIdx);
            end
        end
    end
end

% add reaction formula from biocyc.metabolicRxns to biocyc.rxns 7th col
biocyc.rxns(:,7) = cell(length(biocyc.rxns),1);

for i = 1:length(biocyc.rxns)
    biocyc.rxns{i,7} = '';
    if isempty(biocyc.rxns{i,1})
        continue;
    end

    idx = find(strcmp(biocyc.metabolicRxnNames,biocyc.rxns{i,1}));
    if isempty(idx)
        idx = find(startsWith(biocyc.metabolicRxnNames,biocyc.rxns{i,1}));
    end

    if isempty(idx)
        continue;
    end

    biocyc.rxns{i,7} = char(biocyc.metabolicRxns.Reactions(idx(1),1).Reaction);
end

% Constructing humancyc 'S matrix'
biocyc.S = zeros(length(biocyc.compounds),length(biocyc.rxns));
for i = 1:length(biocyc.rxns)

    if ~isempty(biocyc.rxns{i,3})
%       find reagents from hm.rxns
        reag = split(biocyc.rxns{i,3},',');
%       assign to hmcyc.S
        for j = 1:length(reag)
            idx = find(strcmp(biocyc.compounds(:,1),reag{j}));
            biocyc.S(idx,i) = biocyc.rxns{i,5}(j);
        end
    end

%   same for products
    if ~isempty(biocyc.rxns{i,4})
        prod = split(biocyc.rxns{i,4},',');
        for j = 1:length(prod)
            idx = find(strcmp(biocyc.compounds(:,1),prod{j}));
            biocyc.S(idx,i) = biocyc.rxns{i,6}(j);
        end
    end
end

% Process hmcyc.protcplxs and constructing hmcyc.P matrix, relationship
% between genes and cplx
[~,w] = size(biocyc.protcplxs_raw);

biocyc.protcplxs = [biocyc.protcplxs_raw(2:end,1),biocyc.protcplxs_raw(2:end,w),biocyc.protcplxs_raw(2:end,2)];
biocyc.P = zeros(length(biocyc.genes),length(biocyc.protcplxs));

for i = 1:length(biocyc.protcplxs)

    if isempty(biocyc.protcplxs{i,2}) || any(ismissing(biocyc.protcplxs{i,2}))
        continue;
    end

%   analysing complex stoich from string
    su = split(biocyc.protcplxs{i,2},',');
    for j = 1:length(su)
        dt = split(su{j},'*');
        idx = find(strcmp(biocyc.genes(:,4),dt{2}));
        if isempty(idx)
            continue;
        elseif length(idx) ~= 1
            warning('Duplication in hmcyc.genes col 5: row %d and row %d',idx(1),idx(2));
            biocyc.P(idx(1),i) = str2double(dt{1});
        else
            biocyc.P(idx,i) = str2double(dt{1});
        end
    end
end

% Constructing hmcyc.C matrix, relationship between genes and enzrxns
biocyc.C = zeros(length(biocyc.genes),length(biocyc.enzrxns));

for i = 1:length(biocyc.enzrxns)
    if isempty(biocyc.enzrxns{i,3})
        continue;
    end

%   if enzrxn col 3 is directly found as a gene, record and continue
%   elseif find it as a complex, copy the column from hmcyc.P
%   else return error
    idx = find(strcmp(biocyc.genes(:,4),biocyc.enzrxns{i,3}));
    
    if ~isempty(idx)
        if length(idx) ~= 1
            warning('Duplication in hmcyc.genes col 5: row %d and row %d',idx(1),idx(2));
        end
        biocyc.C(idx(1),i) = 1;
    else
        idx = find(strcmp(biocyc.protcplxs(:,1),biocyc.enzrxns{i,3}));
        if isempty(idx)
            warning('Unrecognized protein or protcplx: %s',biocyc.enzrxns{i,3});
            continue;
        end

        biocyc.C(:,i) = biocyc.P(:,idx);
    end
end

% Constructing hmcyc.K matrix, relationship between enzrxns and rxns
biocyc.K = zeros(length(biocyc.rxns),length(biocyc.enzrxns));

for i = 1:length(biocyc.enzrxns)
    if isempty(biocyc.enzrxns{i,4})
        continue;
    end
    idx = find(strcmp(biocyc.rxns(:,1),biocyc.enzrxns{i,4}));
    if isempty(idx)
        warning('Unrecognized rxns mapping from enzrxns: %s',biocyc.enzrxns{i,4});
    end

    biocyc.K(idx,i) = 1;
end

end
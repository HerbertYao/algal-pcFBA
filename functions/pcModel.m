function [constraintedModel,fullProtein,fullCplx,C_matrix,K_matrix,fullProteinMM] = pcModel(model,fasta,maxWeightFrac)

% (MAIN) A function that produce a proteomically constrainted m-model from 
% a common m-model and protein sequence data. 
% 
% USAGE:
%   new_model = proteinConstraintModel(model,'Paeruginosa.fasta',550);
% 
% INPUTS:
%   model: A functional COBRA model struct with the field 'genes'
%   fasta: The file name of a proteome fasta file containing geneID 
%          (same set of ID as model.genes) as header and respective protein 
%          sequence in one-letter symbol.
%          Alternatively, you can put an estimated average protein length 
%          (avg # of amino acid per protein) if you don't have access to
%          fasta. The modelling accuracy will be compromised.
% 
% OPTIONAL INPUTS:
%   maxWeightFrac: A double denotes the maximum weight fraction of total
%                  protein components, in mg/gDW. Default = 550mg/gDW
%   keff_refTable: A n*2 cell array with enzyme ID in the first column and
%                  keff value in the second column, unit in 1/h. Enzymes 
%                  not found in this table will be assigned keff = 234000
%                  (NOT YET FINISHED)
% 
% OUTPUTS:
%   constraintedModel: A m-model with proteomic constraint
%   fullProtein:       A full list of proteins added to model.mets
%   fullEnzyme:        A full list of enzymes added to model.mets
%   C_matrix:          A matrix which related to enzyme formation
%   K_matrix:          A matrix of enzyme's catalyzing coefficient (k_eff)
%   fullProteinMM:     Molar mass of each protein in fullProtein
% 
% REQUIREMENTS:
%   MATLAB with Bioinformatics and Deep Learning Toolbox installed
%   Configured CobraToolbox with at least one capable solver
% 
% Important Note:
%   New_model.mets will contain all proteins and enzymes (protein
%   complexes). Each protein and enzyme will have an exchange rxn and a
%   enzyme formation rxn. Fluxes of these reactions are not real flux, but
%   representing their concentration in nmol/gDW. Therefore, the
%   concentration of protein_i is -v_{EX_protein_i} and the concentration 
%   of enzyme_j is v_{EX_enzyme_j} = v_{enzymeForm_enzyme_j}

% -------------------------------------------------------------------------
%% Step 0: Parse Inputs

if isa(fasta,'char')
    fastaFile = fastaread(fasta);
    fprintf('Fasta file read successfully\n');
    
elseif isa(fasta,'double')
    avgProteinLen = fasta;
    fprintf('Input protein length: %d\n',avgProteinLen);
    
else
    error('Second input is neither an existing filename nor avg length\n');
end

if ~exist('maxWeightFrac','var')
    maxWeightFrac = 150;
    fprintf('Using default maximum protein weight fraction of 0.15\n');
    
elseif maxWeightFrac > 500
    warning('Maximum proteome weight fraction is over 50%.\n');
    
elseif maxWeightFrac < 5
    warning('Maximum proteome weight fraction is under 5%.\n');
end

%% Step 1: List all unique enzymes and construct k_eff matrix
fprintf('Constructing k_eff matrix...');

% 1.1 List unique cplxes and put in the form of x(655)x(663)x(659)
fullCplx = {};

for i = 1:length(model.rules)
    if ~isempty(model.rules{i})
        try
            cplxList = split(parseGeneRule(model.rules{i}),';');
        catch
            error('Error in geneRule parsing: rxn %d\n',i);
        end
        
%       Add new enzymes to fullEnzyme
        for j = 1:length(cplxList)
            if ~any(strcmp(fullCplx,cplxList{j}))
                fullCplx{length(fullCplx)+1,1} = cplxList{j};
            end
        end
    end
end

% 1.2 Construct k_matrix such that vj <= sum_i (keff_ij*f)*ei
fullRxn = model.rxns;
K_matrix = zeros(length(fullRxn),length(fullCplx));

for i = 1:length(model.rules)
    cplxList = split(parseGeneRule(model.rules{i}),';');
    
    for j = 1:length(cplxList)
        idx = find(strcmp(fullCplx,cplxList{j}));
        K_matrix(i,idx) = 234000;
    end
end
fprintf('done\n')

%% Step 2: Parse the fullCplx to construct C matrix
fprintf('Parsing geneRules and constructing C matrix...');

% Initialize C_matrix
fullProtein = model.genes;
C_matrix = zeros(length(fullProtein),length(fullCplx));

% Parse and assign values
for i = 1:length(fullCplx)
    
%   Remove brackets and the first 'x'
    cplx = erase(fullCplx{i},{'(',')'});
    cplx(1) = [];
    
%   Assign entries in C_matrix
    cplxComp = split(cplx,'x');
    for j = 1:length(cplxComp)
        C_matrix(str2double(cplxComp{j}),i) = 1;
    end
end

fprintf('done\n');

%% Step 3: Add proteins and enzymes to model.mets
%         Add EX_protein, EX_enzyme, and enzymeForm to model.rxns

% 3.1 Add protein and EX_protein first
fprintf('Adding protein metabolites and protein exchange reactions...');

for i = 1:length(fullProtein)
    name = ['protein_',fullProtein{i}];
    if ~any(strcmp(model.mets,name))
        model = addMetabolite(model,name,'csense','G');
        model = addExchangeRxn(model,name,-1000000,0);
    else
        warning('Replicate proteins in model: %s\n',name);
        model = addMetabolite(model,[name,'_1'],'csense','G');
        model = addExchangeRxn(model,[name,'_1'],-1000000,0);
    end
end

fprintf('done\n');

% 3.2 Add cplx and cplxForm reactions
fprintf('Adding complex metabolites and cplx formation reactions...');

for i = 1:length(fullCplx)
%   Metabolite
    cplxName = ['cplx_',fullCplx{i}];
    model = addMetabolite(model,cplxName);
    
%   cplxForm rxns
    proteinIdx = find(C_matrix(:,i));
    
    rxnList = fullProtein(proteinIdx);
    for j = 1:length(rxnList)
        rxnList{j} = ['protein_',rxnList{j}];
    end
    coefList = - C_matrix(proteinIdx,i);
    
    rxnList{length(rxnList)+1,1} = cplxName;
    coefList(length(coefList)+1,1) = 1;
    
    model = addReaction(model,['cplxForm_',fullCplx{i}],...
        'metaboliteList',rxnList,...
        'stoichCoeffList',coefList,...
        'reversible',false);
end

fprintf('done\n');

% 3.3 Add enzyme, EX_enzyme, and enzymeForm reactions
fprintf('Adding enzymes metabolites, enzyme formation, and enzyme exchange reactions...');

for i = 1:length(fullRxn)
    
%   Pass if the rxn doesn't need enzymes
    if all(K_matrix(i,:)==0)
        continue;
    end
    
%   Add 2 enzymes for each reaction, 1 forward and 1 reverse
    enzName = ['enzyme_',fullRxn{i}];
    enzNameRev = [enzName,'_rev'];
    model = addMetabolite(model,enzName);
    model = addMetabolite(model,enzNameRev);
    
%   Find all protein complex that can run this reaction
    cplxList = find(K_matrix(i,:));
    
%   Mapping protein complex to this enzymes
    for j = 1:length(cplxList)
        cplxName = ['cplx_',fullCplx{cplxList(j)}];
        
%       Adding forward enzymeForm reaction
        model = addReaction(model,['enzymeForm_',fullRxn{i},'_',cplxName],...
            'metaboliteList',{cplxName,enzName},...
            'stoichCoeffList',[-1,1],...
            'reversible',false);
        
%       Adding reverse enzymeForm reaction
        model = addReaction(model,['enzymeForm_',fullRxn{i},'_rev_',cplxName],...
            'metaboliteList',{cplxName,enzNameRev},...
            'stoichCoeffList',[-1,1],...
            'reversible',false);
    end
end

% Add all enzyme exchange reactions
for i = 1:length(model.mets)
    if contains(model.mets{i},'enzyme_')
        model = addExchangeRxn(model,model.mets{i},0,1000000);
    end
end

fprintf('done\n');

%% Step 4: Add proteomic constraints to the model using K_matrix
fprintf('Adding enzymatic constraints to the model...');

f = 1/1000000; % conversion factor from nmol to mmol

for i = 1:length(fullRxn)
    
%   Pass if the respective K_matrix row is all 0
    if all(K_matrix(i,:) == 0)
        
        if ~isempty(model.rules{i}) % safeguard
            warning('K matrix construction is wrong for rxn %s\n',model.rxns{i});
        end
        continue;
        
%   Else add forward and reverse enzymes to rxn
    else
        if isempty(model.rules{i}) % safeguard
            warning('K matrix construction is wrong for rxn %s\n',model.rxns{i});
        end
        
%       Forward enzyme
        idx = find(strcmp(model.mets,['enzyme_',fullRxn{i}]));
        model.S(idx,i) = - 1/234000/f;
%       Reverse enzyme
        idx = find(strcmp(model.mets,['enzyme_',fullRxn{i},'_rev']));
        model.S(idx,i) = 1/234000/f;
        
    end
end

fprintf('done\n');

%% Step 5: Add total proteome weight constraint to the model
fprintf('Adding protein weight constraint to the model...');

% 5.1 Calculate each protein's molar mass either using avg length or fasta
fullProteinMM = zeros(length(fullProtein),1);

% If fasta is available, find the sequence first and calc molar mass
if exist('fastaFile','var')
    
%   Calculate the avg protein length
    avgLen = 0;
    for i = 1:length(fastaFile)
        avgLen = avgLen + length(fastaFile(i).Sequence);
    end
    avgLen = avgLen / length(fastaFile);
    
%   Find and calculate each protein's molar mass
    for i = 1:length(fullProteinMM)
        fullProteinMM(i) = calcProteinMM(findProteinSeq(fullProtein{i},fastaFile,avgLen));
    end
    
% If avgProteinLen is input by user, calculate the avgMM and assign to 
% every entry
elseif exist('avgProteinLen','var')
    avgSeq(1:avgProteinLen) = 'Z';
    avgMM = calcProteinMM(avgSeq);
    
    for i = 1:length(fullProteinMM)
        fullProteinMM(i) = avgMM;
    end
    
%   Else return a error message
else
    error('Neither avg protein length or fasta file is found\n');
end

% 5.2 Implement the protein weight constraint

model = addMetabolite(model,'proteinWC','b',maxWeightFrac,'csense','L');

for i = 1:length(fullProtein)
    idx = find(strcmp(model.rxns,['EX_protein_',fullProtein{i}]));
    model.S(length(model.mets),idx) = -fullProteinMM(i)*f;
end

fprintf('done\n');

constraintedModel = model;

end

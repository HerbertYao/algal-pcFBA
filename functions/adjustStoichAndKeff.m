function model_adj = adjustStoichAndKeff(model,C_matrix,K_matrix)

% A function that modify PC model based on C and Keff matrix provided
% 
% USAGE:
%   model_adj = proteinConstraintModel(model,C_matrix,K_matrix);
% 
% INPUTS:
%   model:    A PC model built from the function 'pcModel.m'
%   C_matrix: Denotes stoich relationship between protein and cplx, must be
%             the same dimension as C_matrix output from 'pcModel.m'
%   K_matrix: Denotes catalyzing capability of enzymes to rxns, must be the
%             same demension as K_matrix output from 'pcModel.m'
% 
% OUTPUTS:
%   model_adj: The adjusted PC model

% Checking input format
[ch,cw] = size(C_matrix);
[kh,kw] = size(K_matrix);

proteinIdxList = find(contains(model.mets,'protein_'));

cplxIdxList = find(contains(model.mets,'cplx_'));
cplxList = model.mets(cplxIdxList);

lastRxnIdx = find(strcmp(model.rxns,['EX_',model.mets{proteinIdxList(1)}]))-1; % The last real rxn located before the first EX_protein
rxnIdxList = 1:lastRxnIdx;
rxnList = model.rxns(rxnIdxList);

if cw ~= kw
    error('Width of C_matrix and K_matrix are different, but they must be the same');
elseif ch ~= length(proteinIdxList)
    error('Height of C_matrix must equal to the total number of protein in the model');
elseif kw ~= length(cplxList)
    error('Width of K_matrix must equal to the total number of cplx in the model');
elseif kh ~= length(rxnIdxList)
    error('Height of K_matrix must equal to the total number of rxn in the model');
end

% Adjust cplxForm rxns...
fprintf('Updating complex stoich...');

for i = 1:length(proteinIdxList)
    for j = 1:length(cplxIdxList)
        
        idx = find(strcmp(model.rxns,['cplxForm_',erase(cplxList{j},'cplx_')])); % idx of enzymeForm rxn
        
%       Update entries for all C ~= 0
        if C_matrix(i,j) ~= 0

%           A couple of safeguard
            if C_matrix(i,j) < 0 % Warning for negative C entry but proceed
                warning('Negative stoich coefficient at protein %d, cplx %d\n',i,j);
            end
            if model.S(proteinIdxList(i),idx) == 0 % Warning for adding extra specie of protein to cplx
                warning('Coef of protein %d, cplx %d changed from zero to a non-zero value\n',i,j);
            end
            
            model.S(proteinIdxList(i),idx) = - C_matrix(i,j);
            
%       Another warning for changing non-zero coef to zero
        else
            if model.S(proteinIdxList(i),idx) ~= 0
                warning('Coef of protein %d, cplx %d changed from a non-zero value to zero\n',i,j);
                model.S(proteinIdxList(i),idx) = C_matrix(i,j);
            end
        end

    end
    fprintf('%d out of %d proteins\n',i,length(proteinIdxList));
end

fprintf('done\n');

% Update K_eff by manipulating enzymeForm coefficient (default = 1)
%   basically means converting them to an avg enzyme of Keff = 234000
fprintf('Updating rate constants...');

for i = 1:length(cplxIdxList)
    for j = 1:length(rxnIdxList)
        
        if K_matrix(j,i) ~= 0 
            
%           Finding the respective enzymeForm rxn in model
            idx = find(strcmp(model.rxns,['enzymeForm_',rxnList{j},'_',cplxList{i}]));
            
            if length(idx) == 1 % easiest case: enzymeForm exists
                model.S(cplxIdxList(i),idx) = - 234000/K_matrix(j,i);
                model.S(cplxIdxList(i),idx+1) = - 234000/K_matrix(j,i); % Reverse enzyme
                
            elseif isempty(idx) % doesn't yet support adding new keff entries, but will be later
                error('Keff of cplx %d, rxn %d changed from zero to a non-zero value',i,j);
            else
                error('Something is going wrong when modifying keff: cplx %d, rxn %d',i,j);
            end
            
%       Not yet support removing keff entries, but will be later
        else
            idx = find(strcmp(model.rxns,['enzymeForm_',rxnList{j},'_',cplxList{i}]));
            
            if ~isempty(idx)
                error('Keff of cplx %d, rxn %d changed from a non-zero value to zero',i,j);
            end
        end
    end
    
    fprintf('%d out of %d complex\n',i,length(cplxIdxList));
end

fprintf('done\n');

model_adj = model;

end

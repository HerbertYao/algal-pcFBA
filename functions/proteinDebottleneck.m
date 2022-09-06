function [FBAsol,model_new,relaxProt,relaxLevel] = proteinDebottleneck(model,ptBudget)

% return a FBA solution with protein lb relief slightly to optimize the
% objective function.

model_new = model;
% add a new protein budget
model = addMetabolite(model,'newBudget');
model.b(end) = ptBudget;
model.csense(end) = 'L';

% construct a new LP
lp.A = model.S;
lp.b = model.b;
lp.c = model.c;
lp.lb = model.lb;
lp.ub = model.ub;
lp.osense = -1;
lp.csense = model.csense;

% add a set of new protein exchange reactions under newBudget
proteinExIdx = find(startsWith(model.rxns,'EX_protein_'));

for i = 1:length(proteinExIdx)
    lp.A = [lp.A,lp.A(:,proteinExIdx(i))];
    lp.c(end+1) = 0;
    lp.lb(end+1) = -1000;
    lp.ub(end+1) = 0;
    lp.A(length(lp.b),length(lp.c)) = -1; % new EX_prot consumes ptBudget
end

proteinExIdxNew = length(model.rxns)+1:length(lp.c); % indexing new EX_prot

% solve and return new model
lpSol = solveCobraLP(lp);

relaxProt = {};
relaxLevel = [];

for i = 1:length(proteinExIdxNew)
    if lpSol.full(proteinExIdxNew(i)) ~= 0
        relaxProt{end+1,1} = erase(model.rxns{proteinExIdx(i)},'EX_');
        relaxLevel(end+1,1) = - lpSol.full(proteinExIdxNew(i));
        model_new.lb(proteinExIdx(i)) = model_new.lb(proteinExIdx(i)) + lpSol.full(proteinExIdxNew(i));
        model_new.ub(proteinExIdx(i)) = model_new.ub(proteinExIdx(i)) + lpSol.full(proteinExIdxNew(i));
    end
end

FBAsol = optimizeCbModel(model_new,'max');

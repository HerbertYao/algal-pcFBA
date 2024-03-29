% New Keff Estimation Method
%  Must have Gurobi solver installed

if ~exist('model_ori','var')
    load('pso_input.mat');
    model_ori = model_adj;
end

model_adj = model_ori;
model_adj.lb(1799) = 0.02;

%% Treat and rescale the data

% Count the percent of modelled proteome in the total proteome
%   1. calculate MW for protein with available sequence
%   2. calculate m and u by sum_ T*MW
%   3. budget = m / (m + u) * 550mg/gDW

fasta = struct2cell(fastaread('all_fasta.fasta'))';
rnaseq = readtable('rnaseq.txt');
data_new = readcell('Transcriptome_new.txt');
data_new = cell2mat(data_new);

if ~exist('fastaMW','var')
    fastaMW = zeros(length(fasta),1);
    for i = 1:length(fastaMW)
        fastaMW(i) = calcProteinMM(fasta{i,2})*1e-6;
    end
end

m = 0;
u = 0;
ct = 0;

for i = 1:height(rnaseq)
    idx = find(contains(fasta(:,1),rnaseq.geneName{i}));
    bool = any(contains(model_adj.rxns,rnaseq.geneName{i}));

    if bool
        ct = ct + 1;
        if ~isempty(idx)
            m = m + fastaMW(idx(1))*rnaseq.piReplete(i);
        else
            m = m + mean(fastaMW)*rnaseq.piReplete(i);
        end
    
    else
        if ~isempty(idx)
            u = u + fastaMW(idx(1))*rnaseq.piReplete(i);
        else
            u = u + mean(fastaMW)*rnaseq.piReplete(i);
        end
    end
end

% One more step: scale m based on the fraction of modelled protein not in
% rnaseq data
m = m / ct * length(find(contains(model_adj.rxns,'EX_protein_')));

model_adj.b(8772) = 550 * m / (m + u);

clear i idx ct bool m u;

% Scale data based on protein mass budget
T_rep = cell2mat(overlaidT_rep(:,2));
T_dep = cell2mat(overlaidT_dep(:,2));

proteinExIdx = find(startsWith(model_adj.rxns,'EX_protein_'));
bgt = model_adj.b(find(strcmp(model_adj.mets,'proteinWC'))); % total proteome mass budget in mg
mw = -model_adj.S(find(strcmp(model_adj.mets,'proteinWC')),proteinExIdx); % protein MW vector in mg/nmol

sc = bgt / (mw * T_rep);
T_rep = T_rep * sc;
sc = bgt / (mw * T_dep);
T_dep = T_dep * sc;

for i = 1:5
    sc = bgt / (mw * data_new(:,i));
    data_new(:,i) = data_new(:,i) * sc;
end

%% Formulate base model
model.A = model_adj.S;
model.rhs = model_adj.b;
model.lb = model_adj.lb;
model.ub = model_adj.ub;
model.varnames = model_adj.rxns;
model.constrnames = model_adj.mets;
model.sense = '';

for i = 1:length(model_adj.csense)
    if strcmp(model_adj.csense(i),'E')
        model.sense(i,1) = '=';
    elseif strcmp(model_adj.csense(i),'L')
        model.sense(i,1) = '<';
    elseif strcmp(model_adj.csense(i),'G')
        model.sense(i,1) = '>';
    else
        error('Wrong csense in the model');
    end
end

% Relief the equality constraint on protein-cplx stoich so that un-used
% proteins are allowed
for i = 1:length(model.sense)
    if startsWith(model.constrnames{i},'protein_')
        model.sense(i) = '>';
    end
end

% Break cplx - enz reactions, have separate dilute reactions instead
cplxDLIdx = find(startsWith(model_adj.rxns,'enzymeForm_'));
enzDMIdx = [length(model_adj.rxns)+1:length(model_adj.rxns)+length(cplxDLIdx)]';
proteinExIdx = find(startsWith(model.varnames,'EX_protein_'));

% Change A matrix
for i = 1:length(enzDMIdx)
    e = find(model_adj.S(:,cplxDLIdx(i))>0); % find the enzyme that's produced
    
%   Modify original enzymeForm reactions to cplx dilution
    model.A(e,cplxDLIdx(i)) = 0;
    model.varnames{cplxDLIdx(i)} = ['DL_',erase(model.varnames{cplxDLIdx(i)},'enzymeForm_')]; % change rxn name to DL_
    
%   Add new enzyme DM reactions
    model.A(e,enzDMIdx(i)) = 1;
    model.varnames{enzDMIdx(i)} = ['DM_',erase(model.varnames{cplxDLIdx(i)},'DL_')]; % name new rxn to DM_
    model.lb(enzDMIdx(i)) = 0;
    model.ub(enzDMIdx(i)) = 1000000;
end

clear e bgt sc;

%% Stack problems on top of the original
%   add variables and constraints

varlen = length(model.varnames);
constrlen = length(model.constrnames);

for i = 1:length(model.varnames)
    model.varnames{end+1} = [model.varnames{i},'__2'];
end
for i = 1:length(model.constrnames)
    model.constrnames{end+1} = [model.constrnames{i},'__2'];
end

for i = 1:length(model.varnames)
    model.varnames{end+1} = [model.varnames{i},'__3'];
end
for i = 1:length(model.constrnames)
    model.constrnames{end+1} = [model.constrnames{i},'__3'];
end

for i = 1:length(model.varnames)
    model.varnames{end+1} = [model.varnames{i},'__4'];
end
for i = 1:length(model.constrnames)
    model.constrnames{end+1} = [model.constrnames{i},'__4'];
end

for i = 1:length(model.varnames)
    model.varnames{end+1} = [model.varnames{i},'__5'];
end
for i = 1:length(model.constrnames)
    model.constrnames{end+1} = [model.constrnames{i},'__5'];
end

modelA = model.A;
model.A = [model.A,sparse(constrlen,varlen);sparse(constrlen,varlen),modelA];
model.A = [model.A,sparse(2*constrlen,varlen);sparse(constrlen,2*varlen),modelA];
model.A = [model.A,sparse(3*constrlen,varlen);sparse(constrlen,3*varlen),modelA];
model.A = [model.A,sparse(4*constrlen,varlen);sparse(constrlen,4*varlen),modelA];

model.rhs = [model.rhs;model.rhs;model.rhs;model.rhs;model.rhs];
model.lb = [model.lb;model.lb;model.lb;model.lb;model.lb];
model.ub = [model.ub;model.ub;model.ub;model.ub;model.ub];
model.sense = [model.sense;model.sense;model.sense;model.sense;model.sense];

% record new idx
cplxDLIdx(:,2) = cplxDLIdx(:,1) + varlen;
enzDMIdx(:,2) = enzDMIdx(:,1) + varlen;
cplxDLIdx(:,3) = cplxDLIdx(:,2) + varlen;
enzDMIdx(:,3) = enzDMIdx(:,2) + varlen;
cplxDLIdx(:,4) = cplxDLIdx(:,3) + varlen;
enzDMIdx(:,4) = enzDMIdx(:,3) + varlen;
cplxDLIdx(:,5) = cplxDLIdx(:,4) + varlen;
enzDMIdx(:,5) = enzDMIdx(:,4) + varlen;

%% New Variables: r (ri = enzj / cplxi, any j)
rIdx = [length(model.varnames)+1:length(model.varnames)+length(fullCplx)]';

for i = 1:length(rIdx)
    model.A(1,rIdx(i)) = 0;
    model.varnames{rIdx(i)} = ['R_',fullCplx{i}];
    
%   if the protein is waivered, don't allow change
%   if the protein abundance is too low, don't allow change
%   else r is bounded [0.1, 10] of the original rate constant
    if variedKeff(i)
        model.lb(rIdx(i)) = 0.1;
        model.ub(rIdx(i)) = 10;
    else
        model.lb(rIdx(i)) = 1;
        model.ub(rIdx(i)) = 1;
    end
end

% New linear constraint: mean(r) = 1
model.rhs(end+1) = length(rIdx);

for i = 1:length(rIdx)
    model.A(length(model.rhs),rIdx(i)) = 1;
end

model.sense(end+1) = '=';
model.constrnames{end+1} = 'r_avg';

%% Quadratic constraints: c*r - e = 0

% there's one constraint per original enzymeForm rxn
%   they have identical index to cplxDLIdx

for j = 1:5
    for i = 1:length(cplxDLIdx(:,j))
    %   Look up which r is corresponding to this enzyme
        cplxName = model.constrnames{find(model.A(:,cplxDLIdx(i,1)))};
        idx = find(strcmp(fullCplx,erase(cplxName,'cplx_')));

        model.quadcon(i+(j-1)*length(cplxDLIdx)).Qc = sparse(length(model.varnames),length(model.varnames));
        model.quadcon(i+(j-1)*length(cplxDLIdx)).Qc(cplxDLIdx(i,j),rIdx(idx)) = 1;
        model.quadcon(i+(j-1)*length(cplxDLIdx)).q = sparse(length(model.varnames),1);
        model.quadcon(i+(j-1)*length(cplxDLIdx)).q(enzDMIdx(i,j)) = -1;
        model.quadcon(i+(j-1)*length(cplxDLIdx)).rhs = 0;
        model.quadcon(i+(j-1)*length(cplxDLIdx)).sense = '=';
    end
end

clear cplxName idx;

% Define objective function min p^2 - 2T*p
model.obj = zeros(length(model.varnames),1);
model.Q = sparse(length(model.varnames),length(model.varnames));

for i = 1:length(proteinExIdx)
%   exclude the waiver list
    if any(strcmp(waiverList,erase(model.varnames{proteinExIdx(i)},'EX_protein_')))
        continue;
    end
    
%   EX_p are non-positive fluxes so T is assigned positive (to yield -2Tp)
    model.obj(proteinExIdx(i)) = 2 * T_rep(i);
    model.Q(proteinExIdx(i),proteinExIdx(i)) = 1;
    
    model.obj(proteinExIdx(i)+varlen) = 2* T_dep(i);
    model.Q(proteinExIdx(i)+varlen,proteinExIdx(i)+varlen) = 1;
    
    model.obj(proteinExIdx(i)+varlen*2) = 2* data_new(i,1);
    model.Q(proteinExIdx(i)+varlen*2,proteinExIdx(i)+varlen*2) = 1;
    
    model.obj(proteinExIdx(i)+varlen*3) = 2* data_new(i,2);
    model.Q(proteinExIdx(i)+varlen*3,proteinExIdx(i)+varlen*3) = 1;

    model.obj(proteinExIdx(i)+varlen*4) = 2* data_new(i,3);
    model.Q(proteinExIdx(i)+varlen*4,proteinExIdx(i)+varlen*4) = 1;
end

%% Solve
model.lb(27281) = -0.1;
model.lb(28855) = 0.01;

params.NonConvex = 2;
params.PoolGap = 0.05;
model.modelsense = 'min';
model_gurobi = rmfield(model,{'varnames','constrnames'});
sol = gurobi(model_gurobi,params);

%% Evaluation

r_vector = sol.x(rIdx);
cplxIdx = find(startsWith(model_adj.mets,'cplx_'));

% adjust the model using r
for i = 1:length(r_vector)
    
    if r_vector(i) == 1 % pass unchanged cplx
        continue;
    end
    
    % all enzymeForm rxn thats consuming this particular cplx
    idx = find(contains(model_adj.rxns,['_cplx_',fullCplx{i}]));
    for j = 1:length(idx) % large r means higher keff
        model_adj.S(cplxIdx(i),idx(j)) = model_adj.S(cplxIdx(i),idx(j)) / r_vector(i);
    end
end

[FBAsol,model_fit] = overlayMultiomicsData(model_adj,T_rep,100,waiverList);
% FBAsol = optimizeCbModel(model_adj,'max');

% Plot
X = [];
Y = [];
lbl = {};
for i = 1:length(proteinExIdx)
    if ~any(strcmp(waiverList,erase(model_adj.rxns{proteinExIdx(i)},'EX_protein_')))
        Y(length(Y)+1,1) = model.obj(proteinExIdx(i));
        X(length(X)+1,1) = -sol.x(proteinExIdx(i));
%         Y(length(Y)+1,1) = T_rep(i);
%         X(length(X)+1,1) = - FBAsol.v(proteinExIdx(i));
        lbl{length(lbl)+1,1} = model.varnames{proteinExIdx(i)};
    end
end

X1 = [];
Y1 = [];
lbl1 = {};
for i = 1:length(proteinExIdx)
    if ~any(strcmp(waiverList,erase(model_adj.rxns{proteinExIdx(i)},'EX_protein_')))
        Y1(length(Y1)+1,1) = model.obj(proteinExIdx(i)+varlen);
        X1(length(X1)+1,1) = -sol.x(proteinExIdx(i)+varlen);
        lbl1{length(lbl1)+1,1} = model.varnames{proteinExIdx(i)};
    end
end

figure;
subplot(1,2,1);
scatter(X,Y);
xlabel('Best fitted Proteome Abundance');
ylabel('Transcriptome');
subplot(1,2,2);
scatter(X1,Y1);
xlabel('Best fitted Proteome Abundance');
ylabel('Transcriptome');

% Retrieved solution pool... 7 local optimal obj that is no further than 1%
% away from global optimal obj
metRxnIdx = [1:2647]';

figure;
for i = 1:length(sol.pool)
    subplot(2,3,i);
    scatter(sol.x(metRxnIdx),sol.pool(i).xn(metRxnIdx));
    hold on;
%     scatter(sol.x(metRxnIdx+varlen),sol.pool(i).xn(metRxnIdx+varlen));
%     scatter(sol.x(metRxnIdx+varlen*2),sol.pool(i).xn(metRxnIdx+varlen*2));
%     hold off;
    xlabel('Global Optimum Flux');
    ylabel(['Local Optimum ', num2str(i), ' Flux']);
    xline(0);
    yline(0);
end

% TODOs:
%   lit review to find knowledge gaps in C.reinhardtii

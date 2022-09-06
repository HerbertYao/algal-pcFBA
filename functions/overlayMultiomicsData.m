function [FBAsol,model_fit] = overlayMultiomicsData(varargin)

% A function that find a solution which fits a transcriptome count matrix
% as closely as possible. Reaction constraints have to be set up prior to
% using this function (such as minimal growth rate or production)
% 
% USAGE:
%   FBAsol = overlayMultiomicsData(model_pc,count,5,[0.1,1000]);
% 
% INPUTS:
%   model:    A PC-model produced by function pcModel.m or preferably
%               refined by adjustStoichAndKeff.m
%   data:     Transcriptome or proteome data in N*1 matrix
%   thres:    Lower threshold for an abundance count to be considered
%               relevant. Default = 0
%   waiver:   A N*1 cell contains proteinMets that needs to be 
%               waivered from fitting. E.g. {'protein_b0001'}. Default = {}
%   varargin: optional parameter value pairs
%       keffEstimate: If enzymatic rate constants are treated as variables.
%                     Setting to true will lead to nonconvex QP, which
%                     gurobi solver must be installed. Default = false
%       objWeight: Allows the algorithm to minimize weighted sum of 
%                  proteome fitting. Only works when variableKeff = false. 
%                  Default = []
%       rBounds: Allows setting upper and lower bounds for r variables when
%                variableKeff = true. Default = [0.1,1.9]
%       variedKeff: M*1 binary vector indicating if each keff can be varied
%                   or not. Default = []
% 
% OUTPUTS:
%   FBAsol: FBA solution structure of best possible transcriptome fitting

%% Parser
p = inputParser;
addRequired(p,'model',@isstruct);
addRequired(p,'data',@isnumeric);
addOptional(p,'thres',0,@isnumeric);
addOptional(p,'waiver',{},@iscell);

addParameter(p,'keffEstimate',false,@islogical);
addParameter(p,'objWeight',[],@isnumeric);
addParameter(p,'rBounds',[0.1,1.9],@isnumeric);
addParameter(p,'variedKeff',[],@isnumeric);

parse(p,varargin{:});
model = p.Results.model;
data = p.Results.data;
thres = p.Results.thres;
waiver = p.Results.waiver;
variableKeff = p.Results.keffEstimate;
objWeight = p.Results.objWeight;
rBounds = p.Results.rBounds;
variedKeff = p.Results.variedKeff;

%% Default enzymatic keff -- Quadratic optimization with linear constraints
if ~variableKeff
    % initialize QP problem
    qp.A = model.S;
    qp.b = model.b;
    qp.c = zeros(length(model.c),1);
    qp.lb = model.lb;
    qp.ub = model.ub;
    qp.osense = 1;
    qp.csense = model.csense;
    qp.F = zeros(length(model.c),length(model.c));
    
    % Set qp.c and qp.F
    proteinExIdx = find(contains(model.rxns,'EX_protein_'));
    
    % length of data must equal to length of proteinExIdx, else return error
    if length(proteinExIdx) ~= length(data)
        error('Length of data must equal to the number of exchanged proteins in the model');
    end

    % prepare objective weight
    if isempty(objWeight)
        objWeight = ones(length(proteinExIdx),1);
    end

    % Rescale data so it sums to protein budget
    bgt = model.b(find(strcmp(model.mets,'proteinWC'))); % total proteome mass budget (mg)
    mw = -model.S(find(strcmp(model.mets,'proteinWC')),proteinExIdx); % protein MW vector (mg/nmol protein)
    sc = bgt / (mw * data);
    data = data * sc; % so that: bgt / (data'*mw) = 1
    
    for i = 1:length(proteinExIdx)
        
    %   Waiver
        protName = erase(model.rxns{proteinExIdx(i)},'EX_');
        if any(strcmp(waiver,protName))
            continue;
        end
        
    %   Assigning both linear and quadratic objectives
        if data(i) > thres
            qp.c(proteinExIdx(i)) = data(i) * objWeight(i);
            qp.F(proteinExIdx(i),proteinExIdx(i)) = objWeight(i);
        end
    end
    
    % Solve and prepare the return
    FBAsol = solveCobraQP(qp);
    model_fit = qp;

else
%% Make enzymatic rate constants variable: nonconvex bilinear optimization
%   Gurobi must be installed correctly and on path

%   Prepare
    proteinExIdx = find(startsWith(model.rxns,'EX_protein_'));
    fullProteinIdx = find(startsWith(model.mets,'protein_'));
    fullCplxIdx = find(startsWith(model.mets,'cplx_'));
    fullCplx = model.mets(fullCplxIdx);
    cplxFormIdx = find(startsWith(model.rxns,'cplxForm_'));

%   Rescale all data sets
    bgt = model.b(find(strcmp(model.mets,'proteinWC')));
    mw = -model.S(find(strcmp(model.mets,'proteinWC')),proteinExIdx);

    [~,wid] = size(data);
    for i = 1:wid
        data(:,i) = data(:,i) * (bgt/(mw*data(:,i)));
    end

%   Compile a list of keff that will be varied
%   Extract C matrix
    C_matrix = -model.S(fullProteinIdx,cplxFormIdx);
    
%   if a enzyme is either waivered or too low in abundance, its not varied
%   'too low in abundance' is currently defined as:
%       among all subunits that form the enzyme, there isn't any one 
%       that reach the average transcription level
    if isempty(variedKeff) % only do when no variedKeff is supplied
        variedKeff = zeros(length(fullCplx),1);
        mn = mean(mean(data));
    
        for i = 1:length(variedKeff)
            idx = find(C_matrix(:,i));
            
            for j = 1:length(idx)
                if any(data(idx(j),:) > mn)
                    variedKeff(i) = 1;
                    break;
                end
            end
        end
    end

%   Formulate base model
    model_bp.A = model.S;
    model_bp.rhs = model.b;
    model_bp.lb = model.lb;
    model_bp.ub = model.ub;
    model_bp.varnames = model.rxns;
    model_bp.constrnames = model.mets;
    model_bp.sense = '';

    for i = 1:length(model.csense)
        if strcmp(model.csense(i),'E')
            model_bp.sense(i,1) = '=';
        elseif strcmp(model.csense(i),'L')
            model_bp.sense(i,1) = '<';
        elseif strcmp(model.csense(i),'G')
            model_bp.sense(i,1) = '>';
        else
            error('Wrong csense in the model');
        end
    end

%   Relaxing equality constraints on protein - cplx relationship so that
%   proteins can be translated but not folded
    for i = 1:length(model_bp.sense)
        if startsWith(model_bp.constrnames{i},'protein_')
            model_bp.sense(i) = '>';
        end
    end

%   Break cplx - enz reactions, have separate dilute reactions instead
%       cplxDL: cplx -->
%       enzDM: --> enz
    cplxDLIdx = find(startsWith(model.rxns,'enzymeForm_'));
    enzDMIdx = [length(model.rxns)+1:length(model.rxns)+length(cplxDLIdx)]';

%   Change A matrix
    for i = 1:length(enzDMIdx)
        e = find(model.S(:,cplxDLIdx(i))>0); % find the enzyme that's produced
        
%       Modify original enzymeForm reactions to cplx dilution
        model_bp.A(e,cplxDLIdx(i)) = 0;
        model_bp.varnames{cplxDLIdx(i)} = ['DL_',erase(model_bp.varnames{cplxDLIdx(i)},'enzymeForm_')]; % change rxn name to DL_
        
%       Add new enzyme DM reactions
        model_bp.A(e,enzDMIdx(i)) = 1;
        model_bp.varnames{enzDMIdx(i)} = ['DM_',erase(model_bp.varnames{cplxDLIdx(i)},'DL_')]; % name new rxn to DM_
        model_bp.lb(enzDMIdx(i)) = 0;
        model_bp.ub(enzDMIdx(i)) = 1000000;
    end

%   Stack extra model on the top of the base
%       this part is skipped if only 1 set of data is provided
    varlen = length(model_bp.varnames);
    constrlen = length(model_bp.constrnames);

%   variable names and constraint names
    for i = 1:wid-1
        for j = 1:varlen
            model_bp.varnames{end+1} = [model_bp.varnames{j},'__',num2str(i+1)];
        end
        for j = 1:constrlen
            model_bp.constrnames{end+1} = [model_bp.constrnames{j},'__',num2str(i+1)];
        end
    end

%   modify A, rhs, lb, ub, csense
    A = model_bp.A;
    rhs = model_bp.rhs;
    lb = model_bp.lb;
    ub = model_bp.ub;
    sense = model_bp.sense;

    for i = 1:wid-1
        model_bp.A = [model_bp.A,sparse(constrlen*i,varlen);sparse(constrlen,varlen*i),A];
        model_bp.rhs = [model_bp.rhs;rhs];
        model_bp.lb = [model_bp.lb;lb];
        model_bp.ub = [model_bp.ub;ub];
        model_bp.sense = [model_bp.sense;sense];
    end
    
%   record new idx
    for i = 1:wid-1
        cplxDLIdx(:,i+1) = cplxDLIdx(:,i) + varlen;
        enzDMIdx(:,i+1) = enzDMIdx(:,i) + varlen;
    end

%   New variable: r (r_i = enz_j / cplx_i)
%       r is a ratio for enzymatic rate constant, only a function of cplx
    rIdx = [length(model_bp.varnames)+1:length(model_bp.varnames)+length(fullCplx)]';

    for i = 1:length(rIdx)
        model_bp.A(1,rIdx(i)) = 0;
        model_bp.varnames{rIdx(i)} = ['R_',fullCplx{i}];
        
%       variability of each enz depends on variedKeff bool
%       r is bounded [0.1, 10] of the original rate constant if variable
        if variedKeff(i)
            model_bp.lb(rIdx(i)) = rBounds(1);
            model_bp.ub(rIdx(i)) = rBounds(2);
        else
            model_bp.lb(rIdx(i)) = 1;
            model_bp.ub(rIdx(i)) = 1;
        end
    end

%   New linear constraint: mean(r) = 1
%       this is to prevent bulk keff from starting to deviate from sanity
    model_bp.rhs(end+1) = length(rIdx);
    for i = 1:length(rIdx)
        model_bp.A(length(model_bp.rhs),rIdx(i)) = 1;
    end
    
    model_bp.sense(end+1) = '=';
    model_bp.constrnames{end+1} = 'r_avg';

%   Quadratic Constraints: cplx*r - enz = 0
%       there are one per original enzymeForm rxn

    for j = 1:wid
        for i = 1:length(cplxDLIdx(:,j))
        %   Look up which r is corresponding to this enzyme
            cplxName = model_bp.constrnames{find(model_bp.A(:,cplxDLIdx(i,1)))};
            idx = find(strcmp(fullCplx,cplxName));
    
            model_bp.quadcon(i+(j-1)*length(cplxDLIdx)).Qc = sparse(length(model_bp.varnames),length(model_bp.varnames));
            model_bp.quadcon(i+(j-1)*length(cplxDLIdx)).Qc(cplxDLIdx(i,j),rIdx(idx)) = 1;
            model_bp.quadcon(i+(j-1)*length(cplxDLIdx)).q = sparse(length(model_bp.varnames),1);
            model_bp.quadcon(i+(j-1)*length(cplxDLIdx)).q(enzDMIdx(i,j)) = -1;
            model_bp.quadcon(i+(j-1)*length(cplxDLIdx)).rhs = 0;
            model_bp.quadcon(i+(j-1)*length(cplxDLIdx)).sense = '=';
        end
    end

%   Define objective function min p^2 - 2T*p
%       if wid ~= 1, this is defined for all datasets
    model_bp.obj = zeros(length(model_bp.varnames),1);
    model_bp.Q = sparse(length(model_bp.varnames),length(model_bp.varnames));
    
    for i = 1:length(proteinExIdx)
    %   exclude the waiver list
        if any(strcmp(waiver,erase(model_bp.varnames{proteinExIdx(i)},'EX_')))
            continue;
        end
        
    %   EX_p are non-positive fluxes so T is assigned non-negative (to yield -2Tp)
        for j = 1:wid
            model_bp.obj(proteinExIdx(i) + varlen*(j-1)) = 2 * data(i,j);
            model_bp.Q(proteinExIdx(i)+varlen*(j-1),proteinExIdx(i)+varlen*(j-1)) = 1;
        end
    end

%   Solve and return
    params.NonConvex = 2;
    params.PoolGap = 0.05;
    model_bp.modelsense = 'min';
    model_fit = model_bp;
    model_bp = rmfield(model_bp,{'varnames','constrnames'});
    FBAsol = gurobi(model_bp,params);
    
end

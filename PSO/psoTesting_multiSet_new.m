%% Trying out PSO to estimate keff

if ~exist('model_ori','var')
    load('pso_input');
    model_ori = model_adj;
end

model_adj = model_ori;
data = cell2mat(overlaidT_rep(:,2));
data_dep = cell2mat(overlaidT_dep(:,2));

% slightly modify the PC-model to allow unfolded proteins, or Ce <= p
%   this is done by relieving csense_prot >= 0
for i = 1:length(model_adj.mets)
    if startsWith(model_adj.mets{i},'protein_')
        model_adj.csense(i) = 'G';
    end
end

%% Initiate particle locations
% Take 50? random samples first
noPart = 50;

totalT = 0;

for i = 1:length(overlaidT_rep)
    totalT = totalT + overlaidT_rep{i,2};
end

% Boolean shows keffs to be perturbed
variedKeff = zeros(length(fullCplx),1);

for i = 1:length(variedKeff)
    idx = find(C_matrix_adj(:,i));
    
    for j = 1:length(idx)
        if overlaidT_rep{idx(j),2} > totalT/2000
            variedKeff(i) = 1;
        end
    end
end

clear idx i j;

% ONLY DO THIS WHEN A gBest OF PREVIOUS RESULT IS LOADED
if exist('record_gBest','var')
    clear init_keffScale;
    floored = zeros(length(fullCplx),1);
    
    for i = 1:length(variedKeff)
        if record_gBest(i) == 1
            variedKeff(i) = 0;
        elseif (record_gBest(i)/ksFloor) < 1.001 % enzyme around ksFloor are considered 'floored'
            variedKeff(i) = 1;
            floored(i) = 1;
        else
            variedKeff(i) = 1;
        end
    end
    fprintf('%d enzymes Keff perturbed\n',length(find(variedKeff))-length(find(floored)));
end

% Record the average rate constants of 143 enzymes to be varied. This avg
% must be kept throughout perturbing
K_vector_base = zeros(length(fullCplx),1);

for i = 1:length(K_vector_base)
    K_vector_base(i) = sum(K_matrix_adj(:,i))/length(find(K_matrix_adj(:,i)));
end

avgKeff = mean(K_vector_base(find(variedKeff)));

% Now there are confusions on the objective function...either bm or SSE
%   ideally optimizing for biomass

model_adj = changeRxnBounds(model_adj,'EX_pi_e',-10,'l');

if ~exist('init_keffScale','var')
    % Record
    init_flux = zeros(length(model_adj.rxns),noPart);
    init_flux_dep = zeros(length(model_adj.rxns),noPart);
    init_keffScale = zeros(length(fullCplx),noPart);
    init_obj = zeros(1,noPart);

    % Sampling...
    for i = 1:noPart
        model_alt = model_adj;

    %   Take an random set of 143 values in [0.25,1.75]
    %       0.25 means keff is reduced by 75%
    %       1.75 means keff is increased by 75%
        variedScale = zeros(length(find(variedKeff)),1);
        for j = 1:length(variedScale)
            variedScale(j) = random('Uniform',0.25,1.75);
        end

    %   Rescale the mean of ksScale to 1
        variedScale = variedScale / mean(variedScale);
        
%       If variable 'floored' exists, fix respective keffScale to the floor
%       and distribute their excess uniformly
        if exist('floored','var')
            idx = 1;
            ex = 0;
            for j = 1:length(fullCplx)
                if floored(j)
                    ex = ex + variedScale(idx) - ksFloor;
                    variedScale(idx) = ksFloor;
                    idx = idx + 1;
                elseif variedKeff(j)
                    idx = idx + 1;
                end
            end
            
            for j = 1:length(variedScale)
                if variedScale(j) ~= ksFloor
                    variedScale(j) = variedScale(j) + ex/(length(variedScale)-length(find(floored)));
                end
            end
        end

        if (mean(variedScale) > 1.00001) || (mean(variedScale) < 0.99999) % safeguard
            error('Mistake in keffScale rescaling');
        end

    %   Assign scales to model and record
        idx = 1;
        for j = 1:length(fullCplx)
            if ~variedKeff(j)
                init_keffScale(j,i) = 1;
            else
                init_keffScale(j,i) = variedScale(idx);

    %           First find all enzymeForm rxns
                cplxIdx = find(strcmp(model_alt.mets,['cplx_',fullCplx{j}]));
                rxnIdx = find(model_alt.S(cplxIdx,:)<0);
                for k = 1:length(rxnIdx)
                    model_alt.S(cplxIdx,rxnIdx(k)) = model_alt.S(cplxIdx,rxnIdx(k))/variedScale(idx);
                end

                idx = idx + 1;
            end
        end

    %   Optimize and collect results
        FBAsol_alt = overlayMultiomicsData(model_alt,data,100,waiverList);
        init_flux(:,i) = FBAsol_alt.full;
        init_obj(i) = FBAsol_alt.obj;
        
        model_alt = changeRxnBounds(model_alt,'EX_pi_e',-0.1,'l');
        FBAsol_alt = overlayMultiomicsData(model_alt,data_dep,100,waiverList);
        init_flux_dep(:,i) = FBAsol_alt.full;
        init_obj(i) = init_obj(i) + FBAsol_alt.obj;
    end
end

%% Randomize the original velocity

init_velo = zeros(length(fullCplx),noPart);

for i = 1:noPart
    for j = 1:length(fullCplx)
        
%       Won't waste time on constant keffs
        if ~variedKeff(j)
            continue;
        end
        
%       Two criteria: mean = 0 (no change of avg keff)
%                     mean(abs) ~ 0.05
        init_velo(j,i) = random('Uniform',-0.1,0.1);
    end
    
%   Correcting mean to 0
    mn = mean(init_velo(:,i));
    for j = 1:length(fullCplx)
        if variedKeff(j)
            init_velo(j,i) = init_velo(j,i) - mn*length(fullCplx)/length(find(variedKeff));
        end
    end

%   UNIQUE TO APPROACH 1
%   If contains var 'floored', fix their respective v to 0
%   their residue value is redistributed among non-floored varied keffs
    if exist('floored','var') && false
        sm = 0;
        for j = 1:length(floored)
            if floored(j)
                sm = sm + init_velo(j,i);
                init_velo(j,i) = 0;
            end
        end

        for j = 1:length(floored)
            if variedKeff(j) && ~floored(j)
                init_velo(j,i) = init_velo(j,i) + sm / (length(find(variedKeff))-length(find(floored)));
            end
        end
    end

end

clear cplxIdx i idx j k mn rxnIdx;

%% Start Iteration
tic

% Parameters
maxIt = 150;
inert = 0.8;
pBestWght = 0.1;
gBestWght = 0.2;
ksFloor = 0.1;
convThres = 0.001;

warning('off');

if ~exist('record_keffScale','var')

%   Configure records
    record_flux = init_flux;
    record_flux_dep = init_flux_dep;
    record_keffScale = init_keffScale;
    record_obj = init_obj;
    record_velo = init_velo;
    record_pBest = record_keffScale;

    if ~exist('floored','var') % if floored exists, inherent the previous gBest
        record_gBest = init_keffScale(:,find(init_obj == min(init_obj)));
    end
    
%   Starting
    for k = 2:maxIt
        fprintf('Running iteration %d out of %d\n',k,maxIt);
        
%       Calculate new velo vector for all cplx and record in record_velo
        velo_new = zeros(length(fullCplx),noPart);
        
        for i = 1:noPart
%           For each particle: V = inert * V(t-1) 
%                              + rand([0,1]) * pBestWght * (pBest-current) 
%                              + rand([0,1]) * gBestWght * (gBest-current)
            velo_new(:,i) = inert*record_velo(:,i,k-1) + rand()*pBestWght*(record_pBest(:,i)-record_keffScale(:,i,k-1)) + rand()*gBestWght*(record_gBest-record_keffScale(:,i,k-1));
        end
        
%       Move all particles accordingly
        keffScale_new = record_keffScale(:,:,k-1);
        
        for i = 1:noPart
            keffScale_new(:,i,1) = keffScale_new(:,i,1) + velo_new(:,i);
            
%           Raise all keffscale above the floor
%           the 'dificit' is distributed among all other keffs
            if any(keffScale_new(:,i,1) < ksFloor)
                df = 0;
                toRaise = find(keffScale_new(:,i,1) < ksFloor);
                
                for j = 1:length(toRaise) % calculate deficit
                    df = df + ksFloor - keffScale_new(toRaise(j),i,1);
                    keffScale_new(toRaise(j),i,1) = ksFloor;
                end
                
                for j = 1:length(fullCplx) % reducing all other keffscale
                    if variedKeff(j) && ~any(toRaise == j)
                        keffScale_new(j,i,1) = keffScale_new(j,i,1) - df/(length(variedScale)-length(toRaise));
                    end
                end
            end
        end
        
%       Calculate new objectives after moving
        obj_new = zeros(1,noPart);
        flux_new = zeros(length(model_adj.rxns),noPart);
        flux_new_dep = zeros(length(model_adj.rxns),noPart);
        
%       Pre-allocated cells used to feed and receive info regarding parfor
        models = cell(1,noPart);
        models_dep = cell(1,noPart);
        FBAsols = cell(1,noPart);
        FBAsols_dep = cell(1,noPart);

        for i = 1:noPart
            
%           Setting up model based on each vector in keffScale_new
            model_alt = model_adj;
            
            for j = 1:length(fullCplx)
                if keffScale_new(j,i,1) ~= 1
                    cplxIdx = find(strcmp(model_alt.mets,['cplx_',fullCplx{j}]));
                    rxnIdx = find(model_alt.S(cplxIdx,:)<0);
                    for m = 1:length(rxnIdx)
                        model_alt.S(cplxIdx,rxnIdx(m)) = model_alt.S(cplxIdx,rxnIdx(m))/keffScale_new(j,i,1);
                    end
                end
            end
            
%           Save models in the cell
            models{i} = model_alt;
            models_dep{i} = changeRxnBounds(model_alt,'EX_pi_e',-0.1,'l');
        end
        
%       Solve all using parfor and save in FBAsol map
%         parfor (i = 1:noPart,4)
        for i = 1:noPart
            FBAsols{i} = overlayMultiomicsData(models{i},data,100,waiverList);
            FBAsols_dep{i} = overlayMultiomicsData(models_dep{i},data_dep,100,waiverList);
        end
        
%       Extract results from parfor
        for i = 1:noPart
            FBAsol_alt = FBAsols{i};
            flux_new(:,i) = FBAsol_alt.full;
            obj_new(i) = FBAsol_alt.obj;
            
            FBAsol_alt = FBAsols_dep{i};
            flux_new_dep(:,i) = FBAsol_alt.full;
            obj_new(i) = obj_new(i) + FBAsol_alt.obj;
        end
        
%       Record all
        record_velo(:,:,k) = velo_new;
        record_flux(:,:,k) = flux_new;
        record_keffScale(:,:,k) = keffScale_new;
        record_obj(:,:,k) = obj_new;
        
%       Update gBest and pBest
        gBestObj = 100000;
        
        for i = 1:noPart
            [obj,idx] =  min(record_obj(1,i,:));
            record_pBest(:,i) = record_keffScale(:,i,idx);
            
            if obj < gBestObj
                gBestObj = obj;
                record_gBest = record_pBest(:,i);
            end
        end
        
%       Check if already converged
        if all(all(abs(velo_new(:,:)) < convThres))
            fprintf('Finished on convergence\n');
            lastIt = k;
            break;
        end
        
%       Check if maxIt is reached without converging
        if k == maxIt
            fprintf('Finished unconverged\n');
            lastIt = maxIt;
        end
        
    end
end

toc
clear i j k m cplxIdx rxnIdx df obj_new flux_new velo_new keffScale_new idx obj models FBAsols;
warning('on');

%% Analysis

% Calculate gBestObj as a function of iteration
record_gBestObj = zeros(1,lastIt);
record_gBestObj(1) = min(record_obj(1,:,1));

for i = 2:length(record_gBestObj)
    if min(record_obj(1,:,i)) < record_gBestObj(i-1)
        record_gBestObj(i) = min(record_obj(1,:,i));
    else
        record_gBestObj(i) = record_gBestObj(i-1);
    end
end

% Calculate grant avg speed over iterations
record_avgSpeed = zeros(1,lastIt);

for i = 1:length(record_avgSpeed)
    sm = 0;
    for j = 1:noPart
        sm = sm + sumsqr(record_velo(:,j,i))^0.5;
    end
    record_avgSpeed(i) = sm/noPart;
end

% Apply the result back to model
model_alt = model_adj;

% Pi replete
for i = 1:length(fullCplx)
    if (record_gBest(i) ~= 1) && (record_gBest(i) > 0)
        cplxIdx = find(strcmp(model_alt.mets,['cplx_',fullCplx{i}]));
        rxnIdx = find(model_alt.S(cplxIdx,:)<0);
        for j = 1:length(rxnIdx)
            model_alt.S(cplxIdx,rxnIdx(j)) = model_alt.S(cplxIdx,rxnIdx(j))/record_gBest(i);
        end
    end
end
[FBAsol_alt,model_alt] = overlayTranscriptome(model_alt,overlaidT_rep,100,waiverList,'sqr');

X = [];
Y = [];
lbl = {};
for i = 1:length(model_adj.rxns)
    if contains(model_adj.rxns{i},'EX_protein_') && ~any(strcmp(waiverList,erase(model_adj.rxns{i},'EX_protein_')))
        Y(length(Y)+1,1) = model_alt.c(i);
        X(length(X)+1,1) = -FBAsol_alt.full(i);
        lbl{length(lbl)+1,1} = model_adj.rxns{i};
    end
end

% Pi_dep
model_alt = model_adj;
model_alt = changeRxnBounds(model_alt,'EX_pi_e',-0.1,'l');

for i = 1:length(fullCplx)
    if (record_gBest(i) ~= 1) && (record_gBest(i) > 0)
        cplxIdx = find(strcmp(model_alt.mets,['cplx_',fullCplx{i}]));
        rxnIdx = find(model_alt.S(cplxIdx,:)<0);
        for j = 1:length(rxnIdx)
            model_alt.S(cplxIdx,rxnIdx(j)) = model_alt.S(cplxIdx,rxnIdx(j))/record_gBest(i);
        end
    end
end
[FBAsol_alt,model_alt] = overlayMultiomicsData(model_alt,data_dep,100,waiverList);

X_dep = [];
Y_dep = [];
lbl_dep = {};
for i = 1:length(model_adj.rxns)
    if contains(model_adj.rxns{i},'EX_protein_') && ~any(strcmp(waiverList,erase(model_adj.rxns{i},'EX_protein_')))
        Y_dep(length(Y_dep)+1,1) = model_alt.c(i);
        X_dep(length(X_dep)+1,1) = -FBAsol_alt.full(i);
        lbl_dep{length(lbl_dep)+1,1} = model_adj.rxns{i};
    end
end

% Plotting
figure;

subplot(2,2,1);
plot(1:maxIt,reshape(mean(record_obj(1,:,:)),[1,maxIt]),1:maxIt,record_gBestObj);
title('Avg objective values across particles');
xlabel('iteration');
ylabel('avg objective values');
set(gca,'YScale','log');

subplot(2,2,2);
plot(1:maxIt,record_avgSpeed);
title('Average particle speed over iterations');
xlabel('iteration');
ylabel('particle speed');

subplot(2,2,3);
scatter(X,Y);
title('Predicted proteome vs transcriptome Pi Rep');
xlabel('predicted');
ylabel('transcriptome');

subplot(2,2,4);
scatter(X_dep,Y_dep);
title('Predicted proteome vs transcriptome Pi Dep');
xlabel('predicted');
ylabel('transcriptome');

clear X X_dep Y Y_dep lbl lbl_dep cplxIdx rxnIdx mn i j k;

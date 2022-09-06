% Application of new merged model
changeCobraSolver('gurobi','all',0);
changeCobraSolverParams('LP','feasTol',1e-9);
changeCobraSolverParams('LP','optTol',1e-6);

if ~exist('model_ori','var')
    load('model_merged');
    model_ori = model_merged;
end

if ~exist('fst','var')
    fst = fastaread('/Users/16hy16/Documents/MATLAB/projects/Algal/model_recon/fasta/Chlamydomonas_reinhardtii_full.txt');
    fst = struct2table(fst);
end

model = model_ori;

%% TAG Production Enhancement

if ~exist('data_tag','var')
    data_tag = readtable('GSE55253_Expression_GEO_2014Feb.txt');
    data_tagNew = readtable('output_TAG_transcriptome_nuc.txt');
end

% Add tag production and dilution reactions
model_alt = model;
for i = 1314:1355
    model_alt = addReaction(model_alt,['TAG_prod_',model.mets{i}],...
        'metaboliteList',{model.mets{i},'TAG'},...
        'stoichCoeffList',[-1,1],...
        'reversible',false);
end

model_alt = addExchangeRxn(model_alt,'TAG',0,1000);

% Allow acetate uptake and photon exchange
% No nitrogen of any form
model_alt = changeRxnBounds(model_alt,'EX_ac_e',-1000,'l');
model_alt = changeRxnBounds(model_alt,'EX_TAG',1e-6,'l');
model_alt = changeRxnBounds(model_alt,'EX_nh4_e',0,'l');
% model_alt = changeRxnBounds(model_alt,'Biomass_Chlamy_auto',1e-3,'l');

% Compose the dataset and waiverList
fullProteinIdx = find(startsWith(model_ori.mets,'protein_'));
fullProtein = model_ori.mets(fullProteinIdx);
proteinExIdx = find(startsWith(model_ori.rxns,'EX_protein_'));

if ~exist('Trscpt','var')
    Trscpt = zeros(length(fullProtein),1);
    waiver = {};
    
    for i = 1:length(Trscpt)
        prot = erase(fullProtein{i},'protein_');
        idx = find(contains(data_tagNew.Gene,prot));
    
        if ~isempty(idx)
            Trscpt(i,1) = data_tagNew.SRR1174401(idx);
            Trscpt(i,2) = data_tagNew.SRR1174402(idx);
            Trscpt(i,3) = data_tagNew.SRR1174403(idx);
            Trscpt(i,4) = data_tagNew.SRR1174404(idx);
            Trscpt(i,5) = data_tagNew.SRR1174405(idx);
            Trscpt(i,6) = data_tagNew.SRR1174406(idx);
            Trscpt(i,7) = data_tagNew.SRR1174407(idx);
            Trscpt(i,8) = data_tagNew.SRR1174408(idx);
            Trscpt(i,9) = data_tagNew.SRR1174409(idx);
            
            Trscpt(i,10) = data_tagNew.SRR1174410(idx);
            Trscpt(i,11) = data_tagNew.SRR1174411(idx);
            Trscpt(i,12) = data_tagNew.SRR1174412(idx);
            Trscpt(i,13) = data_tagNew.SRR1174413(idx);
            Trscpt(i,14) = data_tagNew.SRR1174414(idx);
            Trscpt(i,15) = data_tagNew.SRR1174415(idx);
            Trscpt(i,16) = data_tagNew.SRR1174417(idx);
        else
            Trscpt(i,:) = zeros(1,16);
            waiver{end+1,1} = fullProtein{i};
        end
    end
end

% The function can only take on a maximum of 4 datasets but there're 16
%   hierarchical clustering for feature selection
figure;
dendrogram(linkage(pdist(Trscpt')));
title('Hierarchical Clustering of all 16 datasets');

% Define a new dataset of 4 col as distinctive features for keff estimation
%   col 1: avg of set 1 and 2
%   col 2: avg of set 3 to 6
%   col 3: avg of set 7 to 11
%   col 4: avg of set 12 to 16
feat = [];
feat(:,1) = mean(Trscpt(:,1:3),2);
feat(:,2) = Trscpt(:,4);
feat(:,3) = mean(Trscpt(:,5:8),2);
feat(:,4) = mean(Trscpt(:,9:16),2);

% relieving the constraints on proteins
% for i = 1:length(fullProteinIdx)
%     model_alt.csense(fullProteinIdx(i)) = 'G';
% end
for i = 1:length(fullProtein)
    model_alt = addReaction(model_alt,['DL_',fullProtein{i}],...
        'metaboliteList',fullProtein(i),'stoichCoeffList',[-1],...
        'reversible',false);
end

% model_alt = addExchangeRxn(model_alt,'proteinWC',0,150);
% model_alt.b(find(strcmp(model_alt.mets,'proteinWC'))) = 0;
% model_alt.b(find(strcmp(model_alt.mets,'proteinWC'))) = 0;

% estimating rate constants
if ~exist('rValues','var')
    [FBAsol,model_fit] = overlayMultiomicsData(model_alt,feat,0,waiver,...
        'keffEstimate',true);
    rIdx = find(startsWith(model_fit.varnames,'R_'));
    rValues = FBAsol.x(rIdx);
end

% plot results
varlen = length(model_alt.rxns) + length(find(startsWith(model_alt.rxns,'enzymeForm_')));
idx = [proteinExIdx, proteinExIdx+varlen, proteinExIdx+2*varlen, proteinExIdx+3*varlen];

figure;
for i = 1:4
    subplot(2,2,i);
    scatter(-FBAsol.x(idx(:,i)),feat(:,i));
%     scatter(feat(:,2),-FBAsol.x(idx(:,i)));
    title(['Feature ',num2str(i),' Fitting Result']);
end

clear idx varlen;

% Update Keff
fullCplxIdx = find(startsWith(model.mets,'cplx_'));
cplxFormIdx = find(startsWith(model.rxns,'cplxForm_'));
enzymeFormIdx = find(startsWith(model.rxns,'enzymeForm_'));

for i = 1:length(rValues)
    model_alt.S(fullCplxIdx(i),enzymeFormIdx) = model_alt.S(fullCplxIdx(i),enzymeFormIdx) / rValues(i);
end
%}
% Convex QP for proteome abundance
changeCobraSolverParams('QP','feasTol',1e-9);
changeCobraSolverParams('QP','printLevel',0);

record_FBAsol = [];
Trscpt_fit = [];
for i = 1:16
    if i == 1
        model_alt = changeRxnBounds(model_alt,'EX_nh4_e',-1000,'l');
    else
        model_alt = changeRxnBounds(model_alt,'EX_nh4_e',0,'l');
    end

%   Try scaling the data so lower
%   transcripts won't be overlooked
%   A better idea might be to re-define convex-qp objective
    dt = Trscpt(:,i);

    wt = dt.^1;
    wt(find(wt == 0)) = 1;
    wt = (1./wt);
    wt = wt / mean(wt);

%   Figure out the protein abundance vector by convex overlay
    [FBAsol_alt,model_qp] = overlayMultiomicsData(model_alt,dt,0,waiver,'objWeight',wt);
    disp(FBAsol_alt.obj);
%   [FBAsol_alt,~] = overlayMultiomicsData(model_alt,dt,0,waiver);
    record_FBAsol(:,i) = FBAsol_alt.full;
    Trscpt_fit(:,i) = model_qp.c(proteinExIdx)./wt;
end

% Plot fitting results
figure;
for i = 1:16
    subplot(4,4,i);
    ol = find((-record_FBAsol(proteinExIdx,i)./Trscpt_fit(:,i) > 3) + (-record_FBAsol(proteinExIdx,i)./Trscpt_fit(:,i) < 1/3));
    ol2 = find((-record_FBAsol(proteinExIdx,i)./Trscpt_fit(:,i) > 10) + (-record_FBAsol(proteinExIdx,i)./Trscpt_fit(:,i) < 0.1));
    scatter(Trscpt_fit(:,i),-record_FBAsol(proteinExIdx,i));
    hold on;
    scatter(Trscpt_fit(ol,i),-record_FBAsol(proteinExIdx(ol),i),[],[0.9290 0.6940 0.1250]);
    scatter(Trscpt_fit(ol2,i),-record_FBAsol(proteinExIdx(ol2),i),[],'red');
    hold off;

    set(gca,'YScale','log');
    set(gca,'XScale','log');
end

% Find the outlier...
proteinOL = zeros(length(fullProtein),1);
for i = 1:16
    idx = find((-record_FBAsol(proteinExIdx,i)./Trscpt_fit(:,i) > 10) + (-record_FBAsol(proteinExIdx,i)./Trscpt_fit(:,i) < 0.1));
    for j = 1:length(idx)
        if any(strcmp(waiver,fullProtein(idx(j))))
            continue;
        else
            proteinOL(idx(j)) = proteinOL(idx(j)) + 1;
        end
    end
end

% Checking TAG production limits
tag_ub = [];
tag_lb = [];
models = {}; % save all 16 models

for i = 1:16
    if i == 1
        model_alt = changeRxnBounds(model_alt,'EX_nh4_e',-1000,'l');
    else
        model_alt = changeRxnBounds(model_alt,'EX_nh4_e',0,'l');
    end

%   Fix the protein abundance vector and check TAG variability
    model_alt.lb(proteinExIdx) = record_FBAsol(proteinExIdx,i);
    model_alt.ub(proteinExIdx) = record_FBAsol(proteinExIdx,i)*0.98;
    model_alt = changeObjective(model_alt,'EX_TAG');

%   fix any falsely positive protein bounds
    idx = find(record_FBAsol(proteinExIdx,i) > 0);
    model_alt.lb(proteinExIdx(idx)) = 0;
    model_alt.ub(proteinExIdx(idx)) = 0;

%   Relieve protein abundance bound for waivered proteins
    for j = 1:length(waiver)
        idx = find(strcmp(model_alt.rxns,['EX_',waiver{j}]));
        model_alt.lb(idx) = -1000;
        model_alt.ub(idx) = 0;
    end

%   Conduct a thing like FVA
    FBAsol_alt = optimizeCbModel(model_alt,'min');
    tag_lb(i,1) = FBAsol_alt.f;
    FBAsol_alt = optimizeCbModel(model_alt,'max');
    tag_ub(i,1) = FBAsol_alt.f;
    models{i} = model_alt;
end

% Plot TAG production
% v_tag = record_FBAsol(18287,:);
% x = {'0h','0.5h','2h','4h','8h','12h','24h','48h','b','b+0.5h','b+2h','b+4h','b+8h','b+12h','b+24h','b+48h'};
% X = categorical(x);
% X = reordercats(X,x);
% figure;
% bar(X,[tag_lb,tag_ub]);
% set(gca,'YScale','log');
% ylabel('Total TAG Prod. Flux');
% xlabel('Time');

clear idx prot X x dt wt i j model_qp ol ol2;

% Two Parameter: transcriptome exponent and protein UB multiplier
% Now the production rate is incredibly low (~2e-6)

%% De-bottlenecking
%   max production rate is 0.59 without proteomic bounds, so relieving
%   certain bounds gotta be dramastically increasing the production

models_db = cell(16,1);
FBAsols = cell(16,1);
relaxList = zeros(length(fullProtein),16);
% test = [];

for i = 1:16
    model_alt = models{i};
    [FBAsols{i},models_db{i},relaxProt,relaxLevel] = proteinDebottleneck(model_alt,20);

    for j = 1:length(relaxProt)
%         relaxList(find(strcmp(fullProtein,relaxProt{j}))) = relaxList(find(strcmp(fullProtein,relaxProt{j}))) + relaxLevel(j);
        relaxList(find(strcmp(fullProtein,relaxProt{j})),i) = relaxLevel(j);
    end
    disp(FBAsols{i}.f);
end

% plot the 15 most relieved proteins across samples
[~,idx] = maxk(sum(relaxList,2),25);
figure;
h = heatmap(1:16,erase(fullProtein(idx),'protein_'),relaxList(idx,:));
xlabel('Sample');
ylabel('Debottlenecked Protein ID');

%% Unbiased Analysis
% Now do unbiased analysis on debottlenecked models...

% Attempt 1:
% Relieving TAG bound
% model_alt.lb(15069) = FBAsols{15}.f*0.90;
% 
% if ~exist('sampleStruct','var')
%     bs.method = 'uniform';
%     bs.index = 15069;
%     bs.param = [FBAsols{15}.f*0.95,FBAsols{15}.f*0.999];
%     [sampleStruct,mf] = gpSampler(model_alt,10000,[],600);
%     sampleStruct = chrrSampler(model_alt,1,30000,0,{},false,0);
% end
% 
% save('struct_dataset15_ge95_2','sampleStruct','-v7.3');
% 
% figure;
% histogram(sampleStruct.points(15069,:));
% set(gca,'YScale','log');
%}

% Attempt 2:
% Using ACHR sampler...
%   gotta add dilution reactions for proteins and proteinWC
% if ~exist('warmupPts','var')
%     warmupPts = createHRWarmup(model_alt,10000);
%     save('ACHRWarmUpPts','warmupPts','-v7.3');
% end
% 
% for i = 2:5
%     load('ACHR_last_point.mat');
%     ACHRSampler(model_alt,warmupPts,'ACHR_sample15_ge90',1,10000,10,curPoint,i-1);
% end

% Attempt 3:
% FVA to figure out v lb and ub -> set lb and ub to M-model -> sample
% 
% first do FVA only on metabolic reactions
% load('model_alt.mat');
% model_alt.lb(15069) = 0;
% fullReaction = [model_alt.rxns(1:2641);model_alt.rxns(15027:15069)];
% [minV,maxV] = fluxVariability(model_alt,0,'max',fullReaction);
% 
% % remove all PC structures
% rmMets = model_alt.mets([find(startsWith(model_alt.mets,'protein_'));...
%     find(startsWith(model_alt.mets,'cplx_'));...
%     find(startsWith(model_alt.mets,'enzyme_'));...
%     find(strcmp(model_alt.mets,'proteinWC[c]'))]);
% rmRxns = model_alt.rxns([find(startsWith(model_alt.rxns,'EX_protein_'));...
%     find(startsWith(model_alt.rxns,'cplxForm_'));...
%     find(startsWith(model_alt.rxns,'enzymeForm_'));...
%     find(startsWith(model_alt.rxns,'EX_enzyme_'));...
%     find(startsWith(model_alt.rxns,'DL_protein'));...
%     find(startsWith(model_alt.rxns,'EX_proteinWC'))]);
% 
% model_m = removeMetabolites(model_alt,rmMets,false);
% model_m = removeRxns(model_m,rmRxns,'metFlag',false);
% 
% % apply lb and ub and sample
% model_m.lb = minV;
% model_m.ub = maxV;
% 
% [sampleStruct,mf] = gpSampler(model_m,100000);
% figure;
% histogram(sampleStruct.points(2684,:));
% set(gca,'YScale','log');

% Using only FVA
fullReaction = [model_alt.rxns(1:2641);model_alt.rxns(17671:17713)];

% save FBA results: 0%, 50%, 90%, 99%
if ~exist('FVAsols','var')
    FVAsols = zeros(length(fullReaction),16);
    for i = 1:16
        disp(i);
        model_alt = models_db{i};
        model_alt.lb(17713) = 0;

        FBAsol_alt = optimizeCbModel(model_alt,'max');
        optPerc = [0,0.5,0.9,0.99];
        for j = 1:4
            model_alt.lb(17713) = FBAsol_alt.f * optPerc(j);
            for k = 1:length(fullReaction)
                model_alt = changeObjective(model_alt,fullReaction{k});
                FVAsol = optimizeCbModel(model_alt,'min');
                FVAsols(k,i,j) = FVAsol.f;
            end
        end
        for j = 1:4
            model_alt.lb(17713) = FBAsol_alt.f*optPerc(5-j);
            for k = 1:length(fullReaction)
                model_alt = changeObjective(model_alt,fullReaction{k});
                FVAsol = optimizeCbModel(model_alt,'max');
                FVAsols(k,i,j+4) = FVAsol.f;
            end
        end
%         [FVAsols(:,i,1),FVAsols(:,i,8)] = fluxVariability(model_alt,0,'max',fullReaction);
%         [FVAsols(:,i,2),FVAsols(:,i,7)] = fluxVariability(model_alt,50,'max',fullReaction);
%         [FVAsols(:,i,3),FVAsols(:,i,6)] = fluxVariability(model_alt,90,'max',fullReaction);
%         [FVAsols(:,i,4),FVAsols(:,i,5)] = fluxVariability(model_alt,99,'max',fullReaction);
    end
end

% find rxns with few FVA overlaps in sample 3 and 15
variedRxns = {};
overlapRatio = 0.5;
sample1 = 3;
sample2 = 9;

for i = 1:length(fullReaction)
    if (max(FVAsols(i,sample1,:)) < 0.2) && (max(FVAsols(i,sample2,:)) < 0.2)
        continue;
    end

    if FVAsols(i,sample1,1) > FVAsols(i,sample2,8)
        variedRxns(end+1,1) = fullReaction(i);
    elseif FVAsols(i,sample2,1) > FVAsols(i,sample1,8)
        variedRxns(end+1,1) = fullReaction(i);
    elseif FVAsols(i,sample1,8) > FVAsols(i,sample2,8)
        if FVAsols(i,sample1,1) > FVAsols(i,sample2,1)
            if overlapRatio > (FVAsols(i,sample2,8)-FVAsols(i,sample1,1)) / (FVAsols(i,sample1,8)-FVAsols(i,sample2,1))
                variedRxns(end+1,1) = fullReaction(i);
            end
        else
            if overlapRatio > (FVAsols(i,sample2,8)-FVAsols(i,sample2,1)) / (FVAsols(i,sample1,8)-FVAsols(i,sample1,1))
                variedRxns(end+1,1) = fullReaction(i);
            end
        end

    elseif FVAsols(i,sample1,8) < FVAsols(i,sample2,8)
        if FVAsols(i,sample1,1) > FVAsols(i,sample2,1)
            if overlapRatio > (FVAsols(i,sample1,8)-FVAsols(i,sample1,1)) / (FVAsols(i,sample2,8)-FVAsols(i,sample2,1))
                variedRxns(end+1,1) = fullReaction(i);
            end
        else
            if overlapRatio > (FVAsols(i,sample1,8)-FVAsols(i,sample2,1)) / (FVAsols(i,sample2,8)-FVAsols(i,sample1,1))
                variedRxns(end+1,1) = fullReaction(i);
            end
        end
    end
end

% plot FVA results
% figure;
% subplot(2,2,1);
% scatter(abs(FVAsols(:,2,1)),abs(FVAsols(:,15,1)));
% hold on;
% scatter(abs(FVAsols(:,2,8)),abs(FVAsols(:,15,8)));
% hold off;
% legend({'V_{min}','V_{max}'});
% set(gca,'YScale','log');
% set(gca,'XScale','log');
% title('0%');
% 
% subplot(2,2,2);
% scatter(abs(FVAsols(:,2,2)),abs(FVAsols(:,15,2)));
% hold on;
% scatter(abs(FVAsols(:,2,7)),abs(FVAsols(:,15,7)));
% hold off;
% legend({'V_{min}','V_{max}'});
% set(gca,'YScale','log');
% set(gca,'XScale','log');
% title('50%');
% 
% subplot(2,2,3);
% scatter(abs(FVAsols(:,2,3)),abs(FVAsols(:,15,3)));
% hold on;
% scatter(abs(FVAsols(:,2,6)),abs(FVAsols(:,15,6)));
% hold off;
% legend({'V_{min}','V_{max}'});
% set(gca,'YScale','log');
% set(gca,'XScale','log');
% title('90%');
% 
% subplot(2,2,4);
% scatter(abs(FVAsols(:,2,4)),abs(FVAsols(:,15,4)));
% hold on;
% scatter(abs(FVAsols(:,2,5)),abs(FVAsols(:,15,5)));
% hold off;
% legend({'V_{min}','V_{max}'});
% set(gca,'YScale','log');
% set(gca,'XScale','log');
% title('100%');

% plot some reactions...
% 
% sample 1 - 3
% rxnList = {'ENO','PYK','ICL','MALSm','ATPSm','NADHOR','LLDH_ferr_m',...
%     'PFLACTm','CS','CSm','FUMm','SUCDH_q8_m','SUCDHm','PPCKm','SO4t',...
%     'NanoG0198_chl','Tr_ATP(3h)t_h_chl','Tr_PEP(pi)th_c_chl','P680P_chl',...
%     'P700A0_chl','CEF_chl'};
% sample 3 - 9
% rxnList = {'PPC','ENO','ALCDH_pqq_x','ICL','MALSm','ATPSm','CYOR_q8_m',...
%     'ATGDm','ATUDm','CMP','UMPP','UTCY','UTUP','ACALD','LLDH_ferr_m'...
%     'PFLACTm','SUCL_gdp_m','CSm','PCr','FUMm','SUCDH_q8_m','SUCDHm',...
%     'PPCKm','SO4t','NanoG0198_chl','P680P_chl','P700A0_chl','CEF_chl'};
% sample 9 - 15
% rxnList = {'CSm','FUMm','SUCDH_q8_m','SUCDHm','P680P_chl','P700A0_chl'};
% 
% FA and TAG Synthesis
% rxnList = {'COA1819ZD9DS','AGPATCOA1601819Z','PAPA1601819Z','ACOADAGAT1601819Z160',...
%     'PLDAGAT1601819Z1602','EX_TAG'};
% rxnList = {'AGPATCOA1601819Z','PAPA1601819Z','ACOADAGAT1601819Z160',...
%     'PLDAGAT1601819Z1602','EX_TAG'};
% protIdx = {359,[],[352,353,354,355,356,357],395,[]};
% rev = [0,0,0,0,0];
% ACCOA Supply
% rxnList = {'ACALD','MMSD','PFLACTm','ACS','PTArm'};
% protIdx = {149,551,926,523,[927,928]};
% rev = [0,0,0,0,1];
rxnList = {'ACALD','MMSD','PFLACTm','ACS','PTArm','ACKrm'};
protIdx = {149,551,926,523,[927,928],[922,923]};
rev = [0,0,0,0,1,0];
% TCA Cycle
% rxnList = {'CS','CSm','ACONT','ACONTm','ICL','ICDH','ICDH_nad',...
%     'ICDHm','SUCLm','SUCDH_q8_m','SUCDHm','FUMm','MDHCm','MALSm'};
% ATP synthesis
% rxnList = {'PYK','PYKm','ATPSm','NanoG0425_chl','ATPSh_chl'};
% protIdx = {[499,500,501,502,503,504],[499,500,501,502,503,504],...
%     find(C_matrix(:,511)),[504,1419],find(C_matrix(:,1184))};
% rev = [0,0,0,0,0];
% rxnList = {'PYK','PYKm','ATPSm','NanoG0425_chl','ATPSh_chl','PGK'};
% protIdx = {[499,500,501,502,503,504],[499,500,501,502,503,504],...
%     find(C_matrix(:,511)),[504,1419],find(C_matrix(:,1184)),497};
% rev = [0,0,0,0,0,1];
% De novo synthesis
% rxnList = {'HACD2m','ECOAH2m','ACOAR2m','ACACT3m'};
% protIdx = {[142,264],144,127,265};
% rev = [0,0,0,0];
% Chloroplast Activities
% rxnList = {'NanoG0589_chl','NanoG0198_chl','P680P_chl','P700A0_chl',...
%     'CEF_chl','ATPSh_chl','Tr_ATP(3h)t_h_chl','Tr_G3P(pi)thr_c_chl'};
% Miscellaneous
% rxnList = {'ENO','ALDDHm','PPCKm','NA1Hth','SO4t'};

figure;
FVAsols_plt = FVAsols;

for i = 1:length(rxnList)
    subplot(5,5,i);
    yyaxis left;

    idx = find(strcmp(fullReaction,rxnList{i}));
    if rev(i) == 1
        for j = 1:16
            FVAsols_plt(idx,j,:) = -flip(FVAsols_plt(idx,j,:));
        end
    end

    h1 = bar(1:16,reshape([FVAsols_plt(idx,:,1),FVAsols_plt(idx,:,8)-FVAsols_plt(idx,:,1)],[16,2]),...
        0.8,'stacked','BaseValue',min(FVAsols_plt(idx,:,1)));
    h1(1).Visible = 'off';
    h1(2).FaceColor = [0.4660 0.6740 0.1880];
    hold on;
    h2 = bar(1:16,reshape([FVAsols_plt(idx,:,2),FVAsols_plt(idx,:,7)-FVAsols_plt(idx,:,2)],[16,2]),0.65,'stacked');
    h2(1).Visible = 'off';
    h2(2).FaceColor = [0.9290 0.6940 0.1250];
    h3 = bar(1:16,reshape([FVAsols_plt(idx,:,3),FVAsols_plt(idx,:,6)-FVAsols_plt(idx,:,3)],[16,2]),0.5,'stacked');
    h3(1).Visible = 'off';
    h3(2).FaceColor = [0.8500 0.3250 0.0980];
    h4 = bar(1:16,reshape([FVAsols_plt(idx,:,4),FVAsols_plt(idx,:,5)-FVAsols_plt(idx,:,4)],[16,2]),0.35,'stacked');
    h4(1).Visible = 'off';
    h4(2).FaceColor = [1 0 0];
    hold off;
%     ylabel('flux');

    if ~isempty(protIdx{i})
        yyaxis right;
        if length(protIdx{i}) == 1
            p = plot(Trscpt_fit(protIdx{i},:),'-k');
        else
            p = plot(sum(Trscpt_fit(protIdx{i},:)),'-k');
        end
        p.LineWidth = 2;
%         ylabel('expression');
    end
%     xlabel('sample no');

    title(rxnList{i},'FontSize',16,'Interpreter','none');
%     set(gca,'YScale','log');
    xline([1.5,9.5],'--',{'N-free','Ac Boost'});
%     set(gca,'XLim',[min(min(FVAPlot)),max(max(FVAPlot))]);
end
clear sample1 sample2 idx i j k h h1 h2 h3 h4 p protIdx;

%% Classify Enzymatic Reactions by FVA
% automatic classification using Spearman's corr

rhos = zeros(2684,1);
ps = zeros(2684,1);

for i = 1:2641
    x = zeros(16,1);
    y = zeros(16,1);

    if ~any(K_matrix(i,:))
        continue;
    end

    idx = find(K_matrix(i,:));
    for j = 1:length(idx)
%         x = x + (Trscpt_fit + relaxList)' * C_matrix(:,idx(j));
        x = x + (Trscpt_fit)' * C_matrix(:,idx(j));
    end

    for j = 1:16
        y(j) = max(abs(FVAsols(i,j,:)));
    end

    [rhos(i),ps(i)] = corr(x,y,'Type','Spearman');
end

clear i j x y;

% find rho > 0.7
figure;
histogram(rhos);

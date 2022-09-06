function model_new = mapRxnsToBioCyc(model,biocyc,startFrom)
% Semi-auto curator of M-model rules field using BioCyc DB
% Only works when the biocyc file is processed as anticipated

if ~exist('startFrom','var')
    startFrom = 1;
end

% check if fields already exists
if ~isfield(model,'rxnBioCycID')
    model.rxnBioCycID = cell(length(model.rxns),1);
    model.rxnBioCycIDCand = cell(length(model.rxns),1);
    model.rxnBioCycIDEnzRxn = cell(length(model.rxns),1);
end

term = false; % termination flag

for i = startFrom:length(model.rxns)

    if isempty(model.rules{i}) % continue if rule is unavailable
        continue;
    end

%   parse out genes
    gns = regexp(model.rules{i},'\d*','Match');
    gns = unique(gns);

%   scoring biocyc.enzrxns based on geneRules
    sc_enz = zeros(length(biocyc.enzrxns),1);

    for j = 1:length(gns)

%       find the respective biocyc gene
        idx = find(strcmp(biocyc.genes(:,3),model.genes{str2double(gns(j))}));
        if isempty(idx)
            idx = find(startsWith(biocyc.genes(:,3),[model.genes{str2double(gns(j))},',']));
        end
        if isempty(idx)
            idx = find(endsWith(biocyc.genes(:,3),[',',model.genes{str2double(gns(j))}]));
        end
        if isempty(idx)
            idx = find(contains(biocyc.genes(:,3),[',',model.genes{str2double(gns(j))},',']));
        end
        if isempty(idx)
            continue;
        end

%       scoring based on biocyc.C
        for k = 1:length(idx)
            er = find(biocyc.C(idx(k),:));
            if ~isempty(er)
                sc_enz(er) = sc_enz(er) + 1;
            end
        end
    end

%   if no enzrxn scores, no need to proceed
    if ~any(sc_enz)
        continue;
    end

%   scoring biocyc.rxns from biocyc.enzrxns and biocyc.K
    sc_rxn = biocyc.K * sc_enz;

%   rank biocyc.rxns from highest to lowest score
    [~,idx] = maxk(sc_rxn,length(find(sc_rxn)));

%   save all candidates to model.rxnBioCycIDCand
    model.rxnBioCycIDCand{i} = '';
    for j = 1:length(idx)
        model.rxnBioCycIDCand{i} = [model.rxnBioCycIDCand{i},biocyc.rxns{idx(j),1}];
        if j ~= length(idx)
            model.rxnBioCycIDCand{i} = [model.rxnBioCycIDCand{i},','];
        end
    end

%   let user decide which to map to model.rxns
    j = 1;
    while true
        fprintf('\n===================================================\n');
        fprintf('\nMapping reaction %d/%d: %s (%s)\n',i,length(model.rxns),model.rxns{i},model.rxnNames{i});
        fprintf('  <strong>%s</strong>\n',model.reactions{i});
        if isempty(model.rxnBioCycID{i})
            fprintf('  Unmapped\n');
        else
            fprintf('  Currently mapped to %s\n',model.rxnBioCycID{i});
        end
        fprintf('\n=============== TO =============== \n\n');
        fprintf('Candidate %d/%d\n',j,length(idx));
        fprintf('Reaction %d: %s\n',idx(j),biocyc.rxns{idx(j),1});
        fprintf('  <strong>%s</strong>\n',biocyc.rxns{idx(j),7});
        fprintf('  Parent of: \n')

%       print all enzrxn of this rxn
        er = find(biocyc.K(idx(j),:));
        if length(er) > 1
            for k = 1:length(er)
                if contains(biocyc.enzrxns{er(k),3},'CPLX')
                    fprintf('  - %d. %s (%s) (<strong>%s</strong>)\n',k,biocyc.enzrxns{er(k),1},biocyc.enzrxns{er(k),2},biocyc.enzrxns{er(k),3});
                else
                    fprintf('  - %d. %s (%s) (%s)\n',k,biocyc.enzrxns{er(k),1},biocyc.enzrxns{er(k),2},biocyc.enzrxns{er(k),3});
                end
                fprintf('       %s\n',biocyc.enzrxns{er(k),5});
            end
        else
            if contains(biocyc.enzrxns{er,3},'CPLX')
                fprintf('  - %s (%s) (<strong>%s</strong>)\n',biocyc.enzrxns{er,1},biocyc.enzrxns{er,2},biocyc.enzrxns{er,3});
            else
                fprintf('  - %s (%s) (%s)\n',biocyc.enzrxns{er,1},biocyc.enzrxns{er,2},biocyc.enzrxns{er,3});
            end
            fprintf('       %s\n',biocyc.enzrxns{er,5});
        end
        fprintf('\n===================================================\n');
        fprintf('\n');
        
        prompt = '[Y]es, [N]o, [R]estart, [P]ass, or [T]erminate?\n';

        ch = upper(input(prompt,'s'));
        
        switch ch
            case 'Y' % confirm mapping to biocyc.rxn
                model.rxnBioCycID{i} = biocyc.rxns{idx(j),1};

%               choosing EnzRxn only if more than 1 enzrxn is available
                if length(er) > 1
                    chosen = false;

                    while ~chosen
                        fprintf('\n===================================================\n');
                        fprintf('\nMapping enzrxns...\n');
                        prompt = '[A]ll, [N]one, or type in enzrxn indexes separated by comma (E.g. 1,2,4)\n';
                        ch = upper(input(prompt,'s'));

                        if strcmp(ch,'A') % add all enzrxn to model.rxnBioCycIDEnzRxn
                            model.rxnBioCycIDEnzRxn{i} = '';
                            for k = 1:length(er)
                                model.rxnBioCycIDEnzRxn{i} = [model.rxnBioCycIDEnzRxn{i},biocyc.enzrxns{er(k),1}];
                                if k ~= length(er)
                                    model.rxnBioCycIDEnzRxn{i} = [model.rxnBioCycIDEnzRxn{i},','];
                                end
                            end
                            chosen = true;

                        elseif strcmp(ch,'N') % leave a '?'
                            model.rxnBioCycIDEnzRxn{i} = '?';
                            chosen = true;

                        else % parse out and assign input enzrxn
                            ch = split(ch,',');
                            chosen = true;

                            try
                                model.rxnBioCycIDEnzRxn{i} = '';
                                for k = 1:length(ch)
                                    model.rxnBioCycIDEnzRxn{i} = [model.rxnBioCycIDEnzRxn{i},biocyc.enzrxns{er(str2double(ch(k))),1}];
                                    if k ~= length(ch)
                                        model.rxnBioCycIDEnzRxn{i} = [model.rxnBioCycIDEnzRxn{i},','];
                                    end
                                end
                            catch
                                warning('Error in input text');
                                chosen = false;
                            end
                        end
                    end

%               if theres only 1 enzrxn, choose automatically
                elseif length(er) == 1
                    model.rxnBioCycIDEnzRxn{i} = biocyc.enzrxns{er,1};
                else
                    error('Unknown error');
                end

                clc
                break;

            case 'N' % proceed to the next candidate inline
                if j == length(idx)
                    model.rxnBioCycID{i} = '?';
                    clc;
                    break;
                else
                    clc;
                    j = j + 1;
                end

            case 'R' % restart from the first candidate
                j = 1;
                clc;

            case 'P' % continue without making modification
                clc;
                break;

            case '' % directly hitting ENTER is also accepted as 'pass'
                clc;
                break;

            case 'T' % termination
                clc;
                term = true;
                break;

            otherwise
                clc;
                fprintf('Invalid Input\n');
        end
    end

%   termination flag check
    if term
        model_new = model;
        fprintf('Terminated\n');
        break;
    end
end

% Check if finished
if i == length(model.rxns)
    model_new = model;
    fprintf('Finished\n');
end

end
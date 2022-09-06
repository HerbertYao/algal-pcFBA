function enzymeString = parseGeneRule(ruleString)
% This function parse model.rules{i} into a array of all possible enzymes

ruleString = erase(ruleString,' '); % Remove all spaces
enzymeString = ''; % Initialize the return string

% If empty, return
if isempty(ruleString)
    enzymeString = ruleString;
    
% First check if the ruleString is only a single gene (no '&' or '|')
% if so, return it directly
elseif (~contains(ruleString,'&') && ~contains(ruleString,'|'))

%   remove the bracket if it exists
    if strcmp(ruleString(1),'(')
        ruleString(1) = [];
        ruleString(end) = [];
    end

    enzymeString = ruleString;
    
else
%   Else Looking for operators

%   determine if the most outer bracket exists or not by counting brackets
    while true
        numBracket = 0;
        minNumBracket = 10;
        for i = 1:length(ruleString)-1
            if strcmp(ruleString(i),'(')
                numBracket = numBracket + 1;
            elseif strcmp(ruleString(i),')')
                numBracket = numBracket - 1;
            end
    
            if numBracket < minNumBracket
                minNumBracket = numBracket;
            end
        end
    
%       remove the outer bracket if need to. If minNoBracket == 0, break
        if minNumBracket == 0
            break;

        elseif minNumBracket == 1
            ruleString(1) = [];
            ruleString(end) = [];

        elseif minNumBracket ~= 0 % don't allow negative cases
            error('Error in rules: %s',ruleString);
        end
    end

%   Find the first operator which isn't in a bracket
    numBracket = 0;
    operator = '';
    
    for i = 1:length(ruleString)
        
%       Finding exposed '|' or '&'
%       once the operator is determined, all following exposed operators
%       must be the same. Otherwise return error. E.g. (x(1)|x(2)&x(3))
        if (strcmp(ruleString(i),'|') && (numBracket == 0))
            if isempty(operator)
                operator = 'OR';
            elseif ~strcmp(operator,'OR')
                error('Contradiction in rules: %s',ruleString);
            end
            ruleString(i) = '$';
        
        elseif (strcmp(ruleString(i),'&') && (numBracket == 0))
            if isempty(operator)
                operator = 'AND';
            elseif ~strcmp(operator,'AND')
                error('Contradiction in rules: %s',ruleString);
            end
            ruleString(i) = '$';

%       Bracket counting
        elseif strcmp(ruleString(i),'(')
            numBracket = numBracket + 1;
        elseif strcmp(ruleString(i),')')
            numBracket = numBracket - 1;
        end
    end
    
%   First split the string by '$' after counting
    ruleStringFrac = split(ruleString,'$');
    
%   If it's an 'OR', split and each pass down to the next level.
%   Splitted strings are connected by ';'
    if strcmp(operator,'OR')
        for j = 1:length(ruleStringFrac)
            enzymeString = append(enzymeString,parseGeneRule(ruleStringFrac{j}));
            
            if j ~= length(ruleStringFrac)
                enzymeString = append(enzymeString,';');
            end
        end
        
%   If it's an 'AND', split and pass down to the next level, but don't
%       append directly
%   Decode the return of each call and save to the stringRecord, where
%       each return occupies a column
%   If a return contains ';', split it into several rows so each entry 
%       either is a single gene or in 'and' relation
%   Use combntns to generate a full combination and append each 
%       combination to enzymeString.
    elseif strcmp(operator,'AND')
        stringRecord = {};
        noStringRecord = [];
            
        for j = 1:length(ruleStringFrac)
            returnStrings = split(parseGeneRule(ruleStringFrac{j}),';');
            returnStrings(cellfun('isempty',returnStrings)) = [];
            
%           Record all return values and number of return strings
            noStringRecord(j) = length(returnStrings);
            for k = 1:length(returnStrings)
                stringRecord{k,j} = returnStrings(k);
            end
        end
                
%       Generate all combinations using function combvec
        comb = combvec(1:noStringRecord(1), 1:noStringRecord(2));
                
        if length(noStringRecord) >= 3
            for j = 3:length(noStringRecord)
                comb = combvec(comb, 1:noStringRecord(j));
            end
        end
                
        comb = comb.';
        [combLen, combWid] = size(comb);
        
%       Append to enzyme string
        for j = 1:combLen
            for k = 1:combWid
                enzymeString = append(enzymeString,stringRecord{comb(j,k),k});
            end
            
            if j ~= combLen
                enzymeString = append(enzymeString,';');
            end
        end
        
%   Return error if non of the two applies
    else
        error('Error in geneRule parsing.\n');
    end
end

end
function K_matrix_new = estimateKeffFromMW(C_matrix,K_matrix,proteinMW)

% Estimating every rate constant of enzymes using its molar weight while
% keeping the same avg k_eff
% Equation: Keff_cplx_i = Keff_avg * (MW_cplx_i / MW_avg)^0.75
% 
% USAGE:
%   K_matrix_adj = proteinConstraintModel(C_matrix,K_matrix,proteinMW);
% 
% INPUTS:
%   C_matrix:  
%   K_matrix:  
%   proteinMW: 
% 
% OUTPUTS:
%   K_matrix_new: Estimated K_matrix

% Record the original Keff
[~,enzymeLen] = size(C_matrix);
KeffList_old = zeros(enzymeLen,1);

for i = 1:length(KeffList_old)
    KeffList_old(i) = mean(K_matrix(find(K_matrix(:,i)),i)); % Only average non-zero keff
    
    if KeffList_old(i) ~= 234000 % safeguard
        warning('Wrong here\n');
    end
end

Keff_avg = mean(KeffList_old);

% First calculate enzyme molar weight
enzymeMW = zeros(enzymeLen,1);

for i = 1:length(enzymeMW)
    for j = 1:length(proteinMW)
        enzymeMW(i) = enzymeMW(i) + C_matrix(j,i)*proteinMW(j);
    end
end

% Calculate mean molar weight
enzymeMW_avg = mean(enzymeMW);

% Adjust rate constants using equation
%   Keff_cplx_i = Keff_avg * (MW_cplx_i / MW_avg)^0.75

KeffList = zeros(enzymeLen,1);

for i = 1:length(KeffList)
    KeffList(i) = Keff_avg * (enzymeMW(i)/enzymeMW_avg)^0.75;
end

% Rescale
KeffList = KeffList * (Keff_avg/mean(KeffList));

% Updated K matrix
K_matrix_new = K_matrix;

for i = 1:length(KeffList)
    idx = find(K_matrix_new(:,i));
    
    for j = 1:length(idx)
        K_matrix_new(idx(j),i) = KeffList(i);
    end
end

end
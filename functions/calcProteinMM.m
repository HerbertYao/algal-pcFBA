function MM = calcProteinMM(seq)
% The function returns the molar mass of a given protein sequence, which
% must be in one-letter symbols. The return unit is g/mol.

MM = 0;
seqList = num2cell(seq);
for i = 1:length(seqList)
    switch seqList{i}
        case '*'
            continue;
        case 'A'
            MM = MM + 71;
        case 'R'
            MM = MM + 156.1;
        case 'N'
            MM = MM + 114;
        case 'D'
            MM = MM + 115;
        case 'C'
            MM = MM + 103;
        case 'E'
            MM = MM + 129;
        case 'Q'
            MM = MM + 128.1;
        case 'G'
            MM = MM + 57;
        case 'H'
            MM = MM + 137.1;
        case 'I'
            MM = MM + 113.1;
        case 'L'
            MM = MM + 113.1;
        case 'K'
            MM = MM + 128.1;
        case 'M'
            MM = MM + 131;
        case 'F'
            MM = MM + 147.1;
        case 'P'
            MM = MM + 97.1;
        case 'S'
            MM = MM + 87;
        case 'T'
            MM = MM + 101;
        case 'W'
            MM = MM + 186.1;
        case 'Y'
            MM = MM + 163.1;
        case 'V'
            MM = MM + 99.1;
%       Use avg AA weight if unrecognized
        otherwise
            MM = MM + 118.8;
    end
end

% Add the molar mass of a water
MM = MM + 18.02;

end
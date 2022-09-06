function seq = findProteinSeq(geneID, fasta, avgLen)
% The function returns a protein sequence from fasta file. If the ID is not
% found, return an average-lengthed 'ZZZZZZZZ' seq

found = false;

% First try complete match
for i = 1:length(fasta)
    if strcmp(fasta(i).Header,geneID)
        seq = fasta(i).Sequence;
        found = true;
        break;
    end
end

% Then try contains
if ~found
    for i = 1:length(fasta)
        if contains(fasta(i).Header,geneID)
            seq = fasta(i).Sequence;
            found = true;
            break;
        end
    end
end

% If still not found, calculate the average length and return
if ~found
    seq(1:round(avgLen)) = 'Z';
    fprintf('Not found gene: %s\n',geneID);
end

end
    
% pretreat nucleotide fasta

% dt_ori = fastaread('nucleotide_all.txt');
dt_ori = fastaread('nucleotide_nuc.txt');

dt = dt_ori;
rm = [];

for i = 1:length(dt)
    t = split(dt(i).Header,'] [');
    idx = find(contains(t,'protein_id='));

    if ~isempty(idx)
        dt(i).Header = erase(t{idx},'protein_id=');
    else
        rm(end+1,1) = i;
    end
end

fastawrite('nucleotide_nuc_new.txt',dt);

% Produce fasta file from genbank file 

% if ~exist('gb_cre','var')
%     gb_cre = genbankread('GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.gbff');
% end
% if ~exist('gb_chl','var')
%     gb_chl = genbankread('cre_chloroplast.gb');
% end
if ~exist('gb_mit','var')
    gb_mit = genbankread('cre_mitochondrion.gb');
end

fasta_mit = {};

for i = 1:length(gb_mit)
    for j = 1:length(gb_mit(i).CDS)

        gene = '';
        product = '';
        proteinId = '';
        locusTag = '';
        geneId = '';
        if ~isempty(gb_mit(i).CDS(j).translation)
            seq = gb_mit(i).CDS(j).translation;
        else
            seq = '';
        end

        txt = gb_mit(i).CDS(j).text;
        [hei,~] = size(txt);

        for k = 1:hei
            if contains(txt(k,:),'/gene=')
                gene = erase(txt(k,:),'/gene=');
                gene = erase(gene,' ');
                gene = erase(gene,'"');
            elseif contains(txt(k,:),'/locus_tag=')
                locusTag = erase(txt(k,:),'/locus_tag=');
                locusTag = erase(locusTag,' ');
                locusTag = erase(locusTag,'"');
            elseif contains(txt(k,:),'/product=')
                product = erase(txt(k,:),'/product=');
                product = erase(product,' ');
                product = erase(product,'"');
            elseif contains(txt(k,:),'/protein_id=')
                proteinId = erase(txt(k,:),'/protein_id=');
                proteinId = erase(proteinId,' ');
                proteinId = erase(proteinId,'"');
            elseif contains(txt(k,:),'/db_xref="GeneID:')
                geneId = erase(txt(k,:),'/db_xref="GeneID:');
                geneId = erase(geneId,' ');
                geneId = erase(geneId,'"');
            elseif contains(txt(k,:),'/translation=') && isempty(seq)
                seq = erase(txt(k,:),'/translation=');
                seq = erase(seq,' ');
                seq = erase(seq,'"');
            end
        end

        fasta_mit{end+1,1} = '>';
        fasta_mit{end,1} = [fasta_mit{end,1},locusTag];
        
        if ~isempty(gene)
            fasta_mit{end,1} = [fasta_mit{end,1},'[gene=',gene,']'];
        end
        if ~isempty(product)
            fasta_mit{end,1} = [fasta_mit{end,1},'[product=',product,']'];
        end
        if ~isempty(proteinId)
            fasta_mit{end,1} = [fasta_mit{end,1},'[protein_id=',proteinId,']'];
        end
        if ~isempty(geneId)
            fasta_mit{end,1} = [fasta_mit{end,1},'[gene_id=',geneId,']'];
        end

        fasta_mit{end+1,1} = seq;
    end
end

clear txt proteinId locusTag product hei i j k gene geneId seq;
writecell(fasta_mit,'cre_mitochondrion_genomic.txt','QuoteStrings',false);

% =========================================================================
% fasta_chl = {};
% 
% for i = 1:length(gb_chl)
%     for j = 1:length(gb_chl(i).CDS)
% 
%         gene = '';
%         product = '';
%         proteinId = '';
%         locusTag = '';
%         geneId = '';
%         if ~isempty(gb_chl(i).CDS(j).translation)
%             seq = gb_chl(i).CDS(j).translation;
%         else
%             seq = '';
%         end
% 
%         txt = gb_chl(i).CDS(j).text;
%         [hei,~] = size(txt);
% 
%         for k = 1:hei
%             if contains(txt(k,:),'/gene=')
%                 gene = erase(txt(k,:),'/gene=');
%                 gene = erase(gene,' ');
%                 gene = erase(gene,'"');
%             elseif contains(txt(k,:),'/locus_tag=')
%                 locusTag = erase(txt(k,:),'/locus_tag=');
%                 locusTag = erase(locusTag,' ');
%                 locusTag = erase(locusTag,'"');
%             elseif contains(txt(k,:),'/product=')
%                 product = erase(txt(k,:),'/product=');
%                 product = erase(product,' ');
%                 product = erase(product,'"');
%             elseif contains(txt(k,:),'/protein_id=')
%                 proteinId = erase(txt(k,:),'/protein_id=');
%                 proteinId = erase(proteinId,' ');
%                 proteinId = erase(proteinId,'"');
%             elseif contains(txt(k,:),'/db_xref="GeneID:')
%                 geneId = erase(txt(k,:),'/db_xref="GeneID:');
%                 geneId = erase(geneId,' ');
%                 geneId = erase(geneId,'"');
%             elseif contains(txt(k,:),'/translation=') && isempty(seq)
%                 seq = erase(txt(k,:),'/translation=');
%                 seq = erase(seq,' ');
%                 seq = erase(seq,'"');
%             end
%         end
% 
%         fasta_chl{end+1,1} = '>';
%         fasta_chl{end,1} = [fasta_chl{end,1},locusTag];
%         
%         if ~isempty(gene)
%             fasta_chl{end,1} = [fasta_chl{end,1},'[gene=',gene,']'];
%         end
%         if ~isempty(product)
%             fasta_chl{end,1} = [fasta_chl{end,1},'[product=',product,']'];
%         end
%         if ~isempty(proteinId)
%             fasta_chl{end,1} = [fasta_chl{end,1},'[protein_id=',proteinId,']'];
%         end
%         if ~isempty(geneId)
%             fasta_chl{end,1} = [fasta_chl{end,1},'[gene_id=',geneId,']'];
%         end
% 
%         fasta_chl{end+1,1} = seq;
%     end
% end
% 
% clear txt proteinId locusTag product hei i j k gene geneId seq;
% writecell(fasta_chl,'cre_chloroplast_genomic.txt','QuoteStrings',false);

% =========================================================================
% fasta_cre = {};
% 
% for i = 1:length(gb_cre)
%     for j = 1:length(gb_cre(i).CDS)
% 
%         if isempty(gb_cre(i).CDS(j).translation)
%             continue;
%         end
% 
%         dbRef = '';
%         proteinId = '';
%         locusTag = '';
%         oldLocusTag = '';
%         geneId = '';
% 
%         txt = gb_cre(i).CDS(j).text;
%         [hei,~] = size(txt);
% 
%         for k = 1:hei
%             if contains(txt(k,:),'/db_xref="Phytozome:')
%                 dbRef = erase(txt(k,:),'/db_xref="Phytozome:');
%                 dbRef = erase(dbRef,' ');
%                 dbRef = erase(dbRef,'"');
%             elseif contains(txt(k,:),'/locus_tag=')
%                 locusTag = erase(txt(k,:),'/locus_tag=');
%                 locusTag = erase(locusTag,' ');
%                 locusTag = erase(locusTag,'"');
%             elseif contains(txt(k,:),'/old_locus_tag=')
%                 oldLocusTag = erase(txt(k,:),'/old_locus_tag=');
%                 oldLocusTag = erase(oldLocusTag,' ');
%                 oldLocusTag = erase(oldLocusTag,'"');
%             elseif contains(txt(k,:),'/protein_id=')
%                 proteinId = erase(txt(k,:),'/protein_id=');
%                 proteinId = erase(proteinId,' ');
%                 proteinId = erase(proteinId,'"');
%             elseif contains(txt(k,:),'/db_xref="GeneID:')
%                 geneId = erase(txt(k,:),'/db_xref="GeneID:');
%                 geneId = erase(geneId,' ');
%                 geneId = erase(geneId,'"');
%             end
%         end
% 
%         fasta_cre{end+1,1} = '>';
%         fasta_cre{end,1} = [fasta_cre{end,1},dbRef];
%         
%         if ~isempty(locusTag)
%             fasta_cre{end,1} = [fasta_cre{end,1},'[locus_tag=',locusTag,']'];
%         end
%         if ~isempty(oldLocusTag)
%             fasta_cre{end,1} = [fasta_cre{end,1},'[old_locus_tag=',oldLocusTag,']'];
%         end
%         if ~isempty(proteinId)
%             fasta_cre{end,1} = [fasta_cre{end,1},'[protein_id=',proteinId,']'];
%         end
%         if ~isempty(geneId)
%             fasta_cre{end,1} = [fasta_cre{end,1},'[gene_id=',geneId,']'];
%         end
% 
%         fasta_cre{end+1,1} = gb_cre(i).CDS(j).translation;
%     end
% end
% 
% clear txt proteinId locusTag oldLocusTag hei i j k dbRef geneId;
% writecell(fasta_cre,'Chlamydomonas_reinhardtii_v5.5_genomic.txt');

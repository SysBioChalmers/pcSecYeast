function [model,allrxns,geneError] = addModificationforSec(model,gene,protein_info,onlyrxns)
% onlyrxns = 0,1. 0 stands for that you are going to get a new model. 1
% stands for that you only want to export the rxnID list.

if nargin < 4
    onlyrxns = 0;
end
allrxns = [];
geneError = [];
 [~,geneidx] =ismember(gene,protein_info(:,2));
    if geneidx ~= 0 && cell2mat(protein_info(geneidx,3)) == 1
        [~,i] =ismember(gene,protein_info(:,2));
        peptide = cell2mat(protein_info(i,2));
        seq = cell2mat(protein_info(i,11));
        SP=cell2mat(protein_info(i,4));
        DSB=cell2mat(protein_info(i,5));
        NG=cell2mat(protein_info(i,6));
        GPI =cell2mat(protein_info(i,9));
        OG =cell2mat(protein_info(i,7));
        Trans =cell2mat(protein_info(i,8));
        
        compartment =cell2mat(protein_info(i,10)); %peptide compartment
        S = regexp(compartment, ',', 'split');
        compartment = S(1);
        peptide = strrep(peptide,'-','_');
        % The first step is to translocate the translated peptide into [ER]
        [model,rxns{1}] = translocate(model,peptide,length(seq),SP,NG,DSB,GPI,onlyrxns);
        % Modify the transported peptide with DSB, GPI, NG, and OG
        [model, peptide_name,rxns{2}] = addDSB(model,peptide,length(seq),DSB,onlyrxns);
        [model, peptide_name,rxns{3}] = addGPI(model,peptide_name,length(seq),GPI,onlyrxns);
        [model, peptide_name,rxns{4}] = addOG(model,peptide_name,length(seq),OG,onlyrxns);
        [model, peptide_name,rxns{5}] = addNG(model,peptide_name,length(seq),NG,onlyrxns);
        % add misfolding reactions
        [model, peptide_name,rxns{6}] = addGPIMisfold(model,peptide_name,peptide,length(seq),NG,OG,DSB,GPI,onlyrxns);
        [model,peptide_name,rxns{7}] = addNGMisfold(model,peptide_name,peptide,length(seq),NG,OG,DSB,GPI,Trans,compartment,onlyrxns);
        [model,peptide_name,rxns{8}] = addOtherMisfold(model,peptide_name,peptide,length(seq),NG,OG,DSB,GPI,onlyrxns);
        %another pathway for misfolding degradation(the end peptide should be %s
        
        % transporting the peptide from [ER] to [g]
        [model, peptide_name,rxns{9}] = coat_GPI(model,peptide_name,GPI,onlyrxns);
        [model, peptide_name,rxns{10}] = coat_trans_membrane(model,peptide_name,GPI,Trans,onlyrxns);
        [model, peptide_name,rxns{11}] = coat_other(model,peptide_name,GPI,Trans,onlyrxns);
        
        %procssing in [g]: Golgi processing N
        [model, peptide_name,rxns{12}] = golgiProcessing_N(model,peptide_name,length(seq),NG,onlyrxns);
        [model, peptide_name,rxns{13}] = golgiProcessing_O(model,peptide_name,length(seq),OG,onlyrxns);
        
        % mature
        [model, peptide_name,rxns{14}] = mature(model,peptide_name,onlyrxns);
        
        % trnasport to other final destination(change the finale name to be original peptide_mature)
        [model, peptide_name,rxns{15}] = transportToFinal(model,peptide_name,peptide,compartment,onlyrxns);
        
%         % Subunit transport
%         [model,rxns{16}] = subunitTrans(model,peptide,compartment,onlyrxns);
        
        % get all rxns
        for j = 1:length(rxns)
            if ~isempty(rxns{j})
                allrxns = [allrxns;rxns{j}];
            end
        end
    elseif geneidx == 0
        geneError = gene;
    end
    
end

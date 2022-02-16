function model = addfolding(model,peptide,protein_info)
for k = 1:numel(peptide)
     disp(['Adding folding:' num2str(k) '/' num2str(length(peptide))]);
    %gather all information for peptide
    [~,geneidx] =ismember(peptide(k),protein_info(:,2));
    if geneidx == 0
        error([peptide{k},':cannot find gene PTM information'])
    end
    peptide_comp = protein_info(geneidx,10); %peptide compartment
    % proteins processed in the secretory pathway
    if geneidx ~= 0 && cell2mat(protein_info(geneidx,3)) == 1 
        model = addModificationforSec(model,peptide(k),protein_info);
    else % NOT processed in the secretory pathway
        peptide(k) = strrep(peptide(k),'-','_');
        % start to adding reactions
        if strcmp(peptide_comp,'c')==0
            % for other compartment rather than cytosol
            r=sprintf('%s_peptide[c] => %s_peptide[%s]',cell2mat(peptide(k)),cell2mat(peptide(k)),cell2mat(peptide_comp));
            rxnID=sprintf('%s_importing_%s',cell2mat(peptide(k)),cell2mat(peptide_comp));
            rxnID=strrep(rxnID,'-','');
            rxnName=sprintf('Importing %s into %s',cell2mat(peptide(k)),cell2mat(peptide_comp));
            model=addYeastReaction(model,r,rxnID,{rxnName});
            %subunit export to cytosol
            r=sprintf('%s_subunit[%s] => %s_subunit[c]',cell2mat(peptide(k)), cell2mat(peptide_comp), cell2mat(peptide(k)));
            rxnID=sprintf('%s_subunit_export_%s',cell2mat(peptide(k)),cell2mat(peptide_comp));
            rxnID=strrep(rxnID,'-','');
            rxnName=sprintf('export %s subunit from %s',cell2mat(peptide(k)),cell2mat(peptide_comp));
            model=addYeastReaction(model,r,rxnID,{rxnName});        
        end
        %folding
        r=sprintf('%s_peptide[%s] => %s_folding[%s]',cell2mat(peptide(k)),cell2mat(peptide_comp),cell2mat(peptide(k)),cell2mat(peptide_comp));
        rxnID=sprintf('%s_folding_%s',cell2mat(peptide(k)),cell2mat(peptide_comp));
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('%s folding into %s',cell2mat(peptide(k)),cell2mat(peptide_comp));
        model=addYeastReaction(model,r,rxnID,{rxnName});
        
        %misfolding
        r=sprintf('%s_folding[%s] => %s_misfolding[%s]',cell2mat(peptide(k)),cell2mat(peptide_comp),cell2mat(peptide(k)),cell2mat(peptide_comp));
        rxnID=sprintf('%s_misfold_%s',cell2mat(peptide(k)),cell2mat(peptide_comp));
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('%s Misfolding',cell2mat(peptide(k)));
        model=addYeastReaction(model,r,rxnID,{rxnName});
        
        %refolding is not considered for now but can be added if needed in
        % the future
%         r=sprintf('%s_misfolding[%s] => %s_folding[%s]',cell2mat(peptide(k)),cell2mat(peptide_comp),cell2mat(peptide(k)),cell2mat(peptide_comp));
%         rxnID=sprintf('%s_refolding_%s',cell2mat(peptide(k)),cell2mat(peptide_comp));
%         rxnID=strrep(rxnID,'-','');
%         rxnName=sprintf('%s refolding',cell2mat(peptide(k)));
%         model=addYeastReaction(model,r,rxnID,{rxnName});
        
        %misfoldeing degradation
        r=sprintf('%s_misfolding[%s] => %s_subunit[%s]',cell2mat(peptide(k)),cell2mat(peptide_comp),cell2mat(peptide(k)),cell2mat(peptide_comp));
        rxnID=sprintf('%s_degradation_misfolding_%s',cell2mat(peptide(k)),cell2mat(peptide_comp));
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('%s Misfolding',cell2mat(peptide(k)));
        model=addYeastReaction(model,r,rxnID,{rxnName});
        
        %misfolding dilution
        r=sprintf('%s_misfolding[%s] => ',cell2mat(peptide(k)),cell2mat(peptide_comp));
        rxnID=sprintf('%s_dilution_misfolding_%s',cell2mat(peptide(k)),cell2mat(peptide_comp));
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('%s Misfolding dilution',cell2mat(peptide(k)));
        model=addYeastReaction(model,r,rxnID,{rxnName});
        
        % add subunit transport
         model = subunitTrans(model,cell2mat(peptide(k)),peptide_comp);

        
    end
end
end
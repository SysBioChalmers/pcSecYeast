function [model,enzymedata_TP] =addTargetProtein(model,TP,new,fakeProteinInfo)
% new means create a new rxn instead of merging with the original protein
% synthesis rxn
if nargin < 3
    new = false;
    fakeProteinInfo = [];
end
if nargin < 4
    fakeProteinInfo = [];
end
% this function is to add
load('Protein_Sequence.mat');
[~,~,protein_info] = xlsread('TargetProtein.xlsx','protein_info');
TPall.seq = protein_info(2:end,11);
TPall.id = protein_info(2:end,2);
complex_info = protein_info(2:end,[1,10,15:16]);
protein_info = protein_info(2:end,1:14);
[~,~,protein_info2] = xlsread('Protein_Information.xlsx');
protein_info2 = protein_info2(2:end,1:14);
protein_info2(:,1) = protein_info2(:,2);

if new % means that new protein is added with prefix new even thought it is a native protein
    protein_info2(:,1) = strcat(protein_info2(:,1),'new');
    protein_info2(:,2) = strcat(protein_info2(:,2),'new');
    ProteinSequence.id = strcat(ProteinSequence.id,'new');
    protein_info2(:,1) = strrep(protein_info2(:,1),'-','');
    protein_info2(:,2) = strrep(protein_info2(:,2),'-','');
    ProteinSequence.id = strrep(ProteinSequence.id,'-','');
    
end
protein_info = [protein_info;protein_info2];

ProteinSequence.id(end+1:end+length(TPall.id)) = TPall.id;
ProteinSequence.seq(end+1:end+length(TPall.id)) = TPall.seq;
ProteinSequence.fullseq(end+1:end+length(TPall.id)) = TPall.seq;

if ~isempty(fakeProteinInfo)
    protein_info = [protein_info;fakeProteinInfo];
    ProteinSequence.id(end+1:end+length(fakeProteinInfo(:,2))) = fakeProteinInfo(:,2);
    ProteinSequence.seq(end+1:end+length(fakeProteinInfo(:,2))) = fakeProteinInfo(:,11);
    ProteinSequence.fullseq(end+1:end+length(fakeProteinInfo(:,2))) = fakeProteinInfo(:,11);
end

Target = protein_info(ismember(protein_info(:,1),TP),2); % find all proteins in the complex
model = addfolding(model,Target,protein_info);% add folding process for those peptide xx_peptide --> xx_folding
model = addTranslationRxns(model,ProteinSequence,Target);
model = addDegradationRxns(model,ProteinSequence,protein_info,Target);
enzymedata_TP.proteins = Target;
enzymedata_TP.kdeg(1:length(Target),1) = 0;
%enzymedata_TP.enzyme = TP;


% add a exchange rxn for the complex or the TP and add acomplex formation reaction for the complex
for i = 1:length(TP)
    if endsWith(TP(i),'complex')
        [~,geneidx] =ismember(TP(i),complex_info(:,1));
        if geneidx ~= 0
            sub = complex_info{geneidx,4};
            peptide_comp = complex_info(geneidx,2); %peptide compartment
            cmplxid = strcat(cell2mat(TP(i)),'_folding[',peptide_comp,']');
            reaction.eq =  sprintf('%s => %s',sub,cell2mat(cmplxid));
            reaction.rxnName = strcat(TP(i),'_complex_formation');
            reaction.rxns =strcat(TP(i),'_complex_formation');
            model=addYeastReaction(model,reaction.eq,reaction.rxns,reaction.rxnName);
            
%             if ~strcmp(peptide_comp,'e')
%                 metlist_dil = cmplxid;
%                 coeflist_dil = -1;
%                 rxnid_dil = strcat(cell2mat(TP(i)),'_complex_dilution');
%                 model = addReaction(model,rxnid_dil,'metaboliteList',metlist_dil,'stoichCoeffList',coeflist_dil,'reversible',false);
%             else
                disp(['adding exchange rxn for the protein:',num2str(i)])
                rxnid = strcat(cell2mat(TP(i)),' exchange');
                model = addReaction(model,rxnid,'metaboliteList',cmplxid,'stoichCoeffList',-1,'reversible',false);
%            end
        else
            warning('no complex info')
        end
    else
        [~,geneidx] =ismember(TP(i),protein_info(:,2));
        peptide_comp = protein_info(geneidx,10); %peptide compartment
        cmplxid = strcat(cell2mat(TP(i)),'_folding[',peptide_comp,']');
        disp(['adding exchange rxn for the protein:',num2str(i)])
%          if ~strcmp(peptide_comp,'e')
%              metlist_dil = cmplxid;
%              coeflist_dil = -1;
%              rxnid_dil = strcat(cell2mat(TP(i)),'_complex_dilution');
%              model = addReaction(model,rxnid_dil,'metaboliteList',metlist_dil,'stoichCoeffList',coeflist_dil,'reversible',false);
%          else
            rxnid = strcat(cell2mat(TP(i)),' exchange');
            model = addReaction(model,rxnid,'metaboliteList',cmplxid,'stoichCoeffList',-1,'reversible',false);
%         end
    end
    
    
    model.id = [TP{i},'pcSecYeast model'];
    
    % find complex formation
%     enzfmtrxn_id = strcat(cell2mat(TP(i)),'_complex_formation');
%     idx_tmp = ismember(model.rxns,enzfmtrxn_id);
%     if any(idx_tmp)
%         s_tmp = model.S(:,idx_tmp);
%         subunits_tmp = model.mets(s_tmp < 0);
%         na_tmp = repelem({''},max_subunit-length(subunits_tmp));
%         subunits_tmp = cellfun(@(x) strrep(x,'_folding',''),subunits_tmp,'UniformOutput',false);
%         subunits_tmp = cellfun(@(x) x(1:strfind(x,'[')-1),subunits_tmp,'UniformOutput',false);
%         subunits_tmp = cellfun(@(x) strrep(x,'_','-'),subunits_tmp,'UniformOutput',false);
%         enzymedata_TP.subunit(i,:) = [subunits_tmp' na_tmp];
%         
%         % add stoichiometry of subunits
%         stoichi_tmp = abs(full(s_tmp(s_tmp < 0)));
%         na_tmp = repelem(0,max_subunit-length(stoichi_tmp));
%         enzymedata_TP.subunit_stoichiometry(i,:) = [stoichi_tmp' na_tmp];
%         
%     else
%         enzymedata_TP.subunit(i,1) = TP(i);
%         enzymedata_TP.subunit_stoichiometry(i,1) = 1;
%     end
 end

enzymedata_TP = calculateMW(enzymedata_TP,ProteinSequence,protein_info);

end



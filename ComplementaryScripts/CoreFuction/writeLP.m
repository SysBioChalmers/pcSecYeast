%% writeLP
function fileName = writeLP(model,mu,f,f_mito,funmodelER,f_erm,osenseStr,rxnID,enzymedata,factor_k,name)


% f is the fraction (g/gCDW) of the modeled proteins.
% f_mito is the fraction (g/gCDW) of the mitochondrial proteins.
if exist('factor_k', 'var')
    if isempty(factor_k)
        factor_k = 1;
    end
else
    factor_k = 1;
end
model = changeRxnBounds(model,'r_2111',mu);
% if factor_k_withoutmisfolding == 0
%     model.ub(contains(model.rxns,'misfolding')) = 0;
%     model.ub(contains(model.rxns,'misfolding')) = 0;
%     model.ub(contains(model.rxns,'misfolding')) = 0;
%     model.lb(contains(model.rxns,'refolding')) = 0;
% end % it is also not needed to add this but would lead to different results if the precision of the solver is not high.


fileName = sprintf(['Simulation_dilution',name,'.lp']);
fptr = fopen(fileName,'w');

% Set objective function
osenseStr = strcat(osenseStr,'\n');
fprintf(fptr,osenseStr);
index_obj = find(ismember(model.rxns,rxnID));
fprintf(fptr,'obj: X%d\n',index_obj);
fprintf(fptr,'Subject To\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) SV = 0.
% (From Ibrahim)
for i = 1:numel(model.mets)
    j = find(full(model.S(i,:)));
    for m = 1:numel(j)
        s = full(model.S(i,j(m)));
        if mod(m,200) == 0
            sep = newline;
        else
            sep = '';
        end
        if m == 1
           eq = sprintf('%.15f X%d',s,j(m));
        else
           if s>0
               eq = sprintf('%s + %.15f X%d%c',eq,s,j(m),sep);
           else
               eq = sprintf('%s %.15f X%d%c',eq,s,j(m),sep);
           end
        end
    end
    fprintf(fptr,'C%d: %s = 0\n',i,eq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Coupling metabolic reactions and enzymes.

met_enzyme_idx = find(strcmp(enzymedata.label,'Met'));
enzyme_list = enzymedata.enzyme(met_enzyme_idx);

kcat_list = enzymedata.kcat(met_enzyme_idx);

% % newly added pathway reactions
% newidx = contains(enzyme_list,'new_r_');
% expandidx = contains(enzyme_list,'withoutcofactor');
% lowkcat = quantile(kcat_list(~(newidx|expandidx)),0.01,1);
lowkcat = 3600;

for i=1:numel(enzyme_list)
    enzyme = enzyme_list{i};
  	kcat = kcat_list(i);

%     % Change kcats extremely low for original enzymes
%     if ismember(enzyme,enzyme_list(~newidx))
        if kcat < lowkcat
            kcat = lowkcat;
        end
%     end


    %find enzyme formation reaction id
    id_syn = strcat(enzyme,'_formation');
    idx_syn = find(strcmp(model.rxns,id_syn));

    %calculate the coefficient of enzyme formation reaction rate
    coef = kcat/mu;

    %find enzyme used reaction id (metabolic reaction)
    idx_rxn = find(ismember(model.rxns,strrep(enzyme,'_complex','')));

%     fprintf(fptr,'CM%d: X%d - %.15f X%d <= 0\n',...
%                  i,idx_rxn,coef,idx_syn);
    fprintf(fptr,'CM%d: X%d - %.15f X%d = 0\n',...
                 i,idx_rxn,coef,idx_syn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) constrian Secretory part
%kcat_fake = SimulateSecParam(model,mu);
k1=0;
sec_enzyme_idx = find(strcmp(enzymedata.label,'Sec'));
 
complex_list = enzymedata.enzyme(sec_enzyme_idx);
compartment = enzymedata.comp(sec_enzyme_idx);
M_Kcats_modified = enzymedata.kcat(sec_enzyme_idx);% %%%%%remove later
[~,idx]=unique(strcat(complex_list,compartment),'rows');
complex_list=complex_list(idx,:);
compartment=compartment(idx,:);
M_Kcats_modified = M_Kcats_modified(idx,:);

% get all coef for all rxns coef
rxns_coef = enzymedata.rxnscoef;
rxns_coef_id = enzymedata.rxns;
for i=1:numel(complex_list) %check again for 1:n
    k1=k1+1;
    %coupling rxn formation V <= kcat [e]
    complex_name=cell2mat(complex_list(i));
    rxnID=sprintf('%s_formation',complex_name);
    %find enzyme formation reaction id
    id_syn=find(ismember(model.rxns,rxnID));
    if isempty(id_syn)
        warning([rxnID,' is missing',num2str(i)])
    end
    %add the condition collecting collect all catalyzed reactions by this complex
    %find only the metabolic index by setting the rxn index
    %before the index of the complex formation rxn
    rxns_by_this_complex = find(endsWith(model.rxns,['_',complex_list{i}]));
    [~,idx_coef] = ismember(model.rxns(rxns_by_this_complex),rxns_coef_id);
    coef_tmp = rxns_coef(idx_coef);
    %print kcat constriants in lp file
    for r=1:numel(rxns_by_this_complex)
        sep1='';
        if r == 1
            rxns=sprintf('%.15f X%d',coef_tmp(1),rxns_by_this_complex(1));
        else
            if mod(r,150)==0
                sep1=char(10);
            else
                sep1='';
            end
            rxns=sprintf('%s + %.15f X%d%s',rxns, coef_tmp(r),rxns_by_this_complex(r),sep1);
        end
    end

    newKcat = M_Kcats_modified(i,1);
    coef = 1/mu*newKcat;
    fprintf(fptr,'CSEC%d: %s - %.15f X%d = 0\n',k1,rxns,coef,id_syn); %K for average kcat or KK for different kcats, KK1 for GEKO
end

%% 4) Constraint on total enzymes and dummy complex
% total enzymes
dil_rxns = model.rxns(contains(model.rxns,'_dilution') & ~contains(model.rxns,'dummy'));
for i = 1:length(dil_rxns)
    rxn_id = dil_rxns{i};
    comp_name = strrep(rxn_id,'_dilution','');
    idx = find(strcmp(model.rxns,rxn_id));
    
    if mod(i,200) == 0
        sep = newline;
	else
        sep = '';
    end

	MW = enzymedata.enzyme_MW(contains(enzymedata.enzyme,comp_name));
    if isempty(MW)
        MW = enzymedata.enzyme_MW(contains(enzymedata.enzyme,strrep(comp_name,'_complex','')));
    end
     if isempty(MW)
        MW = enzymedata.proteinMWs(contains(enzymedata.proteins,strrep(extractBefore(comp_name,'_misfolding'),'_','-')));
    end
	coeff = MW(1)/1000;
	if i == 1
        eq = sprintf('%.15f X%d',coeff,idx);
	else
        eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
	end
end

% dummy complex
idx = find(strcmp(model.rxns,'dilute_dummy'));
eq = sprintf('%s + %.15f X%d%c',eq,46000/1000,idx,sep); %460 is MW of dummy complex (g/mol)

% dummy complex
idx = find(strcmp(model.rxns,'dilute_dummyER'));
eq = sprintf('%s + %.15f X%d%c',eq,46000/1000,idx,sep); %460 is MW of dummy complex (g/mol)

fprintf(fptr,'Ctotprot: %s = %.15f\n',eq,mu*f); % unmodeled part in the fraction eliminate only the unmodeled protein part in the biomass equation

%unmodeled ER protein
fprintf(fptr,'CERunmodel: %.15f X%d = %.15f\n',46000/1000,idx,mu*funmodelER);

% %kdeg = 0.042; % around 30% nacent protein are degradated 
% trans_rxns = model.rxns(endsWith(model.rxns,'_translation'));
% for i = 1:length(trans_rxns)
%     rxn_id = trans_rxns{i};
%     comp_name = strrep(rxn_id,'_peptide_translation','');
%     comp_name = strrep(comp_name,'r_','');
%     idx = find(strcmp(model.rxns,rxn_id));
%     if mod(i,150) == 0
%         sep = newline;
% 	else
%         sep = '';
%     end
%     kdeg = enzymedata.kdeg(ismember(enzymedata.proteins,strrep(comp_name,'_','-')));
% 	%kdeg = 0.042;
%     MW = enzymedata_all.proteinMWs(ismember(enzymedata_all.proteins,strrep(comp_name,'_','-')));
% 	coeff = MW(1)/1000/mu;
% 	if i == 1
%         eq = sprintf('%.15f X%d',coeff,idx);
% 	else
%         eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
% 	end
% end

% % dummy complex
% idx = find(strcmp(model.rxns,'dilute_dummy'));
% eq = sprintf('%s + %.15f X%d%c',eq,46000/1000*mu,idx,sep); % 46000 is MW of dummy complex (g/mol)  the same as the protein in the biomass 
% 
% fprintf(fptr,'Ctotprot: %s = %.15f\n',eq,f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5) Constraint on mitochondrial proteins.
[~, mitoRxns] = collectCompartment(model,{'m','mm'});
mito_dil_rxns = strcat(mitoRxns,'_complex_dilution');

for i = 1:length(mito_dil_rxns)
    rxn_id = mito_dil_rxns{i};
    comp_name = strrep(rxn_id,'_dilution','');
    idx = find(strcmp(model.rxns,rxn_id));

    if mod(i,200) == 0
        sep = newline;
	else
        sep = '';
    end

	MW = enzymedata.enzyme_MW(ismember(enzymedata.enzyme,comp_name));
	coeff = MW/1000;
	if i == 1
        eq = sprintf('%.15f X%d',coeff,idx);
	else
        eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
	end
end

fprintf(fptr,'Cmito: %s <= %.15f\n',eq,mu*f_mito);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6) Constraint on proteasome.

rxnID=sprintf('Mach_proteasome_complex_formation');
id_syn=find(ismember(model.rxns,rxnID));

%calculate the coefficient of enzyme formation reaction rate
kcat_Proteasome = enzymedata.kcat(ismember(enzymedata.enzyme,'Mach_proteasome_complex'));
coef_syn = kcat_Proteasome/mu;

deg_rxns = model.rxns(endsWith(model.rxns,'_subunit_degradation') | endsWith(model.rxns,'_sp_degradation'));
for i = 1:length(deg_rxns)
    rxn_id = deg_rxns{i};
    if contains(rxn_id,'_subunit_degradation')
        comp_name = strrep(rxn_id,'_subunit_degradation','');
        comp_name = strrep(comp_name,'r_','');
        prot_leng = enzymedata.proteinLength(ismember(enzymedata.proteins,strrep(comp_name,'_','-')));
        coef = prot_leng(1)/467;
    else
        comp_name = strrep(rxn_id,'_sp_degradation','');
        comp_name = strrep(comp_name,'r_','');
        coef = 25/467;
    end

    idx = find(strcmp(model.rxns,rxn_id));

    if mod(i,150) == 0
        sep = newline;
	else
        sep = '';
    end

	if i == 1
        eq = sprintf('%.15f X%d',coef,idx);
	else
        eq = sprintf('%s + %.15f X%d%c',eq,coef,idx,sep);
	end
end
fprintf(fptr,'CPD: %s - %.15f X%d = 0\n',eq,coef_syn,id_syn);


%% 7) Constraint on Ribosome.
rxnID=sprintf('Mach_Ribosome_complex_formation');
id_syn=find(ismember(model.rxns,rxnID));
%kcat_ribo = enzymedataMachine.kcat(ismember(enzymedataMachine.enzyme,'Ribosome'));
kcat_ribo = 3.96E+01*mu/(0.885620915+mu)*3600;

coef = kcat_ribo/mu;

trans_rxns = model.rxns(endsWith(model.rxns,'_translation'));
for i = 1:length(trans_rxns)
    rxn_id = trans_rxns{i};
    comp_name = strrep(rxn_id,'_peptide_translation','');
    comp_name = strrep(comp_name,'r_','');
    idx = find(strcmp(model.rxns,rxn_id));
    prot_leng = enzymedata.proteinLength(ismember(enzymedata.proteins,strrep(comp_name,'_','-')));

    if mod(i,150) == 0
        sep = newline;
	else
        sep = '';
    end

	if i == 1
        eq = sprintf('%d X%d',prot_leng(1),idx);
	else
        eq = sprintf('%s + %d X%d%c',eq,prot_leng(1),idx,sep);
	end
end
idx = find(strcmp(model.rxns,'translate_dummy')); % dummyER has been included by using the rxnID: 
eq = sprintf('%s + %.15f X%d%c',eq,423,idx,sep); % 460 is MW of dummy complex (g/mol)  the same as the protein in the biomass,4.23 is the protein_length 

idx = find(strcmp(model.rxns,'translate_dummyER')); % dummyER has been included by using the rxnID: 
eq = sprintf('%s + %.15f X%d%c',eq,423,idx,sep);

idx = find(strcmp(model.rxnNames,'protein pseudoreaction'));
eq = sprintf('%s + %.15f X%d%c',eq,4.23,idx,sep); % the original 460/4.23

fprintf(fptr,'Cribo: %s - %.15f X%d = 0\n',eq,coef,id_syn);

%% 8) Constraint on Ribosome assembly.
kcat_assembly = enzymedata.kcat(ismember(enzymedata.enzyme,'Mach_Ribosome_Assembly_Factors_complex'));
coef = kcat_assembly/mu;
rxnID=sprintf('Mach_Ribosome_complex_formation');
id_use=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Mach_Ribosome_Assembly_Factors_complex_formation');
id_syn=find(ismember(model.rxns,rxnID));

fprintf(fptr,'Criboassem: X%d - %.15f X%d = 0\n',id_use,coef,id_syn);


%% 9) Constraint on misfolding.

misfold_list = model.rxns(contains(model.rxns,'_misfold_')& ~contains(model.rxns,'dummy'));
for i = 1:length(misfold_list)
    rxn_id = misfold_list{i};
    comp_name = extractBefore(rxn_id,'_misfold_');
    comp_id = ['r_',comp_name,'_peptide_translation'];
    kdegratio = enzymedata.kdeg(ismember(enzymedata.proteins,strrep(comp_name,'_','-')));
    %kdeg = 0.042;
    idx_syn = find(strcmp(model.rxns,comp_id));
    idx_misfold = find(strcmp(model.rxns,rxn_id));
    %coef_misfold = kdeg/(0.36 + kdeg); % bionumber 113409
    coef_misfold = kdegratio;
    coef = coef_misfold;
    fprintf(fptr,'Cmisfold%d: X%d - %.15f X%d = 0\n',i,idx_misfold,coef,idx_syn);

end

%% er machinery protein
dil_rxns = model.rxns(contains(model.rxns,'complex_dilution') & contains(model.rxns,enzymedata.enzyme(sec_enzyme_idx)));
for i = 1:length(dil_rxns)
    rxn_id = dil_rxns{i};
    comp_name = strrep(rxn_id,'_dilution','');
    idx = find(strcmp(model.rxns,rxn_id));
    
    if mod(i,200) == 0
        sep = newline;
	else
        sep = '';
    end

	MW = enzymedata.enzyme_MW(contains(enzymedata.enzyme,comp_name));
    if isempty(MW)
        MW = enzymedata.enzyme_MW(contains(enzymedata.enzyme,strrep(comp_name,'_complex','')));
    end
	coeff = MW(1)/1000;
	if i == 1
        eq = sprintf('%.15f X%d',coeff,idx);
	else
        eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
	end
end

fprintf(fptr,'CtotERMachin: %s <= %.15f\n',eq,mu*0.022); % unmodeled part in the fraction eliminate only the unmodeled protein part in the biomass equation


%% erm constaint

ermprotein = enzymedata.proteins(strcmp(enzymedata.proteinLoc,'erm'));
ermprotein = unique(ermprotein);

[~,idx1] = ismember(ermprotein,enzymedata.proteins);
complexidx1 = enzymedata.enzyme(find(sum(enzymedata.enzymeSubunitMatrix(:,idx1),2)));
[~,idx2] = ismember(complexidx1,enzymedata.enzyme);
matrix = enzymedata.enzymeSubunitMatrix(idx2,idx1);
matrix = (matrix.*enzymedata.proteinMWs(idx1)')./1000;
for i = 1:length(complexidx1)
    complex = complexidx1(i);
    rxn_id = strcat(complex,'_dilution');
    %rxn_id = strrep(rxn_id,'_complex_complex_','_complex_'); % fix the naming issue introduced by the sec complex naming
    idx = find(strcmp(model.rxns,rxn_id));
        if mod(i,50) == 0
            sep = newline;
        else
            sep = '';
        end
        
        coeff = sum(matrix(i,:));
        if i == 1
            eq = sprintf('%.15f X%d',coeff,idx);
        else
            eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
        end
end

fprintf(fptr,'CtotERM: %s <= %.15f\n',eq,mu*f_erm); % 3% from the proteome data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set lower and upper bounds.

fprintf(fptr,'Bounds\n');

for i = 1:length(model.rxns)
	if model.ub(i) >= 100
        fprintf(fptr,'%f <= X%d <= +infinity\n',model.lb(i),i);
    else
        fprintf(fptr,'%f <= X%d <= %f\n',model.lb(i),i,model.ub(i));
	end
end

fprintf(fptr,'End\n');
fclose(fptr);

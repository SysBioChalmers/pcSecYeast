function [constraints,constraints_label] = activeConstraint(model,enzymedata_all,fluxes)
% This function is to calculate which constraints is active.

% calculate the protein abundance
constraints = zeros(5,length(fluxes(1,:)));
constraints_label = {'proteome','ER proteins','ERmembrane','HXTs','mito'};
mu = fluxes(strcmp(model.rxns,'r_2111'),:);

%% total protein
% dummy complex
idx = strcmp(model.rxns,'dilute_dummy');
constraints(1,fluxes(idx,:) == 0) = 1; % unmodeled part in the fraction eliminate only the unmodeled protein part in the biomass equation

%% er machinery protein
sec_enzyme_idx = find(strcmp(enzymedata_all.label,'Sec'));

dil_rxns = model.rxns(contains(model.rxns,'complex_dilution') & contains(model.rxns,enzymedata_all.enzyme(sec_enzyme_idx)));
for i = 1:length(dil_rxns)
    rxn_id = dil_rxns{i};
    comp_name = strrep(rxn_id,'_dilution','');
    idx = find(strcmp(model.rxns,rxn_id));

	MW = enzymedata_all.enzyme_MW(contains(enzymedata_all.enzyme,comp_name));
    if isempty(MW)
        MW = enzymedata_all.enzyme_MW(contains(enzymedata_all.enzyme,strrep(comp_name,'_complex','')));
    end
	coeff = MW(1)/1000;
	ER_machine_abun(i,:) =  coeff.*fluxes(idx,:);
end
constraints(2,sum(ER_machine_abun) - mu*0.0174 >= 0) = 1; % unmodeled part in the fraction eliminate only the unmodeled protein part in the biomass equation



%% erm constaint

ermprotein = enzymedata_all.proteins(strcmp(enzymedata_all.proteinLoc,'erm'));
ermprotein = unique(ermprotein);

[~,idx1] = ismember(ermprotein,enzymedata_all.proteins);
complexidx1 = enzymedata_all.enzyme(find(sum(enzymedata_all.enzymeSubunitMatrix(:,idx1),2))); % index the complex which contains those erm protein as subunits
[~,idx2] = ismember(complexidx1,enzymedata_all.enzyme);
matrix = enzymedata_all.enzymeSubunitMatrix(idx2,idx1);
matrix = (matrix.*enzymedata_all.proteinMWs(idx1)')./1000;
for i = 1:length(complexidx1)
    complex = complexidx1(i);
    rxn_id = strcat(complex,'_dilution');
    %rxn_id = strrep(rxn_id,'_complex_complex_','_complex_'); % fix the naming issue introduced by the sec complex naming
    idx = find(strcmp(model.rxns,rxn_id));
     
   coeff = sum(matrix(i,:));
	ER_membran(i,:) =  coeff.*fluxes(idx,:);
end
constraints(3,sum(ER_membran)- mu*0.008 >= 0) = 1; % 3% from the proteome data

%% hxt constaint

HXTs = {'YHR094C';'YFL011W';'YOL156W';'YEL069C';'YNL318C';'YDL245C';'YJR158W';'YNR072W';'YMR011W';'YDR345C';'YHR092C';'YHR096C';'YDR343C';'YDR342C';'YJL214W';'YJL219W'};
% 'YIL170W' is not in the model
[~,idx1] = ismember(HXTs,enzymedata_all.proteins);
complexidx1 = enzymedata_all.enzyme(find(sum(enzymedata_all.enzymeSubunitMatrix(:,idx1),2))); % index the complex which contains those erm protein as subunits
[~,idx2] = ismember(complexidx1,enzymedata_all.enzyme);
matrix = enzymedata_all.enzymeSubunitMatrix(idx2,idx1);
matrix = (matrix.*enzymedata_all.proteinMWs(idx1)')./1000;
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
     HXT_abun(i,:) =  coeff.*fluxes(idx,:);
end
constraints(4,sum(HXT_abun)- mu*0.0042 >= 0) = 1; % 3% from the proteome data


%% mito
[~, mitoRxns] = collectCompartment(model,{'m','mm'});
mito_dil_rxns = strcat(mitoRxns,'_complex_dilution');

for i = 1:length(mito_dil_rxns)
    rxn_id = mito_dil_rxns{i};
    comp_name = strrep(rxn_id,'_dilution','');
    idx = find(strcmp(model.rxns,rxn_id));

	MW = enzymedata_all.enzyme_MW(ismember(enzymedata_all.enzyme,comp_name));
	coeff = MW/1000;
	mito_abun(i,:) =  coeff.*fluxes(idx,:);
end
f_mito = 0.05;
constraints(5,sum(mito_abun)- mu*f_mito >= 0) = 1; % from proteome abun


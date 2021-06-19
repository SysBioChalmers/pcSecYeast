function [results] = quickcheckComplexFormation(model)
%this function is to quick check whether all complexes can be synthesised
%in the model


complex_list = model.rxns(contains(model.rxns,'complex_formation'));
%quick check for the complex formation
 results.formation = cell(length(complex_list),2); 
 for i = 1:length(complex_list)
     model_tmp = changeObjective(model,complex_list(i),1);
     sol=optimizeCbModel(model_tmp,'max');
     results.formation(i,:) = [complex_list(i),num2cell(sol.f)];
 end
 
 %quick check for the complex degradation
 results.degradation = cell(length(complex_list),2);
 complex_list = model.rxns(contains(model.rxns,'degradation_misfolding'));
  for i = 1:length(complex_list)
     model = changeObjective(model,complex_list(i),1);
     sol=optimizeCbModel(model,'max');
     results.degradation(i,:) = [complex_list(i),num2cell(sol.f)];
 end
 
  %quick check for the complex dilution
 results.dilution = cell(length(complex_list),2);
 complex_list = model.rxns(contains(model.rxns,'complex_dilution'));
  for i = 1:length(complex_list)
     model = changeObjective(model,complex_list(i),1);
     sol=optimizeCbModel(model,'max');
     results.dilution(i,:) = [complex_list(i),num2cell(sol.f)];
 end
 
 %check protein translation reaction 
 results.translation = cell(length(complex_list),2);
 complex_list = model.rxns(contains(model.rxns,'_translation'));
  for i = 1:length(complex_list)
     model = changeObjective(model,complex_list(i),1);
     sol=optimizeCbModel(model,'max');
     results.translation(i,:) = [complex_list(i),num2cell(sol.f)];
     i
 end


function [newModel,rxns] = translocate(model,peptide,seq,SP,NG,DSB,GPI,onlyrxns)
rxns = [];
if SP==1
    Reaction=co_Translation_translocation(peptide);
    for i=1:6
        if onlyrxns == 1
            rxns = [rxns;{Reaction{i}.rxns}];
        else
            model=addYeastReaction(model,Reaction{i}.eq,{Reaction{i}.rxns},{Reaction{i}.rxnNames});
            rxns = [rxns;{Reaction{i}.rxns}];
        end
    end
elseif GPI>0 && SP == 0
    Reaction=Post_Translation_translocation_tail(peptide);
    for i=1:3
        if onlyrxns == 1
            rxns = [rxns;{Reaction{i}.rxns}];
        else
            model=addYeastReaction(model,Reaction{i}.eq,{Reaction{i}.rxns},{Reaction{i}.rxnNames});
            rxns = [rxns;{Reaction{i}.rxns}];
        end
    end      
else  
    Reaction=Post_Translation_translocation(peptide,round(seq/40));
    for i=1:4
        if onlyrxns == 1
            rxns = [rxns;{Reaction{i}.rxns}];
        else
            model=addYeastReaction(model,Reaction{i}.eq,{Reaction{i}.rxns},{Reaction{i}.rxnNames});
            rxns = [rxns;{Reaction{i}.rxns}];
        end
    end

end
newModel=model;

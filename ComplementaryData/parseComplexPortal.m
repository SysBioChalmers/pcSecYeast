[~, ~, raw_complex] = xlsread('Protein_annotation_uniprot.xlsx','complex_portal');
[~, ~, raw_anno] = xlsread('Protein_annotation_uniprot.xlsx','uniprot');

for i = 2:length(raw_complex(:,1)) % first is the head
    proteins = split(raw_complex(i,19),'|');
    proteins = extractBefore(proteins,'(');
    [~,idx] = ismember(proteins,raw_anno(:,2));
    if all(idx)
    complex(i,1) = join(raw_anno(idx,1),'|');
    else
        tmp = proteins(idx==0);
        tmp = extractBefore(tmp,'-');
        proteins(idx == 0) = tmp;
         if all(idx)
        [~,idx] = ismember(proteins,raw_anno(:,2));
            complex(i,1) = join(raw_anno(idx,1),'|');
         else
             warning(['check complex manually: ',raw_complex{i,1}])
         end
    end
end
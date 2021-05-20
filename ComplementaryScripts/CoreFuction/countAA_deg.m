function [sum,energy] = countAA_deg(seq,aa_list,e_list,seq_full)
 % seqfull should be true if it is the full sequence but should be false if
 % it is SP

if iscell(seq)
    seq = cell2mat(seq);
end

sum = struct();
sum.prod = aa_list.prod_deg;
sum.abbr = aa_list.aa;
sum.num = [];

for i = 1:length(sum.abbr)
    AA_abbr = sum.abbr{i};
    sum.num(i,1) = length(strfind(seq,AA_abbr));
end

if seq_full
    % add the first met
    [~,b] = ismember('M',sum.abbr);
    sum.num(b) = sum.num(b) + 1;
end

energy = struct();
seqlen = length(seq);
energy.subs_deg = e_list.subs(1:2);
energy.prod_deg = e_list.prod([1,2,4]);
energy.subs_degnum = repmat(floor(1.3*seqlen),2,1);
energy.prod_degnum = repmat(floor(1.3*seqlen),3,1);

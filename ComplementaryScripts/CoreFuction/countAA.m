%% countAA

function [sum,energy] = countAA(seq,aa_list,e_list)

if iscell(seq)
    seq = cell2mat(seq);
end

sum = struct();
sum.subs = aa_list.subs;
sum.prod = aa_list.prod;
sum.abbr = aa_list.aa;
sum.num = [];

for i = 1:length(sum.abbr)
    AA_abbr = sum.abbr{i};
    sum.num(i,1) = length(strfind(seq,AA_abbr));
end

% add the first met
[~,b] = ismember('M',sum.abbr);
sum.num(b) = sum.num(b) + 1;

% energy
energy = struct();
seqlen = length(seq);
total_atp = 2*seqlen + 2;
total_gtp = 2*seqlen + 4;
total_h2o = total_atp + total_gtp;
total_h = total_atp + total_gtp;
total_pi = total_atp + total_gtp; 
energy.subs = e_list.subs;
energy.prod = e_list.prod;
energy.subnum = [total_h2o total_atp total_gtp]';
energy.prodnum = [total_h total_atp total_gtp total_pi]';

end

%% Extract Protein Sequence

file = 'orf_trans_all_R64-2-1_20150113.fasta';

raw = importdata(file);
n = length(raw);

ProteinSequence = struct();
ProteinSequence.id = cell(n,1);
ProteinSequence.seq = cell(n,1);
ProteinSequence.fullseq = cell(n,1);

for i = 1:n
    pair = raw(i);
    tmp = strfind(pair.Header,' ');
    genename = pair.Header(1:tmp(1)-1);
    sequence = pair.Sequence;
    sequence = strrep(sequence,'*','');
    ProteinSequence.id{i} = genename;
    ProteinSequence.fullseq{i} = sequence;
    reduced_seq = sequence; % Used to remove N-terminus or signal sequences
    ProteinSequence.seq{i} = reduced_seq;
end

save('ProteinSequence.mat','ProteinSequence');
clear;
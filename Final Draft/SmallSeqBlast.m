function SmallSeqBlast(seqmat, protname, NumResults)
% Get and save BLAST search results for a specified Segment, marking after
% finished
    % TAKES ABOUT 3-4 HOURS TO RUN IF YOU GET 2000 HITS OF EACH, DATA SAVED
    % IN FOLDER
for kk = 1:numel(seqmat) % for each sequence in input matrix
    % Parameter values optimized for short sequences
    [requestID, requestTime] = blastncbi(seqmat(1,kk),'blastp', 'Expect',200000,...
        'Word', 2, 'Matrix', 'PAM30', 'MaxNumberSequences', NumResults, 'Database', 'pdb');
    blast_data = getblast(requestID, 'WaitTime', requestTime,'ToFile',"blast_" + protname + "_" +string(kk) + ".xml"); %Get and save BLAST
end
Donefile = fopen('Donefile.txt', 'w'); % write a file confirming that blast has been done before
fwrite(Donefile, '1');
fclose(Donefile);
end

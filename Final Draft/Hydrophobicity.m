function [Score, avgScore] = Hydrophobicity(Seq)
%Calculate Full and Average Hydropathy of input sequence
AminoAcids = struct(...     % Kyte-Doolittle Hydropathy values for each residue
    'Code', {'A','R','N','D','C','Q','E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', ...
    'S', 'T', 'W', 'Y', 'V'}, ...
    'Hydrophobicity', {1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,3.8,-3.9, ...
    1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2});
Score = 0;
Seq = char(Seq); % make sure input is a char if not already
for ii = 1:length(Seq) 
    if ismember(Seq(ii), [AminoAcids.Code]) == true
        ind = [AminoAcids.Code] == Seq(ii); %Get Score 
        Score = Score + AminoAcids(ind).Hydrophobicity; % Add Score to total Score
    else
        continue % Ignore weird other Amino Acids
    end
    avgScore = Score/length(Seq); % Average Score
end
function [percB, percH, percX] = PercentStructure(SecStructure)
% Calculate the Percent Sheet, Helix, and Other content, from array of sec.
% structure
lenTot = length(SecStructure);
numB = 0;
numH = 0;
numX = 0;
for ii = 1:lenTot % Iterate through array
    SS = SecStructure(ii); % Count numbers of each structure
    if SS == 'B'
        numB = numB + 1;
    elseif SS == 'H'
        numH = numH + 1;
    elseif SS == 'X'
        numX = numX + 1;
    end
    percB = numB/lenTot; % return percents
    percH = numH/lenTot;
    percX = numX/lenTot;
end
end
function secStruct = GetSecStruct(StructInstance,PDB)
% Get all ALigned residue indices
try % Make Sure Protein Has Helices
    x = [PDB.Helix.initChainID] == StructInstance.Chain; % Get all chain Helix indices
    PDBHelixes = transpose([PDB.Helix(x).initSeqNum; PDB.Helix(x).endSeqNum]); 
catch
   PDBHelixes = [0,0]; 
end
try % Make Sure Protein Has Sheets
    y = [PDB.Sheet.initChainID] == StructInstance.Chain; % Same formatting for helix and sheets
    PDBSheets = transpose([PDB.Sheet(y).initSeqNum; PDB.Sheet(y).endSeqNum]);
catch
    PDBSheets = [0,0];
end

secStruct = '';
iterator = StructInstance.AlignInd(1):StructInstance.AlignInd(2); % iterate over Aligned indices
for i = iterator
% Record whether the residue in said index position is B-sheet, Helix, or
% other
minIndSheet = PDBSheets(:,1) <= i; % Check within Sheet List
maxIndSheet = PDBSheets(:,2) >= i;
minIndHelix = PDBHelixes(:,1) <= i; % Check within Helix List
maxIndHelix = PDBHelixes(:,2) >= i;
inSheets = minIndSheet == maxIndSheet; % Checks to makes sure Aligned index is within list indices
inHelixes = minIndHelix == maxIndHelix;
if any(inSheets) == 1 % Adds B sheet to output
    secStruct = [secStruct, 'B'];    
elseif any(inHelixes(:)) == 1 % Adds helix to output
    secStruct = [secStruct,'H'];
else % Adds other structure to output
    secStruct = [secStruct,'X'];
end
end
end
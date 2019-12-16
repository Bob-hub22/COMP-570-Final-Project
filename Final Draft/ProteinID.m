function protStruct = ProteinID(filename,StructSize, sizeType, SeqNum)
%Construct Full Structure based on BLAST search Hits
%   filename: xml file BLAST search is saved to; Structsize: Number of Hits
%   you want to add to the structure; sizeType: Big or Small sequence?
%   SeqNum: Which of the three sequences are you using?
    % TAKES 24+ HOURS TO RUN IF YOU ARE GETTING 1000 RESULTS FROM BLAST,
    % DATA SAVED IN FOLDER
blast_struct = blastread(filename);
protStruct = struct('pdbID', {}, 'Chain', {}, 'FullSeq', {}, ...
    'AlignSeq', {}, 'AlignInd', {}, 'AlignScore', {}, 'SegDistOval', {}, ...
    'SegRMSDOval', {}, 'SegDistBlac', {}, 'SegRMSDBlac', {},'SecStruct', {} ,'AlignCoord', {}, 'AvgHydropathy', {}, 'FullHydropathy', {});

correctpdbBlac = getpdb("1BEB"); % Get pdb information for native B-Lactoglobulin and Ovalbumin as controls
correctpdbOval = getpdb("1OVA");
if sizeType == "Small"   % Figure out which B-Lac/Oval Indices to use
    if SeqNum == 1
        Seg1Blac = '20-26:A';
        Seg1Oval = '241-247:A';
    elseif SeqNum == 2
        Seg1Blac = '104-108:A';
        Seg1Oval = '366-370:A';
    elseif SeqNum == 3
        Seg1Blac = '118-121:A';
        Seg1Oval = '379-382:A';
    end
    
elseif sizeType == "Big"
    if SeqNum == 1
        Seg1Blac = '20-26:A';
        Seg1Oval = '240-247:A';
    elseif SeqNum == 2
        Seg1Blac = '102-108:A';
        Seg1Oval = '365-371:A';
    elseif SeqNum == 3
        Seg1Blac = '118-123:A';
        Seg1Oval = '377-383:A';
    end    
end

    blastsize = length([blast_struct.Hits]);  % Iterate through either the entire BLAST result list or the specified size, whichever is smaller
    if blastsize <= StructSize
        iterator = blastsize;
    else
        iterator = StructSize;
    end    
    
for i = 1: iterator
    protStruct(i).pdbID = blast_struct.Hits(i).Accession(1:4); % First 4 characters in BLAST Accession % are the PDB ID
    if blast_struct.Hits(i).Accession(end -1) == blast_struct.Hits(i).Accession(end) % Get the Chain from rest of accession #
        protStruct(i).Chain = lower(blast_struct.Hits(i).Accession(end)); % Some Chain Names are lowercase, shown by two letters in the accession #
    else
        protStruct(i).Chain = blast_struct.Hits(i).Accession(end);
    end
    
    protStruct(i).AlignSeq = blast_struct.Hits(i).Hsps(1).Alignment(3,:); % Get Alignment Sequence
    protStruct(i).AlignInd = blast_struct.Hits(i).Hsps(1).SubjectIndices(:); % Get Alignment Indices
    protStruct(i).AlignScore = blast_struct.Hits(i).Hsps(1).Score; % Get Alignment Score
    [protStruct(i).FullHydropathy, protStruct(i).AvgHydropathy] = Hydrophobicity(protStruct(i).AlignSeq); % Calculate Full and Average Hydropathies
    
    try % Pretty much any field can be missing from a PDB struct with no real pattern to it, so multiple try/catch mechanisms are required
        pdbstruct = getpdb(protStruct(i).pdbID); % Get PDB information
        protStruct(i).AlignCoord = GetMiddleCoordinates(protStruct(i).AlignInd, pdbstruct); % Get Centered Coordinates for Aligned Sequence
        protStruct(i).SecStruct = GetSecStruct(protStruct(i),pdbstruct); % Get Secondary Structure of Aligned Sequence
        
        X = [pdbstruct.Sequence.ChainID] == protStruct(i).Chain; % Get Chain Index
        protStruct(i).FullSeq = pdbstruct.Sequence(X).Sequence; % Get Full Sequence of Chain
        Seg2 = char(string(protStruct(i).AlignInd(1))+"-"... 
            +string(protStruct(i).AlignInd(2))+":"+string(protStruct(i).Chain));  % Put Indices and Chain in correct format for pdbsuperpose
    try
        [protStruct(i).SegDistOval, protStruct(i).SegRMSDOval] = ... 
            pdbsuperpose(correctpdbOval, pdbstruct, 'Display', false, 'Segment',... % Get superposition distances between Oval and Aligned Segments
            {Seg1Oval, Seg2});
        
    catch % Calculate second distance even if first one fails
        [protStruct(i).SegDistBlac, protStruct(i).SegRMSDBlac] = ...
        pdbsuperpose(correctpdbBlac, pdbstruct, 'Display', false, ... % Get superposition distances between Blac and Aligned Segments
        'Segment', {Seg1Blac, Seg2});
        continue
    end

    try
        [protStruct(i).SegDistBlac, protStruct(i).SegRMSDBlac] = ...
            pdbsuperpose(correctpdbBlac, pdbstruct, 'Display', false, ... % Get superposition distances between Blac and Aligned Segments
            'Segment', {Seg1Blac, Seg2});
    catch
        continue
    end

    
    catch
        continue
    end
end
end
function [oneTwoMatch, oneThreeMatch, twoThreeMatch, oneTwoThreeMatch] = BuildMatchStructs(Struct1, Struct2, Struct3)
onelist = getpdbs(Struct1); % Getting List of every PDBID from Structs for cross referencing
twolist = getpdbs(Struct2);
threelist = getpdbs(Struct3);

%Construct Structurea
oneTwoMatch = struct('pdbID',{},'AlignSeq1',{}, 'AlignSeq2', {}, 'AvgHydropathy', {},'SDHydropathy', {}, 'AvgPercB', {}, 'SDPercB',{}, 'SegDist', {},'AvgDistanceB', {}, 'AvgDistanceO', {});
oneThreeMatch = struct('pdbID',{},'AlignSeq1',{}, 'AlignSeq3', {}, 'AvgHydropathy', {},'SDHydropathy', {}, 'AvgPercB', {}, 'SDPercB',{}, 'SegDist', {},'AvgDistanceB', {}, 'AvgDistanceO', {});
twoThreeMatch = struct('pdbID',{},'AlignSeq2',{}, 'AlignSeq3', {}, 'AvgHydropathy', {},'SDHydropathy', {}, 'AvgPercB', {}, 'SDPercB',{},'SegDist', {},'AvgDistanceB', {}, 'AvgDistanceO', {});
oneTwoThreeMatch = struct('pdbID',{},'AlignSeq1',{}, 'AlignSeq2', {}, 'AlignSeq3', {}, 'AvgHydropathy', {},'SDHydropathy', {}, 'AvgPercB', {}, 'SDPercB',{},'SegDist12', {}, 'SegDist13', {}, 'SegDist23', {},'AvgDistanceB', {}, 'AvgDistanceO', {});

%While Not the Most Aesthetic way to fill in the structures, it is one of
%the more efficient ways to do it %%%%% May condense partially into a
%function
counter = [1,1,1,1];
for ii = 1:length(onelist)
    [Present, Loc12] = ismember(onelist(ii), twolist);
    if Present == 1 && isempty(Struct1(ii).SecStruct) == 0 && isempty(Struct2(Loc12).SecStruct) == 0 ... % Check All the necessary Conditions
            && isempty(Struct1(ii).AvgHydropathy) == 0 && isempty(Struct2(Loc12).AvgHydropathy) == 0 ...
            && isnan(Struct1(ii).AlignCoord(1)) == 0 && isnan(Struct2(Loc12).AlignCoord(1)) == 0
            
            % Fill out Structure
            [oneTwoMatch(counter(1)).pdbID, oneTwoMatch(counter(1)).AlignSeq1,oneTwoMatch(counter(1)).AlignSeq2, ...
                oneTwoMatch(counter(1)).AvgHydropathy,oneTwoMatch(counter(1)).SDHydropathy,...
                oneTwoMatch(counter(1)).AvgPercB,oneTwoMatch(counter(1)).SDPercB,~,Coords2,...
                ~,percB2,oneTwoMatch(counter(1)).SegDist,oneTwoMatch(counter(1)).AvgDistanceB, ...
                oneTwoMatch(counter(1)).AvgDistanceO] = FillStruct(Struct1,Struct2, ii, Loc12);
            counter(1) = counter(1) + 1; % Add to counter 

    else
        continue
    end
    [Present,Loc13] = ismember(onelist(ii), threelist);
    if Present == 1 && isempty(Struct1(ii).SecStruct) == 0 && isempty(Struct3(Loc13).SecStruct) == 0 ... % Check Necessary Conditions
            && isempty(Struct1(ii).AvgHydropathy) == 0 && isempty(Struct3(Loc13).AvgHydropathy) == 0 ...
            && isnan(Struct1(ii).AlignCoord(1)) == 0 && isnan(Struct3(Loc13).AlignCoord(1)) == 0
            
            % Fill Out Structure
            [oneThreeMatch(counter(2)).pdbID, oneThreeMatch(counter(2)).AlignSeq1,oneThreeMatch(counter(2)).AlignSeq3, ...
                oneThreeMatch(counter(2)).AvgHydropathy,oneThreeMatch(counter(2)).SDHydropathy,...
                oneThreeMatch(counter(2)).AvgPercB,oneThreeMatch(counter(2)).SDPercB,Coords1,Coords3,...
                percB1,percB3,oneThreeMatch(counter(2)).SegDist,oneThreeMatch(counter(2)).AvgDistanceB, ...
                oneThreeMatch(counter(2)).AvgDistanceO] = FillStruct(Struct1,Struct3, ii, Loc13);
            counter(2) = counter(2) + 1; % Add to Counter 
    else
        continue
    end
    
    if ismember(onelist(ii), twolist) == 1 % This one has a slightly different form from the others, so I didn't use my other function to fill it 
        [Present,Loc123] = ismember(onelist(ii), threelist);
        if Present == 1 && isempty(Struct1(ii).SecStruct) == 0 && isempty(Struct2(Loc12).SecStruct) == 0 ... % Check Necessary Conditions
            && isempty(Struct1(ii).AvgHydropathy) == 0 && isempty(Struct2(Loc12).AvgHydropathy) == 0 ...
            && isnan(Struct1(ii).AlignCoord(1)) == 0 && isnan(Struct2(Loc12).AlignCoord(1)) == 0 ...
            && isempty(Struct3(Loc123).AvgHydropathy) == 0 && isnan(Struct3(Loc123).AlignCoord(1)) == 0
            
            oneTwoThreeMatch(counter(4)).pdbID = Struct1(ii).pdbID;
            oneTwoThreeMatch(counter(4)).AlignSeq1 = string(Struct1(ii).AlignSeq);
            oneTwoThreeMatch(counter(4)).AlignSeq2 = string(Struct2(Loc12).AlignSeq);
            oneTwoThreeMatch(counter(4)).AlignSeq3 = string(Struct3(Loc123).AlignSeq);
                
            oneTwoThreeMatch(counter(4)).AvgHydropathy = mean([Struct1(ii).AvgHydropathy, Struct2(Loc12).AvgHydropathy, Struct3(Loc123).AvgHydropathy]);
            oneTwoThreeMatch(counter(4)).SDHydropathy = std([Struct1(ii).AvgHydropathy, Struct2(Loc12).AvgHydropathy,  Struct3(Loc123).AvgHydropathy]);
            
            oneTwoThreeMatch(counter(4)).AvgPercB = mean([percB1, percB2, percB3]);
            oneTwoThreeMatch(counter(4)).SDPercB = std([percB1, percB2, percB3]);
            
            oneTwoThreeMatch(counter(4)).SegDist12 = RMSDistance(Coords1, Coords2);
            oneTwoThreeMatch(counter(4)).SegDist13 = RMSDistance(Coords1, Coords3);
            oneTwoThreeMatch(counter(4)).SegDist23 = RMSDistance(Coords2, Coords3);
            oneTwoThreeMatch(counter(4)).AvgDistanceB = mean([Struct1(ii).SegDistBlac, Struct2(Loc12).SegDistBlac, Struct3(Loc123).SegDistBlac]);
            oneTwoThreeMatch(counter(4)).AvgDistanceO = mean([Struct1(ii).SegDistOval, Struct2(Loc12).SegDistOval, Struct3(Loc123).SegDistOval]);
            counter(4) = counter(4) + 1; % Add to Counter 
        else
            continue
        end
    end
end
for jj = 1:length(twolist)
    [Present, Loc23] = ismember(twolist(jj), threelist);
        if Present == 1 && isempty(Struct2(jj).SecStruct) == 0 && isempty(Struct3(Loc23).SecStruct) == 0 ... % Check Necessary Conditions
            && isempty(Struct2(jj).AvgHydropathy) == 0 && isempty(Struct3(Loc23).AvgHydropathy) == 0 ...
            && isnan(Struct2(jj).AlignCoord(1)) == 0 && isnan(Struct3(Loc23).AlignCoord(1)) == 0
            
            % Fill out Structure
            [twoThreeMatch(counter(3)).pdbID, twoThreeMatch(counter(3)).AlignSeq2,twoThreeMatch(counter(3)).AlignSeq3, ...
                twoThreeMatch(counter(3)).AvgHydropathy,twoThreeMatch(counter(3)).SDHydropathy,...
                twoThreeMatch(counter(3)).AvgPercB,twoThreeMatch(counter(3)).SDPercB,~,~,...
                ~,~,twoThreeMatch(counter(3)).SegDist,twoThreeMatch(counter(3)).AvgDistanceB, ...
                twoThreeMatch(counter(3)).AvgDistanceO] = FillStruct(Struct2,Struct3, jj, Loc23);
            counter(3) = counter(3) + 1; % Add to Counter
        else
            continue
        end
end
end

function pdblist = getpdbs(Struct) 
%Get each pdbID in string form
pdblist = strings(1,length([Struct.Chain]));
for kk = 1:length(pdblist)
    pdblist(1,kk) = string(Struct(kk).pdbID);
end
end

function [ID, AlignSeq1,AlignSeq2,AvgHydro,StDHydro,AvgPercB,SDPercB,Coords1,Coords2,percB1,percB2,SegDist, AvgDistanceB, AvgDistanceO] = FillStruct(Struct1,Struct2, ii, Loc12)
%Fill out Structure with all Necessary Components
ID = Struct1(ii).pdbID; % Get ID 
AlignSeq1 = string(Struct1(ii).AlignSeq); % Get First Align Sequence
AlignSeq2= string(Struct2(Loc12).AlignSeq); % Get Second Align Sequence
AvgHydro = mean([Struct1(ii).AvgHydropathy, Struct2(Loc12).AvgHydropathy]); % Calc Avg Hydropathy of Pair
StDHydro = std([Struct1(ii).AvgHydropathy, Struct2(Loc12).AvgHydropathy]); % Calc StD of Hydropathy of Pair
          
[percB1,~,~] = PercentStructure(Struct1(ii).SecStruct);
[percB2,~,~] = PercentStructure(Struct2(Loc12).SecStruct);
AvgPercB = mean([percB1, percB2]); % Calc Mean % B-Sheet of pair
SDPercB = std([percB1, percB2]); % Calc StD of % B-Sheet in Pair
Coords1 = Struct1(ii).AlignCoord; % Get First Coordinates
Coords2 = Struct2(Loc12).AlignCoord; %Get Second Coordinates
SegDist = RMSDistance(Coords1, Coords2); % Calculate RMSD
AvgDistanceB = mean([Struct1(ii).SegDistBlac, Struct2(Loc12).SegDistBlac]);
AvgDistanceO = mean([Struct1(ii).SegDistOval, Struct2(Loc12).SegDistOval]);
end

function RMS = RMSDistance(Coords1, Coords2) 
%Use RMSD formula
xcomp = (Coords2(1)-Coords1(1))*(Coords2(1)-Coords1(1));
ycomp = (Coords2(2)-Coords1(2))*(Coords2(2)-Coords1(2));
zcomp = (Coords2(3)-Coords1(3))*(Coords2(3)-Coords1(3));
RMS = sqrt(xcomp + ycomp + zcomp);
end
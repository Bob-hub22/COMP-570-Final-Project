function SSBODist = GetDist_SS_Pairs(Struct)
% Return Matrix of paired distance & % Secondary Structure Measurements for
% each hit (Skipping empty fields)
lenStruct = length(Struct);    
SSBODist = [];
    for ii = 1:lenStruct
        %Skip Empty Structs
        if isempty(Struct(ii).SegDistOval) == 0 && isempty(Struct(ii).SegRMSDOval) == 0 && isempty(Struct(ii).SegDistBlac) == 0 ...
                && isempty(Struct(ii).SegRMSDBlac) == 0 && isempty(Struct(ii).SecStruct) == 0 
            
            [percB, percH, percX] = PercentStructure(Struct(ii).SecStruct); % Get % Secondary Structures
            % Get Distance Measurements
            SSBODist = [SSBODist;[Struct(ii).SegDistOval, Struct(ii).SegRMSDOval, Struct(ii).SegDistBlac, Struct(ii).SegRMSDBlac, percB, percH, percX]];
        else
            continue
        end
    end
end
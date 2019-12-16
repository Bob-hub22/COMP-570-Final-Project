function coords = GetMiddleCoordinates(SeqInds,PDB)
% Get Centroid Coordinates for Segments
m = mean(SeqInds);
coords = zeros(1,3);
try %Return NaNs if Coordinates cannot be found
    if floor(m) == m  %If # residues is odd
        coords = CoordinateGet(m, PDB); % return coordinates for middle residue
    else % if # residues is even
        mlow = floor(m); % Get coordinates for the 2 center residues
        mhigh = ceil(m);
        coordlow = CoordinateGet(mlow, PDB);
        coordhigh = CoordinateGet(mhigh, PDB);
        %Return Average Coordinates
        coords = [(coordlow(1)+coordhigh(1))/2, (coordlow(2)+coordhigh(2))/2, (coordlow(3)+coordhigh(3))/2]; 
    end
catch
    coords = [NaN, NaN, NaN];
end
end

function coordinates = CoordinateGet(meanInd,PDB)
%Get Coordinates of a specified index
res = [PDB.Model.Atom.resSeq] == meanInd;
residueXYZ = transpose([PDB.Model.Atom(res).X; PDB.Model.Atom(res).Y; PDB.Model.Atom(res).Z]);
coordinates = residueXYZ(2,:);
end
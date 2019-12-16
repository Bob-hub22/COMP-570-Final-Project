function [onetwoSeq, onethreeSeq, twothreeSeq, onetwothreeSeq, onetwo, onethree, twothree, onetwothree] = OneTwoThreePair(Struct1, Struct2, Struct3)
% return full lists of proteins that contain at least 2 of the segments
% (Not necessary in final draft, but the BuildMatchStructs Function was
% built using the algorithm from this original function);
onelist = getpdbs(Struct1);
twolist = getpdbs(Struct2);
threelist = getpdbs(Struct3);
onetwo = false(1,length(onelist));
onethree = false(1,length(onelist));
twothree = false(1,length(twolist));
onetwothree = false(1,length(onelist));

for ii = 1:length(onelist)
    [onetwo(1,ii), Loc12] = ismember(onelist(ii), twolist);
    [onethree(1,ii),Loc13] = ismember(onelist(ii), threelist);
    
    if ismember(onelist(ii), twolist) == 1
    [onetwothree(1,ii),Loc123] = ismember(onelist(ii), threelist);
    end
end
for jj = 1:length(twolist)
    [twothree(1,jj), Loc23] = ismember(twolist(jj), threelist);
end
    onetwoSeq = onelist(onetwo);
    onethreeSeq = onelist(onethree);
    twothreeSeq = twolist(twothree);
    onetwothreeSeq = onelist(onetwothree);
end

function pdblist = getpdbs(Struct)
pdblist = strings(1,length([Struct.Chain]));
for kk = 1:length(pdblist)
    pdblist(1,kk) = string(Struct(kk).pdbID);
end
end
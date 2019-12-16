function [PercBhydro, PercHhydro, PercXhydro, AvgHydros, FullHydros, seqInd] = SSHydroCompare(Struct)
% Building Lists for the % Sec.Structure concent and Hydropathy
% measurements for direct comparison
lenTot = length(Struct);
uniqueMat = strings;
seqInd = {};
AvgHydros = [];
FullHydros = [];
AvgPercB = {};
AvgPercH = {};
AvgPercX = {};
for ii = 1:lenTot
   if isempty(Struct(ii).SecStruct) == 0 && isempty(Struct(ii).AvgHydropathy) == 0 && isempty(Struct(ii).FullHydropathy) == 0       
       [percB, percH, percX] = PercentStructure(Struct(ii).SecStruct);
       strseq = string(Struct(ii).AlignSeq);
       if  ismember(strseq, uniqueMat) == 0   
           uniqueMat = [uniqueMat, strseq];
           seqInd{1,end+1} = {ii};          
           AvgPercB{end+1} = {percB};
           AvgPercH{end+1} = {percH};
           AvgPercX{end+1} = {percX};
           AvgHydros = [AvgHydros, Struct(ii).AvgHydropathy];
           FullHydros = [FullHydros, Struct(ii).FullHydropathy];
       elseif ismember(strseq, uniqueMat) == 1
           x = find(uniqueMat == strseq)-1;
           seqInd{x} = [seqInd{x}, ii];
           AvgPercB{x} = [AvgPercB{x}, percB];
           AvgPercH{x} = [AvgPercH{x}, percH];
           AvgPercX{x} = [AvgPercX{x}, percX];
       end     
   else
       continue
   end
   

end
   PercBhydro = zeros(length(AvgPercB),1);
   PercHhydro = zeros(length(AvgPercH),1);
   PercXhydro = zeros(length(AvgPercX),1);
   
   for jj = 1:length(FullHydros)
       PercBhydro(jj,1) = mean(cell2mat(AvgPercB{jj}));
       PercHhydro(jj,1) = mean(cell2mat(AvgPercH{jj}));
       PercXhydro(jj,1) = mean(cell2mat(AvgPercX{jj}));
   end
end

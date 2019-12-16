function [mergedBlac, mergedOval, BlacSmall1, BlacSmall2, BlacSmall3, BlacBig2, BlacBig3, OvalSmall1, OvalSmall2, OvalSmall3,OvalBig1, OvalBig2, OvalBig3] = InitiallizeData(NumResults)
% Initiallize Data, building BLAST results and Protein Data if not
% previously collected
    % TAKES VERY LONG TO RUN THE FIRST TIME, DATA SAVED IN FOLDER
seqBlac = ["YSLAMAA", "LFCME", "VCQC"]; %% important matching sequences from experimental results
seqOval = ["SMLVLLP", "LFCIK", "FFGR"];

seqBlacfull= ["YSLAMAA", "YLLFCME", "VCQCLV"]; %% Larger encapulating sequences of above, that show the same Sec. Structure in vivo
seqOvalfull= ["MSMLVLLP","FLFCIKH", "VLFFGRC"];
mergedBlac = [seqBlac,seqBlacfull(2:3)];   % Merge Short and long Sequences into the same matrices
mergedOval = [seqOval,seqOvalfull];           %Note that The first Blac Seq has the same short and long sequences


Donefile = fopen('Donefile.txt', 'r'); % See if blast was performed previously (BlASTS of all results take about 2-3 hours to run)
Doneb4 = fread(Donefile);
fclose(Donefile);
if Doneb4 == 48 % If not done before
    SmallSeqBlast(mergedBlac, 'Blac', NumResults);  %Get Blast results of all sequences in matrix, save to separate xml files
    SmallSeqBlast(mergedOval, 'Oval', NumResults);
    disp("BLAST searches Obtained; Obtaining Protein Info...");
else
    disp('Blast searches already performed, now loading Protein info...'); % Cute little message to show that blast searches have already been done
end
[BlacSmall1, BlacSmall2, BlacSmall3, BlacBig2, BlacBig3, OvalSmall1, OvalSmall2, OvalSmall3,OvalBig1, OvalBig2, OvalBig3] = InitiallizePDB(NumResults);
end


function [ BlacSmall1, BlacSmall2, BlacSmall3, BlacBig2, BlacBig3, OvalSmall1, OvalSmall2, OvalSmall3,OvalBig1, OvalBig2, OvalBig3] = InitiallizePDB(StructSize)
DoneP = fopen('ProtInfoDone.txt'); % See if Structures have already been constructed and saved into files
ProtDoneb4 = fread(DoneP);
fclose(DoneP);
if ProtDoneb4 == 48
     BlacSmall1 = ProteinID('blast_BlacSmall_1.xml', StructSize, "Small", 1);  %Create and Save Structures for all sequences
     save('BlacSmall1.mat', 'BlacSmall1');

    BlacSmall2 = ProteinID('blast_BlacSmall_2.xml',StructSize, "Small", 2);
    save('BlacSmall2.mat', 'BlacSmall2');

    BlacSmall3 = ProteinID('blast_BlacSmall_3.xml',StructSize, "Small", 3);
    save('BlacSmall3.mat', 'BlacSmall3');

    BlacBig2 = ProteinID('blast_BlacBig_2.xml',StructSize, "Big", 2);
    save('BlacBig2.mat', 'BlacBig2');

    BlacBig3 = ProteinID('blast_BlacBig_3.xml',StructSize, "Big", 3);
    save('BlacBig3.mat', 'BlacBig3');

    OvalSmall1 = ProteinID('blast_OvalSmall_1.xml',StructSize, "Small", 1);
    save('OvalSmall1.mat', 'OvalSmall1');

    OvalSmall2 = ProteinID('blast_OvalSmall_2.xml',StructSize, "Small", 2);
    save('OvalSmall2.mat', 'OvalSmall2');

    OvalSmall3 = ProteinID('blast_OvalSmall_3.xml',StructSize, "Small", 3);
    save('OvalSmall3.mat', 'OvalSmall3');

    OvalBig1 = ProteinID('blast_OvalBig_1.xml',StructSize, "Big", 1);
    save('OvalBig1.mat', 'OvalBig1');

    OvalBig2 = ProteinID('blast_OvalBig_2.xml',StructSize, "Big", 2);
    save('OvalBig2.mat', 'OvalBig2');

    OvalBig3 = ProteinID('blast_OvalBig_3.xml',StructSize, "Big", 3);
    save('OvalBig3.mat', 'OvalBig3');

    disp('Protein Info Loaded');
    fileID = fopen('ProtInfoDone.txt', 'w');  % Mark that structures have been constructed
    fwrite(fileID, '1');
    fclose(fileID);
else
    load('BlacSmall1.mat'); load('BlacSmall2.mat'); load('BlacSmall3.mat');  % Load all structuress
    load('BlacBig2.mat'); load('BlacBig3.mat');

    load('OvalSmall1.mat'); load('OvalSmall2.mat'); load('OvalSmall3.mat');
    load('OvalBig1.mat'); load('OvalBig2.mat'); load('OvalBig3.mat');
    disp('Protein Info Loaded');
end
end
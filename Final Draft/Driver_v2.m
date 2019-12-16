%% Get/Load Protein Data
[mergedBlac, mergedOval, BlacSmall1, BlacSmall2, BlacSmall3, BlacBig2, BlacBig3, OvalSmall1, OvalSmall2, OvalSmall3,OvalBig1, OvalBig2, OvalBig3] =InitiallizeData(1000);
%% Driver for Comparing Distance measurements to Secondary structure percent measurements
SSDistBLS1 = GetDist_SS_Pairs(BlacSmall1); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistBLS2 = GetDist_SS_Pairs(BlacSmall2); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistBLS3 = GetDist_SS_Pairs(BlacSmall3); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistBLB2 = GetDist_SS_Pairs(BlacBig2); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistBLB3 = GetDist_SS_Pairs(BlacBig3); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistOvS1 = GetDist_SS_Pairs(OvalSmall1); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistOvS2 = GetDist_SS_Pairs(OvalSmall2); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistOvS3 = GetDist_SS_Pairs(OvalSmall3); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistOvB1 = GetDist_SS_Pairs(OvalBig1); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistOvB2 = GetDist_SS_Pairs(OvalBig2); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDistOvB3 = GetDist_SS_Pairs(OvalBig3); %OvalD, OvalR, BlacD, BlacR, percB, percH, percX
SSDist = [SSDistBLS1; SSDistBLS2; SSDistBLS3;SSDistBLB2 ;SSDistBLB3;SSDistOvS1;SSDistOvS2;SSDistOvS3;SSDistOvB1;SSDistOvB2;SSDistOvB3];
%% Plot Above
scatter(SSDist(:,1),SSDist(:,2));
subplot(4,6,2);
scatter(SSDist(:,1),SSDist(:,3));
subplot(4,6,3);
scatter(SSDist(:,1),SSDist(:,4));
subplot(4,6,4);
scatter(SSDist(:,1),SSDist(:,5));
subplot(4,6,5);
scatter(SSDist(:,1),SSDist(:,6));
subplot(4,6,6);
scatter(SSDist(:,1),SSDist(:,7));

subplot(4,6,8);
scatter(SSDist(:,2),SSDist(:,3));
subplot(4,6,9);
scatter(SSDist(:,2),SSDist(:,4));
subplot(4,6,10);
scatter(SSDist(:,2),SSDist(:,5));
subplot(4,6,11);
scatter(SSDist(:,2),SSDist(:,6));
subplot(4,6,12);
scatter(SSDist(:,2),SSDist(:,7));

subplot(4,6,15);
scatter(SSDist(:,3),SSDist(:,4));
subplot(4,6,16);
scatter(SSDist(:,3),SSDist(:,5));
subplot(4,6,17);
scatter(SSDist(:,3),SSDist(:,6));
subplot(4,6,18);
scatter(SSDist(:,3),SSDist(:,7));

subplot(4,6,22);
scatter(SSDist(:,4),SSDist(:,5));
subplot(4,6,23);
scatter(SSDist(:,4),SSDist(:,6));
subplot(4,6,24);
scatter(SSDist(:,4),SSDist(:,7));
%% Important Ones from Above
subplot(2,3,1);
scatter(SSDist(:,5),SSDist(:,1));
[fitout, gof] = fit(SSDist(:,5),SSDist(:,1), 'poly1');
hold on 
plot(fitout);
axis([0,max(SSDist(:,5)),0,max(SSDist(:,1))]);
title("R^2 = " + string(gof.rsquare));
xlabel("% Sheet", 'FontWeight', 'bold');
ylabel("Superpose Distance from Oval", 'FontWeight', 'bold');
hold off

subplot(2,3,2);
scatter(SSDist(:,6),SSDist(:,1));
[fitout, gof] = fit(SSDist(:,6),SSDist(:,1), 'poly1');
hold on 
plot(fitout);
axis([0,max(SSDist(:,6)),0,max(SSDist(:,1))]);
title("R^2 = " + string(gof.rsquare));
xlabel("% Helix", 'FontWeight', 'bold');
ylabel("Superpose Distance Oval", 'FontWeight', 'bold');
hold off

subplot(2,3,3);
scatter(SSDist(:,7),SSDist(:,1));
[fitout, gof] = fit(SSDist(:,7),SSDist(:,1), 'poly1');
hold on 
plot(fitout);
axis([0,max(SSDist(:,7)),0,max(SSDist(:,1))]);
title("R^2 = " + string(gof.rsquare));
xlabel("% Other", 'FontWeight', 'bold');
ylabel("Superpose Distance Oval", 'FontWeight', 'bold');
hold off

subplot(2,3,4);
scatter(SSDist(:,5),SSDist(:,3));
[fitout, gof] = fit(SSDist(:,5),SSDist(:,3), 'poly1');
hold on 
plot(fitout);
axis([0,max(SSDist(:,5)),0,max(SSDist(:,3))]);
title("R^2 = " + string(gof.rsquare));
xlabel("% Sheet", 'FontWeight', 'bold');
ylabel("Superpose Distance from BLac", 'FontWeight', 'bold');
hold off

subplot(2,3,5);
scatter(SSDist(:,6),SSDist(:,3));
[fitout, gof] = fit(SSDist(:,6),SSDist(:,3), 'poly1');
hold on 
plot(fitout);
axis([0,max(SSDist(:,6)),0,max(SSDist(:,3))]);
title("R^2 = " + string(gof.rsquare));
xlabel("% Helix", 'FontWeight', 'bold');
ylabel("Superpose Distance from BLac", 'FontWeight', 'bold');
hold off

subplot(2,3,6);
scatter(SSDist(:,7),SSDist(:,3));
[fitout, gof] = fit(SSDist(:,7),SSDist(:,3), 'poly1');
hold on 
plot(fitout);
axis([0,max(SSDist(:,7)),0,max(SSDist(:,3))]);
title("R^2 = " + string(gof.rsquare));
xlabel("% Other", 'FontWeight', 'bold');
ylabel("Superpose Distance from BLac", 'FontWeight', 'bold');
hold off
sgtitle('Correlations Between Secondary Structure and pdbsuperpose Distances');
%% Driver for Comparing average SS frequency to hydrophobicity
[PercBBS1,PercHBS1,PercXBS1, AvgHydroBS1, FullHydroBS1, seqIndBS1] = SSHydroCompare(BlacSmall1);
[PercBBS2,PercHBS2,PercXBS2, AvgHydroBS2, FullHydroBS2, seqIndBS2] = SSHydroCompare(BlacSmall2);
[PercBBS3,PercHBS3,PercXBS3, AvgHydroBS3, FullHydroBS3, seqIndBS3] = SSHydroCompare(BlacSmall3);
[PercBBB2,PercHBB2,PercXBB2, AvgHydroBB2, FullHydroBB2, seqIndBB2] = SSHydroCompare(BlacBig2);
[PercBBB3,PercHBB3,PercXBB3, AvgHydroBB3, FullHydroBB3, seqIndBB3] = SSHydroCompare(BlacBig3);
[PercBOS1,PercHOS1,PercXOS1, AvgHydroOS1, FullHydroOS1, seqIndOS1] = SSHydroCompare(OvalSmall1);
[PercBOS2,PercHOS2,PercXOS2, AvgHydroOS2, FullHydroOS2, seqIndOS2] = SSHydroCompare(OvalSmall2);
[PercBOS3,PercHOS3,PercXOS3, AvgHydroOS3, FullHydroOS3, seqIndOS3] = SSHydroCompare(OvalSmall3);
[PercBOB1,PercHOB1,PercXOB1, AvgHydroOB1, FullHydroOB1, seqIndOB1] = SSHydroCompare(OvalBig1);
[PercBOB2,PercHOB2,PercXOB2, AvgHydroOB2, FullHydroOB2, seqIndOB2] = SSHydroCompare(OvalBig2);
[PercBOB3,PercHOB3,PercXOB3, AvgHydroOB3, FullHydroOB3, seqIndOB3] = SSHydroCompare(OvalBig3);
%% Plot Above
subplot(2,3,1);
scatter(AvgHydroBS1,PercBBS1,'filled');
hold on
scatter(AvgHydroBS1, PercBBS1, 'r');
xline(1.4036,'--','LineWidth', 0.75);
xline(1.6659,'LineWidth', 0.75);
xlabel('Average Hydropathy','FontWeight', 'bold')
ylabel('% Sheet Content','FontWeight', 'bold')
legend('Short Sequence', 'Long Sequence','BLac Small Seg. Hydro.','BLac Long Seg. Hydro.', 'Location', 'best');
title('B-Lactoglobulin, Segment 1');
hold off

subplot(2,3,2);
scatter(AvgHydroBS2,PercBBS2,'filled');
hold on 
scatter(AvgHydroBB2,PercBBB2,'r');
xline(1.4036,'--','LineWidth', 0.75);
xline(1.6659,'LineWidth', 0.75);
xlabel('Average Hydropathy','FontWeight', 'bold')
ylabel('% Sheet Content','FontWeight', 'bold')
legend('Short Sequence', 'Long Sequence','BLac Small Seg. Hydro.','BLac Long Seg. Hydro','Location', 'best');
title('B-Lactoglobulin, Segment 2');
hold off

subplot(2,3,3);
scatter(AvgHydroBS3,PercBBS3,'filled');
hold on
scatter(AvgHydroBB3,PercBBB3,'r');
xline(1.4036,'--','LineWidth', 0.75);
xline(1.6659,'LineWidth', 0.75);
xlabel('Average Hydropathy','FontWeight', 'bold')
ylabel('% Sheet Content','FontWeight', 'bold')
legend('Short Sequence', 'Long Sequence','BLac Small Seg. Hydro.','BLac Long Seg. Hydro.','Location', 'best');
title('B-Lactoglobulin, Segment 3');
hold off

subplot(2,3,4)
scatter(AvgHydroOS1,PercBOS1,'filled');
hold on
scatter(AvgHydroOB1,PercBOB1,'r');
xline(1.4240,'--','LineWidth', 0.75);
xline(1.6845,'LineWidth', 0.75);
xlabel('Average Hydropathy','FontWeight', 'bold')
ylabel('% Sheet Content','FontWeight', 'bold')
legend('Short Sequence', 'Long Sequence','Oval Small Seg. Hydro.','Oval Long Seg. Hydro.','Location', 'best');
title('Ovalbumin, Segment 1');
hold off

subplot(2,3,5)
scatter(AvgHydroOS2,PercBOS2,'filled');
hold on
scatter(AvgHydroOB2,PercBOB2,'r');
xline(1.4240,'--','LineWidth', 0.75);
xline(1.6845,'LineWidth', 0.75);
xlabel('Average Hydropathy','FontWeight', 'bold')
ylabel('% Sheet Content','FontWeight', 'bold')
legend('Short Sequence', 'Long Sequence','Oval Small Seg. Hydro.','Oval Long Seg. Hydro.', 'Location', 'best');
title('Ovalbumin, Segment 2');
hold off

subplot(2,3,6)
scatter(AvgHydroOS3,PercBOS3,'filled');
hold on
scatter(AvgHydroOB3,PercBOB3,'r');
xline(1.4240,'--','LineWidth', 0.75);
xline(1.6845,'LineWidth', 0.75);
xlabel('Average Hydropathy','FontWeight', 'bold')
ylabel('% Sheet Content','FontWeight', 'bold')
legend('Short Sequence', 'Long Sequence','Oval Small Seg. Hydro.','Oval Long Seg. Hydro.','Location', 'best');
title('Ovalbumin, Segment 3');
hold off
sgtitle('Average Hydrophobicity vs % Sheet Content for Independent Segments');
%% Same measurement as above and more, limited to the average values for proteins that contain >1 segment
% Get pairs as structs containing pdbID, Sequences, Average Hydropathies &
% % B-sheet (with associated StD. for error), and Distance between pair
% based on Crystal Structure (pair is used loosely, as it also includes the
% triplets of all three sequences appearing in the same protein)
[oneTwoMatchBS, oneThreeMatchBS, twoThreeMatchBS, oneTwoThreeMatchBS] = BuildMatchStructs(BlacSmall1, BlacSmall2, BlacSmall3);
[oneTwoMatchBB, oneThreeMatchBB, twoThreeMatchBB, oneTwoThreeMatchBB] = BuildMatchStructs(BlacSmall1, BlacBig2, BlacBig3);
[oneTwoMatchOS, oneThreeMatchOS, twoThreeMatchOS, oneTwoThreeMatchOS] = BuildMatchStructs(OvalSmall1, OvalSmall2, OvalSmall3);
[oneTwoMatchOB, oneThreeMatchOB, twoThreeMatchOB, oneTwoThreeMatchOB] = BuildMatchStructs(OvalBig1, OvalBig2, OvalBig3);
%% Comparing Distance Measurements, Clearly Inconsistent, no real correlation not useful at all Supports Data comparing Suprpose vs SS content
DistB = [[oneTwoMatchBS.AvgDistanceB],[oneThreeMatchBB.AvgDistanceB],[twoThreeMatchOS.AvgDistanceB],[oneTwoThreeMatchOB.AvgDistanceB]];
DistO = [[oneTwoMatchBS.AvgDistanceO],[oneThreeMatchBB.AvgDistanceO],[twoThreeMatchOS.AvgDistanceO],[oneTwoThreeMatchOB.AvgDistanceO]];
scatter(DistB,DistO);
%% 3D plots of Hydropathy, % BSheet, and Segment Distance for Proteins containing all three segments

subplot(2,2,1)
z = [[oneTwoThreeMatchBS.AvgHydropathy]; [oneTwoThreeMatchBS.AvgPercB]; [oneTwoThreeMatchBS.SegDist12]];
surf(z,'FaceAlpha',0.4, 'FaceColor', [0.9290 0.6940 0.1250]);
hold on 
z2 =[[oneTwoThreeMatchBS.AvgHydropathy]; [oneTwoThreeMatchBS.AvgPercB]; [oneTwoThreeMatchBS.SegDist13]];
surf(z2,'FaceAlpha',0.5, 'FaceColor', [0.8500 0.3250 0.0980]);
z3 =[[oneTwoThreeMatchBS.AvgHydropathy]; [oneTwoThreeMatchBS.AvgPercB]; [oneTwoThreeMatchBS.SegDist23]];
surf(z3,'FaceAlpha',0.5, 'FaceColor', [0 0.4470 0.7410]);
xlabel('Average Hydropathy');
ylabel('% B-Sheet Content');
zlabel('Distance Beween Segments');
legend('Seg 1 + Seg 2', 'Seg 1 + Seg 3', 'Seg 2 + Seg 3','Location','northeast')
hold off
subplot(2,2,2)
z = [[oneTwoThreeMatchBB.AvgHydropathy]; [oneTwoThreeMatchBB.AvgPercB]; [oneTwoThreeMatchBB.SegDist12]];
surf(z,'FaceAlpha',0.5,'FaceColor', [0.9290 0.6940 0.1250]);
hold on 
z2 =[[oneTwoThreeMatchBB.AvgHydropathy]; [oneTwoThreeMatchBB.AvgPercB]; [oneTwoThreeMatchBB.SegDist13]];
surf(z2,'FaceAlpha',0.5, 'FaceColor', [0.8500 0.3250 0.0980]);
z3 =[[oneTwoThreeMatchBB.AvgHydropathy]; [oneTwoThreeMatchBB.AvgPercB]; [oneTwoThreeMatchBB.SegDist23]];
surf(z3,'FaceAlpha',0.5, 'FaceColor', [0 0.4470 0.7410]);
xlabel('Average Hydropathy');
ylabel('% B-Sheet Content');
zlabel('Distance Beween Segments');
legend('Seg 1 + Seg 2', 'Seg 1 + Seg 3', 'Seg 2 + Seg 3','Location','northeast')
hold off
subplot(2,2,3)
z = [[oneTwoThreeMatchOS.AvgHydropathy]; [oneTwoThreeMatchOS.AvgPercB]; [oneTwoThreeMatchOS.SegDist12]];
surf(z,'FaceAlpha',0.5,'FaceColor', [0.9290 0.6940 0.1250]);
hold on 
z2 =[[oneTwoThreeMatchOS.AvgHydropathy]; [oneTwoThreeMatchOS.AvgPercB]; [oneTwoThreeMatchOS.SegDist13]];
surf(z2,'FaceAlpha',0.5, 'FaceColor', [0.8500 0.3250 0.0980]);
z3 =[[oneTwoThreeMatchOS.AvgHydropathy]; [oneTwoThreeMatchOS.AvgPercB]; [oneTwoThreeMatchOS.SegDist23]];
surf(z3,'FaceAlpha',0.5, 'FaceColor', [0 0.4470 0.7410]);
xlabel('Average Hydropathy','FontWeight', 'bold');
ylabel('% B-Sheet Content','FontWeight', 'bold');
zlabel('Distance Beween Segments','FontWeight', 'bold');
legend('Seg 1 + Seg 2', 'Seg 1 + Seg 3', 'Seg 2 + Seg 3','Location','northeast')
hold off
subplot(2,2,4)
z4 = [[oneTwoThreeMatchOB.AvgHydropathy]; [oneTwoThreeMatchOB.AvgPercB]; [oneTwoThreeMatchOB.SegDist12]];
surf(z4,'FaceAlpha',0.5,'FaceColor', [0.9290 0.6940 0.1250]);
hold on 
z5 =[[oneTwoThreeMatchOB.AvgHydropathy]; [oneTwoThreeMatchOB.AvgPercB]; [oneTwoThreeMatchOB.SegDist13]];
surf(z5,'FaceAlpha',0.5, 'FaceColor', [0.8500 0.3250 0.0980]);
z6 =[[oneTwoThreeMatchOB.AvgHydropathy]; [oneTwoThreeMatchOB.AvgPercB]; [oneTwoThreeMatchOB.SegDist23]];
surf(z6,'FaceAlpha',0.5, 'FaceColor', [0 0.4470 0.7410]);
xlabel('Average Hydropathy','FontWeight', 'bold');
ylabel('% B-Sheet Content','FontWeight', 'bold');
zlabel('Distance Beween Segments','FontWeight', 'bold');
legend('Seg 1 + Seg 2', 'Seg 1 + Seg 3', 'Seg 2 + Seg 3','Location','northeast')
hold off
%% Looking at Relationship between Hydropathy and Distance
subplot(1,2,1);
hydros = [[oneTwoThreeMatchBS.AvgHydropathy],[oneTwoThreeMatchBB.AvgHydropathy],[oneTwoThreeMatchOS.AvgHydropathy],[oneTwoThreeMatchOB.AvgHydropathy]];
dist12 = [[oneTwoThreeMatchBS.SegDist12],[oneTwoThreeMatchBB.SegDist12],[oneTwoThreeMatchOS.SegDist12],[oneTwoThreeMatchOB.SegDist12]];
dist13 = [[oneTwoThreeMatchBS.SegDist13],[oneTwoThreeMatchBB.SegDist13],[oneTwoThreeMatchOS.SegDist13],[oneTwoThreeMatchOB.SegDist13]];
dist23 = [[oneTwoThreeMatchBS.SegDist23],[oneTwoThreeMatchBB.SegDist23],[oneTwoThreeMatchOS.SegDist23],[oneTwoThreeMatchOB.SegDist23]];
[fitoutpoly, gof] = fit(transpose([hydros,hydros,hydros]), transpose([dist12,dist13,dist23]), 'poly2');
scatter(hydros,dist12);
hold on
c = [0.8500 0.3250 0.0980];
scatter(hydros,dist13, '^', 'MarkerEdgeColor', c);
c = [0.9290 0.6940 0.1250];
scatter(hydros,dist23,'+','MarkerEdgeColor', c);
xline(mean([1.4036, 1.4240]), '--', 'LineWidth', 0.75);
xline(mean([1.6845,1.6659]),'LineWidth', 0.75);
legend('Seg 1 + Seg 2', 'Seg 1 + Seg 3','Seg 2 + Seg 3', 'Avg Hydro, Short Seg', 'Avg Hydro, Long Seg', 'Location','northwest');
plot(fitoutpoly);
xlabel('Average Hydropathy','FontWeight', 'bold');
ylabel('Distance Between Segments','FontWeight', 'bold');
title('Hydropathy vs RMSD between segments, R^2 = ' + string(gof.rsquare));
hold off
subplot(1,2,2)
hydros = [[oneTwoThreeMatchBS.AvgHydropathy],[oneTwoThreeMatchBB.AvgHydropathy],[oneTwoThreeMatchOS.AvgHydropathy],[oneTwoThreeMatchOB.AvgHydropathy]];
fullPercB = [[oneTwoThreeMatchBS.AvgPercB],[oneTwoThreeMatchBB.AvgPercB],[oneTwoThreeMatchOS.AvgPercB],[oneTwoThreeMatchOB.AvgPercB]];
[fitoutpoly, gof] = fit(transpose(hydros), transpose(fullPercB), 'poly2');
scatter(hydros,fullPercB);
hold on
xline(mean([1.4036, 1.4240]), '--', 'LineWidth', 0.75);
xline(mean([1.6845,1.6659]),'LineWidth', 0.75);
legend('% Sheet Content',  'Avg Hydro, Short Seg', 'Avg Hydro, Long Seg');
plot(fitoutpoly);
xlabel('Average Hydropathy','FontWeight', 'bold');
ylabel('Average % Sheet Content','FontWeight', 'bold');
title('Hydropathy vs % Sheet Content, R^2 = ' + string(gof.rsquare));
hold off
sgtitle('Proteins Containing All Three Segments');
%% Just Small and Just Big, Same results as above, not necessary
subplot(2,2,1);
hydros = [[oneTwoThreeMatchBS.AvgHydropathy],[oneTwoThreeMatchOS.AvgHydropathy]];
dist12 = [[oneTwoThreeMatchBS.SegDist12],[oneTwoThreeMatchOS.SegDist12]];
dist13 = [[oneTwoThreeMatchBS.SegDist13],[oneTwoThreeMatchOS.SegDist13]];
dist23 = [[oneTwoThreeMatchBS.SegDist23],[oneTwoThreeMatchOS.SegDist23]];
[fitoutpoly, gof2] = fit(transpose([hydros,hydros,hydros]), transpose([dist12,dist13,dist23]), 'poly2');
scatter(hydros,dist12, 'filled');
hold on
scatter(hydros,dist13, 'filled');
scatter(hydros,dist23,'^','filled');
plot(fitoutpoly);
hold off
subplot(2,2,2)
hydros = [[oneTwoThreeMatchBS.AvgHydropathy],[oneTwoThreeMatchOS.AvgHydropathy]];
fullPercB = [[oneTwoThreeMatchBS.AvgPercB],[oneTwoThreeMatchOS.AvgPercB]];
[fitoutpoly, gof] = fit(transpose(hydros), transpose(fullPercB), 'poly2');
scatter(hydros,fullPercB);
hold on
plot(fitoutpoly);
hold off
subplot(2,2,3);
hydros = [[oneTwoThreeMatchBB.AvgHydropathy],[oneTwoThreeMatchOB.AvgHydropathy]];
dist12 = [[oneTwoThreeMatchBB.SegDist12],[oneTwoThreeMatchOB.SegDist12]];
dist13 = [[oneTwoThreeMatchBB.SegDist13],[oneTwoThreeMatchOB.SegDist13]];
dist23 = [[oneTwoThreeMatchBB.SegDist23],[oneTwoThreeMatchOB.SegDist23]];
[fitoutpoly, gof3] = fit(transpose([hydros,hydros,hydros]), transpose([dist12,dist13,dist23]), 'poly2');
scatter(hydros,dist12);
hold on
scatter(hydros,dist13, 'r');
scatter(hydros,dist23, 'g');
plot(fitoutpoly);
hold off
subplot(2,2,4)
hydros = [[oneTwoThreeMatchBB.AvgHydropathy],[oneTwoThreeMatchOB.AvgHydropathy]];
fullPercB = [[oneTwoThreeMatchBB.AvgPercB],[oneTwoThreeMatchOB.AvgPercB]];
[fitoutpoly, gof4] = fit(transpose(hydros), transpose(fullPercB), 'poly2');
scatter(hydros,fullPercB);
hold on
plot(fitoutpoly);
hold off
%% For Proteins Containing at least 2 Segments
subplot(1,2,1);
hydros12 = [[oneTwoMatchBS.AvgHydropathy],[oneTwoMatchOS.AvgHydropathy],[oneTwoMatchBB.AvgHydropathy],[oneTwoMatchOB.AvgHydropathy]];
hydros13 = [[oneThreeMatchBS.AvgHydropathy],[oneThreeMatchOS.AvgHydropathy],[oneThreeMatchBB.AvgHydropathy],[oneThreeMatchOB.AvgHydropathy]];
hydros23 = [[twoThreeMatchBS.AvgHydropathy],[twoThreeMatchOS.AvgHydropathy],[twoThreeMatchBB.AvgHydropathy],[twoThreeMatchOB.AvgHydropathy]];
dist12 = [[oneTwoMatchBS.SegDist],[oneTwoMatchOS.SegDist],[oneTwoMatchBB.SegDist],[oneTwoMatchOB.SegDist]];
dist13 = [[oneThreeMatchBS.SegDist],[oneThreeMatchOS.SegDist],[oneThreeMatchBB.SegDist],[oneThreeMatchOB.SegDist]];
dist23 = [[twoThreeMatchBS.SegDist],[twoThreeMatchOS.SegDist],[twoThreeMatchBB.SegDist],[twoThreeMatchOB.SegDist]];
[fitoutpoly, gof] = fit(transpose([hydros12,hydros13,hydros23]), transpose([dist12,dist13,dist23]), 'poly2');
scatter(hydros12,dist12);
hold on
c = [0.8500 0.3250 0.0980];
scatter(hydros13,dist13, '^', 'MarkerEdgeColor', c);
c = [0.9290 0.6940 0.1250];
scatter(hydros23,dist23,'+','MarkerEdgeColor', c);
xline(mean([1.4036, 1.4240]), '--', 'LineWidth', 0.75);
xline(mean([1.6845,1.6659]),'LineWidth', 0.75);
legend('Seg 1 + Seg 2', 'Seg 1 + Seg 3','Seg 2 + Seg 3', 'Avg Hydro, Short Seg', 'Avg Hydro, Long Seg', 'Location','northeast');
plot(fitoutpoly);
xlabel('Average Hydropathy','FontWeight', 'bold');
ylabel('Distance Between Segments','FontWeight', 'bold');
title('Hydropathy vs RMSD between segments, R^2 = ' + string(gof.rsquare));
hold off
subplot(1,2,2)
sheet12 = [[oneTwoMatchBS.AvgPercB],[oneTwoMatchOS.AvgPercB],[oneTwoMatchBB.AvgPercB],[oneTwoMatchOB.AvgPercB]];
sheet13 = [[oneThreeMatchBS.AvgPercB],[oneThreeMatchOS.AvgPercB],[oneThreeMatchBB.AvgPercB],[oneThreeMatchOB.AvgPercB]];
sheet23 = [[twoThreeMatchBS.AvgPercB],[twoThreeMatchOS.AvgPercB],[twoThreeMatchBB.AvgPercB],[twoThreeMatchOB.AvgPercB]];
[fitoutpoly, gof] = fit(transpose([hydros12,hydros13,hydros23]), transpose([sheet12,sheet13,sheet23]), 'poly2');
scatter(hydros12,sheet12);
hold on
c = [0.8500 0.3250 0.0980];
scatter(hydros13,sheet13,'^','MarkerEdgeColor', c);
c = [0.9290 0.6940 0.1250];
scatter(hydros23,sheet23,'+','MarkerEdgeColor', c);
xline(mean([1.4036, 1.4240]), '--', 'LineWidth', 0.75);
xline(mean([1.6845,1.6659]),'LineWidth', 0.75);
legend('Seg 1 + Seg 2', 'Seg 1 + Seg 3','Seg 2 + Seg 3', 'Avg Hydro, Short Seg', 'Avg Hydro, Long Seg', 'Location','northeast');
plot(fitoutpoly);
xlabel('Average Hydropathy','FontWeight', 'bold');
ylabel('Average % Sheet Content','FontWeight', 'bold');
title('Hydropathy vs % Sheet Content, R^2 = ' + string(gof.rsquare));
axis([min([hydros12,hydros13,hydros23]),max([hydros12,hydros13,hydros23]),0,1]);
hold off
sgtitle('Proteins Containing at Least 2 Segments');
%% Average Hydropathy vs SD of Hydropathy 
hydros = [[oneTwoThreeMatchBS.AvgHydropathy],[oneTwoThreeMatchBB.AvgHydropathy],[oneTwoThreeMatchOS.AvgHydropathy],[oneTwoThreeMatchOB.AvgHydropathy]];
SDhydros = [[oneTwoThreeMatchBS.SDHydropathy],[oneTwoThreeMatchBB.SDHydropathy],[oneTwoThreeMatchOS.SDHydropathy],[oneTwoThreeMatchOB.SDHydropathy]];
scatter(hydros, SDhydros, 'filled');
hold on 
[fitoutpoly1, gof1] = fit(transpose(hydros),transpose(SDhydros), 'poly1');
[fitoutpoly2, gof2] = fit(transpose(hydros),transpose(SDhydros), 'poly2');
plot(fitoutpoly1, 'm')
plot(fitoutpoly2);
xline(mean([1.4036, 1.4240]), '--', 'LineWidth', 0.75);
xline(mean([1.6845,1.6659]),'LineWidth', 0.75);
xlabel('Average Hydropathy');
ylabel('StD. of Average Hydropathy');
axis([1.4,2.4,0,2]);
title("Avg vs StD Hydropathy, Poly1 R^2 = " + string(gof1.rsquare) + ', Poly2 R^2 = ' + string(gof2.rsquare));
legend('Location','northwest');
hold off;
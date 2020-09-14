mainFolder = 'D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Mississippi_River_11_years_20200511\Kelin_model\NEW_BR_2009_2019\BR_New_CS\';
saveFolder = 'D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Mississippi_River_11_years_20200511\equivalentTrapoCS_multipleType\';
fileNamePart = 'Br_';
step = 10;
crossSecLowerLim=(step-1)*28+1;
crossSecHigherLim=step*28;

for i = crossSecLowerLim:crossSecHigherLim
    crossSection = load([mainFolder fileNamePart sprintf('%04d', i) '.dat']);
    lines = load([mainFolder fileNamePart sprintf('%04d', i) '_lines']);
    tables = importCStable([mainFolder fileNamePart sprintf('%04d', i) '_tab'],2,52);
    CStable = table2array(tables);
    CStable(:,1) = CStable(:,1) - min(CStable(:,1));
    allCStable(:,:,i)=CStable;

    figure(1);
    plot(crossSection(:,1),crossSection(:,2)-min(crossSection(:,2)),'k-');
    hold on;
%     for j= 1:length(lines)/3
%         plot(lines( (j-1)*3+1: (j-1)*3+3 ,1),lines((j-1)*3+1: (j-1)*3+3,2),'g-');
%         pause(0.1);
%     end
    
%     hold off;
    bottom(i) = min(crossSection(:,2));
    figure(2);
    plot(CStable(:,1),CStable(:,2),'b');
    hold on;
end

averagTable=mean(allCStable,3);
plot(averagTable(:,1),averagTable(:,2),'r--','linewidth',3);
z = 15.5;
Bw = 80;
Tw = 850;
TwCC=850;
bfd =  (Tw - Bw)/(2.0*z);
for i = 1:length(averagTable(:,1))
    if averagTable(i,1)<=bfd
        trapiArea(i) = (Bw + averagTable(i,1) * z) * averagTable(i,1);
    else
        trapiArea(i) =  (Bw + bfd * z) * bfd + TwCC * (averagTable(i,1) -bfd);
    end
end
plot(averagTable(:,1),trapiArea,'k','linewidth',3);
 ylabel('Area (m^{2})')
xlabel('Depth (m)')
title(['From CS ' num2str(crossSecLowerLim) ' to ' num2str(crossSecHigherLim)])
    figure(2);
    hold off;
    figure(1);
    hold off;


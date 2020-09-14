mainFolder = 'D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Mississippi_River_11_years_20200511\Kelin_model\NEW_BR_2009_2019\BR_New_CS\';
saveFolder = 'D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Mississippi_River_11_years_20200511\equivalentTrapoCS_milderSlope\';
fileNamePart = 'Br_';

for i = 1:280
    crossSection = load([mainFolder fileNamePart sprintf('%04d', i) '.dat']);
    lines = load([mainFolder fileNamePart sprintf('%04d', i) '_lines']);
    tables = importCStable([mainFolder fileNamePart sprintf('%04d', i) '_tab'],2,52);
    CStable = table2array(tables);
    CStable(:,1) = CStable(:,1) - min(CStable(:,1));
    allCStable(:,:,i)=CStable;

    figure(1);
%     hold on;
%     plot(crossSection(:,1),crossSection(:,2)-min(crossSection(:,2)),'k-');
%     hold on;
%     for j= 1:length(lines)/3
%         plot(lines( (j-1)*3+1: (j-1)*3+3 ,1),lines((j-1)*3+1: (j-1)*3+3,2),'g-');
%         pause(0.1);
%     end
    
%     hold off;
    bottom(i) = min(crossSection(:,2));
%     figure(2);
%     plot(CStable(:,1),CStable(:,2),'b');
%     hold on;
end
averagTable=mean(allCStable,3);
plot(averagTable(:,1),averagTable(:,2),'r--','linewidth',5);
z = 15;
Bw = 50;
Tw = 800;
TwCC=1200;
bfd =  (Tw - Bw)/(2.0*z);
for i = 1:length(averagTable(:,1))
    if averagTable(i,1)<=bfd
        trapiArea(i) = (Bw + averagTable(i,1) * z) * averagTable(i,1);
    else
        trapiArea(i) =  (Bw + bfd * z) * bfd + TwCC * (averagTable(i,1) -bfd);
    end
end
plot(averagTable(:,1),trapiArea,'k','linewidth',5);
 ylabel('Area (m^{2})')
xlabel('Depth (m)')

xx = load('D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Mississippi_River_11_years_20200511\Kelin_model\NEW_BR_2009_2019\input\dx4.txt');

generalGradient = (max(bottom) - min(bottom))/xx(end,1);

for i = 1:length(xx)
    newBottom(i)=max(bottom)-generalGradient*xx(i,1);
end
figure(3)
plot(xx(:,1),bottom,'k');
hold on;
% plot(xx(:,1),newBottom,'r');
generalGradient = (-24 - (-30.48))/xx(end,1);

for i = 1:length(xx)
    newBottom(i)=-24-generalGradient*xx(i,1);
end

plot(xx(:,1),newBottom,'b','lineWidth',1.5);
hold off;
legend('Original channel bottom','New applied channel bottom');
ylabel('Bottom Elevation (m^{2})')
xlabel('Distance from U/S (m)')

newX = [-0.5*TwCC -0.5*TwCC -0.5*Tw -0.5*Bw 0.5*Bw 0.5*Tw 0.5*TwCC 0.5*TwCC]+0.5*TwCC;
newY = [2*bfd     bfd       bfd     0       0      bfd    bfd      2*bfd   ];
figure(1)
hold on;
plot(newX,newY,'r');
xlabel('Distance (m)');
ylabel('Depth (m)');

for i = 1:280
    
   fileNameNew = [saveFolder fileNamePart sprintf('%04d', i) ];
   fid = fopen(fileNameNew,'w');
   fprintf(fid,'%s\t%s\n','x','y');
   for j = 1:length(newX)
       fprintf(fid,'%12.3f\t%12.3f\n',newX(j),newY(j)+newBottom(i));
   end
   fclose(fid);
end


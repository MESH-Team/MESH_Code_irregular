mainFolder = 'D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Mississippi_River_11_years_20200511\Kelin_model\NEW_BR_2009_2019\BR_New_CS\';
saveFolder = 'D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Mississippi_River_11_years_20200511\equivalentTrapoCS_multipleType\';
fileNamePart = 'Br_';

xx = load('D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Mississippi_River_11_years_20200511\Kelin_model\NEW_BR_2009_2019\input\dx4.txt');

% % % z = 13;	Bw = 40;	Tw = 700;	TwCC=1200;
% % % z = 13.5;	Bw = 50;	Tw = 700;	TwCC=1050;
% % % z = 13.7;	Bw = 50;	Tw = 720;	TwCC=1020;
% % % z = 14;	Bw = 55;	Tw = 750;	TwCC=1000;
% % % z = 14.2;	Bw = 55;	Tw = 800;	TwCC=900;
% % % z = 14.8;	Bw = 60;	Tw = 820;	TwCC=880;
% % % z = 15.2;	Bw = 70;	Tw = 820;	TwCC=820;
% % % z = 15.4;	Bw = 70;	Tw = 860;	TwCC=880;
% % % z = 15;	Bw = 80;	Tw = 850;	TwCC=850;
% % % z = 15.5;	Bw = 80;	Tw = 850;	TwCC=850;

for i = 1:280
    crossSection = load([mainFolder fileNamePart sprintf('%04d', i) '.dat']);
    lines = load([mainFolder fileNamePart sprintf('%04d', i) '_lines']);
    tables = importCStable([mainFolder fileNamePart sprintf('%04d', i) '_tab'],2,52);
    CStable = table2array(tables);
    CStable(:,1) = CStable(:,1) - min(CStable(:,1));
    allCStable(:,:,i)=CStable;

%     figure(1);
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

generalGradient = (max(bottom) - min(bottom))/xx(end,1);

for i = 1:length(xx)
    newBottom(i)=max(bottom)-generalGradient*xx(i,1);
end
figure(3)
plot(xx(:,1),bottom,'k');
hold on;
plot(xx(:,1),newBottom,'r');
generalGradient = (-24 - (-30.48))/xx(end,1);

for i = 1:length(xx)
    newBottom(i)=-24-generalGradient*xx(i,1);
end
plot(xx(:,1),newBottom,'b');
hold off;
legend('Original channel bottom','Old applied channel bottom', 'New applied channel bottom');
ylabel('Bottom Elevation (m^{2})')
xlabel('Distance from U/S (m)')

BwAll  =  [40	50	50	55	55	60	70	70	80	80];
TwAll  = [700	700	720	750	800	820	820	860	850	850];
zAll   = [13	13.5	13.7	14	14.2	14.8	15.2	15.4	15	15.5];
TwCCAll= [1200	1050	1020	1000	900	880	820	880	850	850];

maxStep = 10;
lgd = cell(maxStep,1);

for step = 1:maxStep
crossSecLowerLim=(step-1)*28+1;
crossSecHigherLim=step*28;

Bw = BwAll(step);
Tw = TwAll(step);
z  = zAll(step);
TwCC = TwCCAll(step);

bfd =  (Tw - Bw)/(2.0*z);
% for i = 1:length(averagTable(:,1))
%     if averagTable(i,1)<=bfd
%         trapiArea(i) = (Bw + averagTable(i,1) * z) * averagTable(i,1);
%     else
%         trapiArea(i) =  (Bw + bfd * z) * bfd + TwCC * (averagTable(i,1) -bfd);
%     end
% end
% plot(averagTable(:,1),trapiArea,'k','linewidth',3);
%  ylabel('Area (m^{2})')
% xlabel('Depth (m)')
% title(['From CS ' num2str(crossSecLowerLim) ' to ' num2str(crossSecHigherLim)])
%     figure(2);
%     hold off;
%     figure(1);
%     hold off;








newX = [-0.5*TwCC -0.5*TwCC -0.5*Tw -0.5*Bw 0.5*Bw 0.5*Tw 0.5*TwCC 0.5*TwCC]+0.5*TwCC;
newY = [2*bfd     bfd       bfd     0       0      bfd    bfd      2*bfd   ];
figure(1)
hold on;
plot(newX,newY);
lgd{step}=['CS no ' num2str(crossSecLowerLim)];
xlabel('Distance (m)');
ylabel('Depth (m)');

% STOP

for i = crossSecLowerLim:crossSecHigherLim
    
   fileNameNew = [saveFolder fileNamePart sprintf('%04d', i) ];
   fid = fopen(fileNameNew,'w');
   fprintf(fid,'%s\t%s\n','x','y');
   for j = 1:length(newX)
       fprintf(fid,'%12.3f\t%12.3f\n',newX(j),newY(j)+newBottom(i));
   end
   fclose(fid);
end

end
figure(1);
legend(lgd);

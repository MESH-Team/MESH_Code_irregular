%clear; close all;

resultReferenceDate = datenum(2009,1,1);

gage_name={'Baton Rouge (RM 228.4)','Donaldsonville (RM 173.6)','Reserve (RM 138.7)','Bonne Carre - North of Spillway (RM 129.2)', ...
           'Bonne Carre (RM 126.9)','New Orleans (RM  102.8)','Harvey Lock (RM 98.3)', 'IHNC Lock (RM 92.7)','Algiers Lock (RM 88.3)'};
file_name={'Baton_rouge','Donaldsonville','Reserve','BC_north', ...
           'BC','NO','Harvey', 'IHNC','Algiers'};
datum_off=[-0.71 -0.88 -0.76 -0.55 -0.7 -0.82 -0.83 -0.11 -0.11];

node_id=[1	73	117	128	131	156	161	166	171	185	198	211	232	252	258	263	269	280];

xxtick=[datenum(2009,1,1) datenum(2010,1,1) datenum(2011,1,1) datenum(2012,1,1) datenum(2013,1,1) datenum(2014,1,1)...
        datenum(2015,1,1) datenum(2016,1,1) datenum(2017,1,1) datenum(2018,1,1) datenum(2019,1,1) datenum(2020,1,1)];
xxticklabel={'','2010','','2012','','2014','','2016','','2018','','2020'};

f2m=0.3048;
resultFOlder = 'C:\Users\mbeg\Desktop\Mississippi_River_11_years_20200511\output_Diffusive\';

mesh=load([resultFOlder 'output_wl.txt']);

  %% Keeping only the desired time range
 timeMin = datenum(2009,01,01);
 timeMax = datenum(2020,01,01);
 
 mesh = mesh ( resultReferenceDate+mesh(:,1)/3600/24 >= timeMin,: );
 mesh = mesh ( resultReferenceDate+mesh(:,1)/3600/24 <= timeMax,: );
 
 %calculating the date time from the results
 cal_day=resultReferenceDate+mesh(:,1)/3600/24;
 
 
 observed = ones(length(cal_day),length(node_id)+1)*NaN;
 observed(:,1)=cal_day;
 simulated = observed;
 
 %% 
dataFolder = 'C:\Users\mbeg\Desktop\Mississippi_River_11_years_20200511\Observed_Data\USACE\';
figure(7)
for k=1:9
    subplot(6,3,k)
  % val=load(['D:\Klhu_work\projects\MP2023\data\USACE\STG_' file_name{k} '_2009_2019.txt']);
   val=load([dataFolder 'STG_' file_name{k} '_2009_2019.txt']);
   tob=datenum(val(:,3),val(:,1),val(:,2));
   wlob=(val(:,4)+datum_off(k))*f2m;
   
   if k==1
       h1=plot(tob,wlob,'r');
   else
       plot(tob,wlob,'r');
   end
   
   D = wlob( ismember( tob, cal_day ) );
   E = ismember( cal_day, tob );
   observed(E,k+1)=D;

   
   hold on
   
      datetick('x','yyyy')
   xlim([datenum(2009,1,1) datenum(2020,1,1)])
   
   ylim([0 15])
   if(k>3)
       ylim([0 8])
   end
    title(gage_name{k});
    ylabel('WL (m), NAVD88)')
   clear val tob wlob
   set(gca,'xtick',xxtick,'xticklabel',xxticklabel);

end

% Belle Chasse
val=load('C:\Users\mbeg\Desktop\Mississippi_River_11_years_20200511\Observed_Data\USGS\BC_2009_2019_raw.dat');


time_BC=datenum(val(:,1),val(:,2),val(:,3));
dis_BC=val(:,5)*f2m^3;
wl_BC=(val(:,4)-6.58)*f2m;

subplot(6,3,10)
plot(time_BC,wl_BC,'r')


D = wl_BC( ismember( time_BC, cal_day ) );
E = ismember( cal_day, time_BC );
observed(E,11)=D;

   hold on
   xlim([datenum(2009,1,1) datenum(2020,1,1)]);
   
   
   ylim([-.5 5])
    title('Belle Chasse (RM 76.4)');
    ylabel('WL (m), NAVD88)')
   set(gca,'xtick',xxtick,'xticklabel',xxticklabel);
   
   
%%%%%%%
gage_name2={'Alliance (RM 62.5)','West Point a la Hache (RM 48.7)','Empire (RM 29.5)','Venice (RM 10.7)', ...
           'West Bay (RM 4.7)','Head of Pass (RM -0.6)','Southwest Pass (BHP 7.5)', 'East Jetty (BHP 17.9)'};
file_name2={'Alliance','WPLH','Empire','Venice', ...
           'WestBay','HoP','SWP7.5', 'East_Jetty'};
datum_off1=[0  -0.21 -0.6 -1.84 -0.923 -0.853 -0.958 -0.896];
datum_off2=[0  0 -0.12 0 0 0 0 0];
datum_off3=[0  0 0 0 -0.26 0 -0.26 0];
datum_time1=[datenum(2000,1,1) datenum(2015,7,1) datenum(2015,7,1) datenum(2015,6,28) ...
             datenum(2015,6,28) datenum(2015,6,30) datenum(2015,6,29) datenum(2015,6,29)];
datum_time2=[datenum(2100,1,1) datenum(2100,1,1) datenum(2018,9,26) datenum(2100,1,1) ...
             datenum(2100,1,1) datenum(2100,1,1) datenum(2100,1,1) datenum(2100,1,1)];
datum_time3=[datenum(2000,1,1) datenum(2000,1,1) datenum(2000,1,1) datenum(2000,1,1) ...
             datenum(2009,1,30) datenum(2000,1,1) datenum(2010,8,19) datenum(2000,1,1)];



for k=1:8
    subplot(6,3,k+10)
    
   val=load([dataFolder 'STG_' file_name2{k} '_2009_2019.txt']);
%    val=load(['D:\Klhu_work\projects\MP2023\data\USACE\STG_' file_name2{k} '_2009_2019.txt']);
   tob=datenum(val(:,3),val(:,1),val(:,2));
   wlob=val(:,4);
   wlob(tob<datum_time1(k))=wlob(tob<datum_time1(k))+datum_off1(k);
   wlob(tob>=datum_time2(k))=wlob(tob>=datum_time2(k))+datum_off2(k);
   wlob(tob<datum_time3(k))=wlob(tob<datum_time3(k))+datum_off3(k);

   wlob=wlob*f2m;
   
   plot(tob,wlob,'r')
   
   D = wlob( ismember( tob, cal_day ) );
   E = ismember( cal_day, tob );
   observed(E,k+11)=D;
   
   hold on
   xlim([datenum(2009,1,1) datenum(2020,1,1)])
   
   
   ylim([-.5 2])
   if(k<=2)
       ylim([-0.5 5])
   end
    title(gage_name2{k});
    ylabel('WL (m), NAVD88)')
   set(gca,'xtick',xxtick,'xticklabel',xxticklabel);

   clear val tob wlob

end

% plot mesh results

for k=1:18
    subplot(6,3,k)
   if k==1
       h2=plot(cal_day,mesh(:,node_id(k)+1),'k');
   else
       plot(cal_day,mesh(:,node_id(k)+1),'k');
   end
    
    simulated(:,k+1)=mesh(:,node_id(k)+1);
%     xlim([datenum(2018,1,1) datenum(2019,1,1)])
end
subplot(6,3,1);
legend([h1 h2],'USGS Observed','MESH Simulated')
lgd.FontSize = 10;
% sgtitle('Datum in NAVD88');


% plotting Q
mesh=load([resultFOlder 'q.txt']);
figure(10)
h3 = plot(time_BC,dis_BC,'ro');
hold on;
k=10; % location of node id for Belle Chasse
h4 = plot(cal_day,mesh(:,node_id(k)+1),'k');
legend([h3 h4],'USACE Observed','MESH Simulated');
lgd.FontSize = 10;
ylabel('Discharge (m^{3}/s)');
set(gca,'xtick',xxtick,'xticklabel',xxticklabel);
xxticklabel={'','2010','','2012','','2014','','2016','','2018','','2020'};
xlim([datenum(2009,1,1) datenum(2020,1,1)]);
title(gage_name{k});

mainFolder =  'D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Multijunction_Network\';
caseFolder1 = 'Large_NHD_Result\';
% caseFolder2 = 'DT_var\';
caseFolder3 = 'Large_NHD_Geometry\';
videoNeeded = 1;
colorPallet = jet;
maxVal = 40;
minVal = 0;
valLine=linspace(minVal,maxVal,length(colorPallet(:,1)));
% hh = colorbar;
% set(hh, 'ylim', [minVal maxVal])

nlinks = 29;
nx1 = [3	3	2	3	2	3	4	3	5	3	2	3	2	2	4	2	3	3	4	2	2	2	3	2	2	2	2	2	2];
startPointX = [-75.4018021	-75.3983002	-75.4139023	-75.43	-75.4097977	-75.3824997	-75.4023972	-75.3770981	-75.3554993	-75.3676987	-75.3871002	-75.3813019	-75.3998032	-75.4030991	-75.4057999	-75.4330978	-75.4448013	-75.4425964	-75.4600983	-75.4720993	-75.4757004	-75.4850998	-75.5006027	-75.4933014	-75.4833984	-75.4793015	-75.4590988	-75.4826965	-75.4934006];
startPointY = [41.2523003	41.2412987	41.233799	41.23	41.2244987	41.239399	41.217701	41.2127991	41.2178993	41.1790009	41.1850014	41.1946983	41.1851997	41.1904984	41.1828995	41.1778984	41.1871986	41.1649017	41.1609001	41.1867981	41.1610985	41.1655006	41.1419983	41.1543999	41.1529007	41.1554985	41.1417999	41.1425018	41.1277008];
% channel 4 start point X and start point y has been changed for plotting
% purpose, to avoind overlapping. Original x = -75.4160995, original 
% y = 41.2503014

endPointX =   [-75.4139023	-75.4139023	-75.4097977	-75.4097977	-75.4023972	-75.4023972	-75.4030991	-75.4030991	-75.3871002	-75.3871002	-75.3998032	-75.3998032	-75.4057999	-75.4057999	-75.4425964	-75.4425964	-75.4600983	-75.4600983	-75.4757004	-75.4757004	-75.4793015	-75.4793015	-75.4833984	-75.4833984	-75.4826965	-75.4826965	-75.4934006	-75.4934006	-75.499];
endPointY =   [41.233799	41.233799	41.2244987	41.2244987	41.217701	41.217701	41.1904984	41.1904984	41.1850014	41.1850014	41.1851997	41.1851997	41.1828995	41.1828995	41.1649017	41.1649017	41.1609001	41.1609001	41.1610985	41.1610985	41.1554985	41.1554985	41.1529007	41.1529007	41.1425018	41.1425018	41.1277008	41.1277008	41.126];


dxAll = zeros(max(nx1),nlinks);

zAll=[   635.830017,626.229980,610.369995,0.00000000,0.00000000;
634.570007,626.859985,610.369995,0.00000000,0.00000000;
610.369995,607.109985,0.00000000,0.00000000,0.00000000;
613.580017,613.580017,607.109985,0.00000000,0.00000000;
607.109985,597.340027,0.00000000,0.00000000,0.00000000;
628.900024,627.719971,597.340027,0.00000000,0.00000000;
597.340027,595.719971,595.080017,589.219971,0.00000000;
633.820007,626.580017,589.219971,0.00000000,0.00000000;
635.320007,635.299988,606.539978,603.400024,597.270020;
607.830017,601.219971,597.270020,0.00000000,0.00000000;
597.270020,589.200012,0.00000000,0.00000000,0.00000000;
603.099976,603.099976,589.200012,0.00000000,0.00000000;
589.200012,582.849976,0.00000000,0.00000000,0.00000000;
589.219971,582.849976,0.00000000,0.00000000,0.00000000;
582.849976,578.979980,578.599976,555.320007,0.00000000;
581.070007,555.320007,0.00000000,0.00000000,0.00000000;
584.429993,580.080017,536.239990,0.00000000,0.00000000;
555.320007,548.119995,536.239990,0.00000000,0.00000000;
536.239990,532.710022,530.960022,523.289978,0.00000000;
581.679993,523.289978,0.00000000,0.00000000,0.00000000;
523.289978,515.859985,0.00000000,0.00000000,0.00000000;
546.739990,515.859985,0.00000000,0.00000000,0.00000000;
541.969971,536.570007,521.479980,0.00000000,0.00000000;
539.859985,521.479980,0.00000000,0.00000000,0.00000000;
521.479980,511.839996,0.00000000,0.00000000,0.00000000;
515.859985,511.839996,0.00000000,0.00000000,0.00000000;
548.119995,509.920013,0.00000000,0.00000000,0.00000000;
511.839996,509.920013,0.00000000,0.00000000,0.00000000;
509.920013,497.390015,0.00000000,0.00000000,0.00000000;];

zAll=zAll*0;

Q = load([mainFolder caseFolder1 'q.txt']);

[a,b]= size(Q);

for j = 1:nlinks
%     dx = load([mainFolder caseFolder3 'dx' num2str(j)]);
    dx = load([mainFolder caseFolder3 'dx' num2str(j,'%04d')]);
    dx = dx(:,2);
    ncomp=nx1(j);
    for i = 1:ncomp
        dxAll(i,j)=dx(i);
    end
end

t = 1;
% dx = load([mainFolder caseFolder3 'dx']);
figure (1);
hold on;
for j = 1:nlinks
    ncomp=nx1(j);
    XX = createLinespace(startPointX(j),endPointX(j),dxAll(1:ncomp-1,j));
    YY = createLinespace(startPointY(j),endPointY(j),dxAll(1:ncomp-1,j));

        hlines{j} = plot(XX,YY);

    colorbar();colormap(colorPallet)
    set(gca, 'CLim', [minVal, maxVal]);
    colorTitleHandle = get(colorbar,'Title');
%     xlim([-5000 15000]); ylim([0 20000]);
    set(colorTitleHandle ,'String','Q (m^{3}/s)');
%     set(colorTitleHandle ,'String','Depth (m)');
end
hold off;

ii = 1;

for t=1:a/nlinks
    for j=1:nlinks
        ncomp=nx1(j);
        colrReq = Q((t-1)*nlinks+j,3:ncomp+2)-zAll(j,1:ncomp);
%         if j == 1
            for i = 1:ncomp-1
                val = (colrReq(i)+colrReq(i+1))/2;
                %rr = interp1(valLine,colorPallet,50);
                if val<minVal
                    temp = interp1(valLine,colorPallet,maxVal-(maxVal-minVal)*1.02,'linear','extrap');
                elseif val>maxVal
                    temp = interp1(valLine,colorPallet,minVal+(maxVal-minVal)*1.02,'linear','extrap');
                else
                    temp = interp1(valLine,colorPallet,val);
                end
                temp(temp<0)=0; temp(temp>1)=1;
                hlines{j}(i).Color(1:3) = temp;
            end
            [hlines{j}.LineWidth] = deal(8);
    end

    
    title(['Time = ' num2str(Q((t-1)*nlinks+1,1)/3600,'%2.2f') ' hr']);


    Fra(ii)= getframe(gcf);
    ii = ii + 1;
end

if videoNeeded ==1
    
   writerObj = VideoWriter([mainFolder caseFolder1 'Map_Animation_Q.avi']);
   writerObj.FrameRate = 10;
  open(writerObj);
    % write the frames to the video
    for j=2:length(Fra)-1
        % convert the image to a frame
        frame = Fra(j) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);

end



        



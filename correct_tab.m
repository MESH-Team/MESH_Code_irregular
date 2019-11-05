caseFolder = 'D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Sim2\CS\';
nel=101;
ncomp=41;
for i = 1:ncomp
    filename=strcat('Sim2_', sprintf('%03d',i));

    filename=[caseFolder filename '_tab'];
    dataTable=importfile_nab(filename,2,nel+1);
    pp=table2array(dataTable);
    pp(1,6)=pp(2,6);

    fileID = fopen(filename, 'w');
    fprintf(fileID,'%s\n','Elev(m)    Area(m2)     Peri(m)      Radi(m)   Conv(m3/s)    topWidth(m)    newI1(m3)    dPdA(1/m)');
    for j=1:nel
        fprintf(fileID,'%f %f %f %f %f %f %f %f\n',pp(j,1),pp(j,2),pp(j,3),pp(j,4),pp(j,5),pp(j,6),pp(j,7),pp(j,8));
    end
    fclose(fileID);

end
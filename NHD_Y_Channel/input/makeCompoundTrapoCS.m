folderName = 'D:\Project_Works\JTTI\codeTest_rectangle\compoundTrapoCS\';
saveFolder = folderName;

BwAll = load([folderName 'bw.txt']);
TwAll = load([folderName 'Tw.txt']);
TwCCAll = load([folderName 'TwCC.txt']);
zAll = load([folderName 'Bed.txt']);

fileNameNew = [saveFolder '\Bank.txt' ];
fid2 = fopen(fileNameNew,'w');

for i=1:length(zAll)
    Bw = BwAll(i);
    Tw = TwAll(i);
    bfd = 5;
    TwCC = TwCCAll(i);
    bed = zAll(i);

    newX = [-0.5*TwCC -0.5*TwCC -0.5*Tw -0.5*Bw 0.5*Bw 0.5*Tw 0.5*TwCC 0.5*TwCC] +0.5*TwCC;
    newY = [4*bfd     bfd       bfd     0       0      bfd    bfd      4*bfd   ];
    
    fileNameNew = [saveFolder '\Test_' sprintf('%04d', i) ];
    fid = fopen(fileNameNew,'w');
    fprintf(fid,'%s\t%s\n','x','y');
    for k = 1:length(newX)
       fprintf(fid,'%12.3f\t%12.3f\n',newX(k),newY(k)+bed);
    end
    fclose(fid);
    
    fprintf(fid2,'%12.3f\t%12.3f\t%12.3f\n',i, (TwCC-Tw)/2,(TwCC-Tw)/2+Tw);
    
end
fclose(fid2);

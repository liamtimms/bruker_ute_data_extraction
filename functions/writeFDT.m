
%----------------------------------------------
%writeFDT
%This function writes an image volume in FDT file format.
%example: writeFDT('DTI.fdt', A)
%A can be of size 1 for 1D data, 2 for a 2D image,3 for a volume, 4 for a set of volume images.
%
%----------------------------------------------
%
%Downloaded from Angelos Barmpoutis' web-page.
%
%This program is freely distributable without 
%licensing fees and is provided without guarantee 
%or warrantee expressed or implied. This program 
%is NOT in the public domain. 
%
%Copyright (c) 2007 Angelos Barmpoutis. 
%
%VERSION : 20070730
%
%-----------------------------------------------
function writeFDT(name,data)
sz=size(data);
sz1=sz(1);

if length(sz)>=2
    sz2=sz(2);
else
    sz2=1;
end

if length(sz)>=3
    sz3=sz(3);
else
    sz3=1;
end

if length(sz)>=4
    sz4=sz(4);
else
    sz4=1;
end

fid=fopen(name,'w','b');
sum=1;

fwrite(fid,sz1,'int');
fwrite(fid,sz2,'int');
fwrite(fid,sz3,'int');
fwrite(fid,sz4,'int');


for i=1:sz4
   for j=1:sz3
      for x=1:sz2
          fwrite(fid,data(:,x,j,i),'float');
      end
   end
end

fclose(fid);
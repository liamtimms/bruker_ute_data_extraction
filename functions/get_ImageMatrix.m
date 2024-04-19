function [Image, dim] = get_ImageMatrix(pathname,scanNum, recoNum, NumImages, Matrix_size, Slice_orientation, datatype)

%Liam Timms NEU, 7/19/2017

dseqFile = strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\2dseq');
fid2dseq = fopen(dseqFile,'rb');
if strcmp('axial',Slice_orientation)==1
    Image=zeros([Matrix_size(1) Matrix_size(2) Matrix_size(3)*NumImages]);
    dim=[4 [Matrix_size(1) Matrix_size(2) Matrix_size(3)*NumImages] 1 1 1 1];
    for i=1:1:Matrix_size(3)*NumImages
        Slice(:,:)=fread(fid2dseq,[Matrix_size(1) Matrix_size(2)],datatype);
        Image(:,:,i)= rot90(Slice,2);
    end
elseif strcmp('sagittal',Slice_orientation)==1
    Image=zeros([Matrix_size(1) Matrix_size(2)*NumImages Matrix_size(3)]);
    dim=[4 [Matrix_size(1) Matrix_size(2)*NumImages Matrix_size(3)] 1 1 1 1];
    for i=1:1:Matrix_size(2)*NumImages
        Slice(:,:)=fread(fid2dseq,[Matrix_size(1) Matrix_size(2)],datatype);
        Image(:,:,i)= rot90(Slice,2);
    end
elseif strcmp('coronal',Slice_orientation)==1
    Image=zeros([Matrix_size(1)*NumImages Matrix_size(2) Matrix_size(3)]);
    dim=[4 [Matrix_size(1)*NumImages Matrix_size(2) Matrix_size(3)] 1 1 1 1];
    for i=1:1:Matrix_size(1)*NumImages
        Slice(:,:)=fread(fid2dseq,[Matrix_size(1) Matrix_size(2)],datatype);
        Image(:,:,i)= rot90(Slice,2);
    end
end

fclose(fid2dseq);

end
function [nii] = get_NiiFromRawData(pathname,scanNum, recoNum, NumImages, RECO_size, NSlices, Slice_orientation, DataType)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Original 2dseq data and build NIFTI %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Matrix_size=RECO_size;
if NSlices~=1
    Matrix_size(3)=NSlices;
end
datatype=DataType;
[Image, dim] = get_ImageMatrix(pathname,scanNum, recoNum, NumImages, Matrix_size, Slice_orientation, datatype);

% Normalize intensity to number of averages since data is stored by
% accumulative signal (i.e. Stot=S1+S2+S3...+SN, where N is numAvg)
% thus, Snew here is Snew=Stot/NumAvg
Image=Image./NumAvg;
% Scale by RG and Slope
Image=Image.*(SLOPE1/RG);

% Assign NIFTI header identifiers
origin = [0 0 0];
dataSaveType = 16; %always use default float32 for storing
voxel_size=RECO_fov./RECO_size;
ss=size(voxel_size);
if ss(1)==2
    voxel_size(3)=SliceDistance;
end
nii=make_nii(Image, voxel_size, origin, dataSaveType);
nii.hdr.dime.dim = dim;


end
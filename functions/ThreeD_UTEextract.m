function [nii] = ThreeD_UTEextract(pathname,scanNum, savenametype,UTE_Count)
%3DUTE_extract pre-processes 3D UTE scans from raw Bruker data File
%   pathname is the filepath with patient name
%   scanNum is the folder number containing UTE3D raw data from Bruker
%   savenametype : 1 for scan number, 2 for progression
%   UTE_Count : used for progression savenametype 2
%   example: ThreeD_UTEextract('C:\...\Patientname.xyz','5',1,1)
z=1
% ACQP FILE: Get RG
scanFileAcqp=strcat(pathname,'\',scanNum,'\acqp');
[fidAcqp temp] = fopen(scanFileAcqp,'rt');
while feof(fidAcqp) == 0
    acqpline = fgetl(fidAcqp);
    RGline=  regexpi(acqpline,'RG=');
    if ~isempty(RGline)
        acqplineSplit = strsplit(acqpline,'=');
        RG = str2num(char(acqplineSplit(2)));
    end
end

% DETERMINE NUMBER OF RECOS
recoData=dir(strcat(pathname,'\',scanNum,'\pdata\'));
recoIndex=[recoData.isdir];
[numRecos,temp]=size(nonzeros(recoIndex));
numRecos=numRecos-2; % there are always 2 additional

% FOR EACH UTE RECO: Export scan and scale by RG and SLOPE
for j=1:1:numRecos
    recoNum=sprintf('%d',j);
    
    % RECO FILE: Get reco type
    scanFileReco=strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\reco');
    [fidReco temp] = fopen(scanFileReco,'rt');
    while feof(fidReco) == 0
        recoline = fgetl(fidReco);
        RecoTypeLine=  regexpi(recoline,'RECO_image_type=');
        if ~isempty(RecoTypeLine)
            RecoLineSplit = strsplit(recoline,'=');
            recoType = char(RecoLineSplit(2)); % COMPLEX_IMAGE MAGNITUDE_IMAGE etc
        end
    end
    
    % RECO FILE: Get RECO_fov in mm (1x3 matrix)
    % RECO_fov./RECO_size is [x_dim y_dim z_dim] in mm
    scanFileReco=strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\reco');
    [fidReco temp] = fopen(scanFileReco,'rt');
    match=0;
    while feof(fidReco) == 0
        recoline = fgetl(fidReco);
        RecoTypeLine=  regexpi(recoline,'RECO_fov=');
        if match == 1;
            match=0; % Reset Match
            RecoLineSplit = strsplit(recoline,' ');
            RECO_fov = [str2num(char(RecoLineSplit))]'; %3x1 matrix
            RECO_fov = RECO_fov.*10; % FOV in mm
        end
        if ~isempty(RecoTypeLine)
            match=1; % Move on to access next line
        end
    end
    
    % RECO FILE: Get RECO_size (matrix size, 3x1 matrix)
    % RECO_fov./RECO_size is [x_dim y_dim z_dim] in mm
    scanFileReco=strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\reco');
    [fidReco temp] = fopen(scanFileReco,'rt');
    match=0;
    while feof(fidReco) == 0
        recoline = fgetl(fidReco);
        RecoTypeLine=  regexpi(recoline,'RECO_size=');
        if match == 1;
            match=0; % Reset Match
            RecoLineSplit = strsplit(recoline,' ');
            RECO_size = [str2num(char(RecoLineSplit))]'; %1x3 matrix
        end
        if ~isempty(RecoTypeLine)
            match=1; % Move on to access next line
        end
    end
    
    % VISU PARS FILE: Get SLOPE1 (possibly SLOPE2 if complex)
    scanFileVisuPars=strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\visu_pars');
    [fidVisuPars temp] = fopen(scanFileVisuPars,'rt');
    match=0;
    while feof(fidVisuPars) == 0
        visuParsline = fgetl(fidVisuPars);
        SLOPEline=  regexpi(visuParsline,'VisuCoreDataSlope=');
        if match == 1
            match=0; % Reset Match
            slopelineSplit = strsplit(visuParsline,' ');
            SLOPE1 = str2num(char(slopelineSplit(1)));
            if strcmp('COMPLEX_IMAGE',recoType)==1
                SLOPE2=str2num(char(slopelineSplit(2)));% Handles complex recos too
            end
        end
        if ~isempty(SLOPEline)
            match=1; % Move on to access next line
        end
    end
    
    % VISU PARS FILE: Get datatype
    scanFileVisuPars=strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\visu_pars');
    [fidVisuPars temp] = fopen(scanFileVisuPars,'rt');
    while feof(fidVisuPars) == 0
        visuParsline = fgetl(fidVisuPars);
        DataLine=  regexpi(visuParsline,'VisuCoreWordType=');
        if ~isempty(DataLine)
            Datatypesplitline = strsplit(visuParsline,'=_');
            datatype=char(Datatypesplitline(2));
            if strcmp('16BIT_SGN_INT',datatype)==1
                datatype = 'bit16';
            elseif strcmp('32BIT_SGN_INT',datatype)==1
                datatype = 'bit32';
            end
            
        end
    end
    
    % METHOD: Get Slice Orientation
    MethodFile=strcat(pathname,'\',scanNum,'\method');
    [fidMethod temp] = fopen(MethodFile,'rt');
    match=0;
    while feof(fidMethod) == 0
        methodLine = fgetl(fidMethod);
        sliceOrient=  regexpi(methodLine,'PVM_SPackArrSliceOrient=');
        if match == 1
            match=0; % Reset Match
            Slice_orientation = char(methodLine); % axial sagittal coronal
        end
        if ~isempty(sliceOrient)
            match=1; % Move on to access next line
        end
    end
    
    % METHOD FILE: Get number of averages
    MethodFile=strcat(pathname,'\',scanNum,'\method');
    [fidMethod temp] = fopen(MethodFile,'rt');
    while feof(fidMethod) == 0
        methodLine = fgetl(fidMethod);
        MethodTypeLine=  regexpi(methodLine,'PVM_NAverages=');
        if ~isempty(MethodTypeLine)
            MethodLineSplit = strsplit(methodLine,'=');
            NumAvg = char(MethodLineSplit(2)); % The number of averages taken
            NumAvg = str2num(NumAvg);
        end
    end
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Original 2dseq data and build NIFTI
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dseqFile = strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\2dseq');
    fid2dseq = fopen(dseqFile,'rb');
    if strcmp('COMPLEX_IMAGE',recoType)==1 % For complex Images real and imaginary both stored in 2dseq
        if strcmp('axial',Slice_orientation)==1
            Image=zeros([RECO_size(1) RECO_size(2) RECO_size(3)*2]);
            dim=[4 [RECO_size(1) RECO_size(2) RECO_size(3)*2] 1 1 1 1];
            for i=1:1:RECO_size(3)*2
                Image(:,:,i)= fread(fid2dseq,[RECO_size(1) RECO_size(2)],datatype);
            end
        elseif strcmp('sagittal',Slice_orientation)==1
            Image=zeros([RECO_size(1) RECO_size(2)*2 RECO_size(3)]);
            dim=[4 [RECO_size(1) RECO_size(2)*2 RECO_size(3)] 1 1 1 1];
            for i=1:1:RECO_size(2)*2
                Image(:,i,:)= fread(fid2dseq,[RECO_size(1) RECO_size(3)],datatype);
            end
        elseif strcmp('coronal',Slice_orientation)==1
            Image=zeros([RECO_size(1)*2 RECO_size(2) RECO_size(3)]);
            dim=[4 [RECO_size(1)*2 RECO_size(2) RECO_size(3)] 1 1 1 1];
            for i=1:1:RECO_size(1)*2
                Image(i,:,:)= fread(fid2dseq,[RECO_size(2) RECO_size(3)],datatype);
            end
        end
    else % For everything but complex images
        Image=zeros(RECO_size); dim=[4 RECO_size 1 1 1 1];
        if strcmp('axial',Slice_orientation)==1
            for i=1:1:RECO_size(3)
                Image(:,:,i)= fread(fid2dseq,[RECO_size(1) RECO_size(2)],datatype);
            end
        elseif strcmp('sagittal',Slice_orientation)==1
            for i=1:1:RECO_size(2)
                Image(:,i,:)= fread(fid2dseq,[RECO_size(1) RECO_size(3)],datatype);
            end
        elseif strcmp('coronal',Slice_orientation)==1
            for i=1:1:RECO_size(1)
                Image(i,:,:)= fread(fid2dseq,[RECO_size(2) RECO_size(3)],datatype);
            end
        end
    end
    fclose(fid2dseq);
    
    % Normalize intensity to number of averages since data is stored by
    % accumulative signal (i.e. Stot=S1+S2+S3...+SN, where N is numAvg)
    % thus, Snew here is Snew=Stot/NumAvg
    Image=Image./NumAvg;
    
    
    % Assign NIFTI header identifiers
    origin = [0 0 0];
    dataSaveType = 16; %always use default float32 for storing
    voxel_size=RECO_fov./RECO_size;
    nii=make_nii(Image, voxel_size, origin, dataSaveType);
    nii.hdr.dime.dim = dim;
    
    
    % Get StudyName for SaveName
    pathsplit = strsplit(pathname,'\');
    [temp nameNum]=size(pathsplit);
    pathsplit=pathsplit(nameNum);
    pathsplit = strsplit(char(pathsplit),'.');
    StudyName=pathsplit(1);
    
    % Saving and Scaling
    
    
    pathsplit=strsplit(pathname, strcat(StudyName,'.'));
    
    if strcmp('COMPLEX_IMAGE',recoType)==1
        if savenametype==1
            SaveName1=strcat(pathsplit(1),StudyName,'_UTE3D_Scan',scanNum,'_','REAL_IMAGE.nii');
            SaveName2=strcat(pathsplit(1),StudyName,'_UTE3D_Scan',scanNum,'_','IMAGINARY_IMAGE.nii');
            SaveName3=strcat(pathsplit(1),StudyName,'_UTE3D_Scan',scanNum,'_',recoType,'.nii');
        elseif savenametype==2
            counter=sprintf('%d',UTE_Count);
            SaveName1=strcat(pathsplit(1),StudyName,'_UTE3D_UTEscan',counter,'_','REAL_IMAGE.nii');
            SaveName2=strcat(pathsplit(1),StudyName,'_UTE3D_UTEscan',counter,'_','IMAGINARY_IMAGE.nii');
            SaveName3=strcat(pathsplit(1),StudyName,'_UTE3D_UTEscan',counter,'_',recoType,'.nii');
        end
        [nii_R, nii_I] = SplitComplex(nii, Slice_orientation); % Split
        nii_R.img = nii_R.img.*(SLOPE1/RG); % Scale
        nii_I.img = nii_I.img.*(SLOPE2/RG);
        nii.img = nii.img.*(SLOPE1/RG);
        
        if savenametype~=3
            save_nii(nii_R,char(SaveName1)) % Save
            save_nii(nii_I,char(SaveName2))
            save_nii(nii,char(SaveName3))
        end
        
    else
        if savenametype==1
            SaveName=strcat(pathsplit(1),StudyName,'_UTE3D_Scan',scanNum,'_',recoType,'.nii');
        elseif savenametype==2
            counter=sprintf('%d',UTE_Count);
            SaveName=strcat(pathsplit(1),StudyName,'_UTE3D_UTEscan',counter,'_',recoType,'.nii');
        end
        nii.img = nii.img.*(SLOPE1/RG);
        
        if savenametype~=3
            save_nii(nii,char(SaveName))
        end
    end
    
    
end % end this reco

end %end reco loop



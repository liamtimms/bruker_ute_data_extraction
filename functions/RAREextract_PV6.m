function [nii, nii_TR1, nii_TR2] = RAREextract_PV6(pathname,scanNum, savenametype, RARE_Count,NSlices)
%RAREextract pre-processes RARE scans from raw Bruker data File
%   pathname is the filepath with patient name
%   scanNum is the folder number containing RARE raw data from Bruker
%   savenametype : 1 for scan number, 2 for progression
%   Count : used for progression savenametype 2
%   example: RAREextract('C:\...\Patientname.xyz','5',1,1)

% while feof(fid) == 0
%     matches=  regexpi(tline,'PVM_SpatResol=');
%     if ~isempty(matches)
%         numDim =  returnNum(tline);
%         tline = fgetl(fid);
%         ImDim= str2num(tline);
%         if numDim ==3
%             avw.hdr.dime.pixdim(1)= 4;
%             avw.hdr.dime.pixdim(2)= ImDim(1);
%             avw.hdr.dime.pixdim(3)= ImDim(2);
%             avw.hdr.dime.pixdim(4)= ImDim(3);
%             avw.hdr.dime.pixdim(5)= 0;
%             NsliceFlag=0; % Data is 3D%
%         end
%         %avw.hdr.dime.vox_units='mm'
%     end
%
% end


nii=[];
z=1
% ACQP FILE: Get RG
scanFileAcqp=strcat(pathname,'\',scanNum,'\acqp');
[fidAcqp temp] = fopen(scanFileAcqp,'rt');
if strcmp(temp,'')
    FileString=fscanf(fidAcqp, '%c');
    RG=get_Param(FileString, 'RG=');
end

% DETERMINE NUMBER OF RECOS
recoData=dir(strcat(pathname,'\',scanNum,'\pdata\'));
recoIndex=[recoData.isdir];
[numRecos,temp]=size(nonzeros(recoIndex));
numRecos=numRecos-2; % there are always 2 additional (according to Codi but wtf)



% FOR EACH RARE RECO: Export scan and scale by RG and SLOPE
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
    
    % VISU PARS FILE: Get SLOPE1 and SLOPE2
    scanFileVisuPars=strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\visu_pars');
    [fidVisuPars temp] = fopen(scanFileVisuPars,'rt');
    match=0;
    if fidVisuPars==-1
        return
    else
        while feof(fidVisuPars) == 0
            visuParsline = fgetl(fidVisuPars);
            SLOPEline=  regexpi(visuParsline,'VisuCoreDataSlope=');
            if match == 1
                match=0; % Reset Match
                slopelineSplit = strsplit(visuParsline,' ');
                SLOPE1 = str2num(char(slopelineSplit(1)));
                SLOPE1
                SLOPE2 = str2num(char(slopelineSplit(2)));
            end
            if ~isempty(SLOPEline)
                match=1; % Move on to access next line
            end
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
    FileString=fscanf(fidMethod, '%c');
    
    ParamName='PVM_SPackArrSliceDistance=';
    SliceDistStr = get_Param(FileString, ParamName);
    Spaces=strfind(SliceDistStr, ' ');
    SliceDistance=str2num(SliceDistStr(Spaces(2)+3:end));
    
    ParamName='PVM_NAverages=';
    NumAvg = get_Param(FileString, ParamName);
    %
    %     while feof(fidMethod) == 0
    %         methodLine = fgetl(fidMethod);
    %         MethodTypeLine=  regexpi(methodLine,'PVM_NAverages=');
    %         if ~isempty(MethodTypeLine)
    %             MethodLineSplit = strsplit(methodLine,'=');
    %             NumAvg = char(MethodLineSplit(2)); % The number of averages taken
    %             NumAvg = str2num(NumAvg);
    %         end
    %     end
    
    
    
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Original 2dseq data and build NIFTI
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NumImages=1;
    Matrix_size=RECO_size;
    Matrix_size(3)=NSlices;
    [Image, dim] = get_ImageMatrix(pathname,scanNum, recoNum, NumImages, Matrix_size, Slice_orientation, datatype);
    
    
    % Normalize intensity to number of averages since data is stored by
    % accumulative signal (i.e. Stot=S1+S2+S3...+SN, where N is numAvg)
    % thus, Snew here is Snew=Stot/NumAvg
    Image=Image./NumAvg;
    
    
    % Assign NIFTI header identifiers
    origin = [0 0 0];
    dataSaveType = 16; %always use default float32 for storing
    voxel_size=RECO_fov./RECO_size;
    voxel_size(3)=SliceDistance;
    nii=make_nii(Image, voxel_size, origin, dataSaveType);
    nii.hdr.dime.dim = dim;
    
    
    % Get StudyName for SaveName
    pathsplit = strsplit(pathname,'\');
    [temp nameNum]=size(pathsplit);
    pathsplit=pathsplit(nameNum);
    pathsplit = strsplit(char(pathsplit),'.');
    StudyName=pathsplit(1);
    
    
    pathsplit=strsplit(pathname, strcat(StudyName,'.'));
    
    
    
    if savenametype==1
        SaveName=strcat(pathsplit(1),'_Scan',scanNum,'_RARE_',recoType,'.nii');
    elseif savenametype==2
        counter=sprintf('%d',RARE_Count);
        SaveName=strcat(pathsplit(1),'_RARE_RAREscan',counter,'_',recoType,'_', n,'.nii');
    end
    nii.img = nii.img.*(SLOPE1/RG);
    
    if savenametype~=3
        save_nii(nii,char(SaveName))
    end
    
    
    
end

end
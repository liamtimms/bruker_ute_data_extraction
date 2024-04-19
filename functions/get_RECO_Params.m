
function [RG, RecoType, Reco_size, voxel_size, DataType, SLOPE1, SLOPE2] = get_RECO_Params(pathname, scanNum)

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
numRecos=numRecos-2; % there are always 2 additional

% For each RECO: Extract RG and SLOPE, etc.
for j=1:1:numRecos
    recoNum=sprintf('%d',j);
    
    % Search RECO FILE...
    scanFileReco=strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\reco');
    [fidReco temp] = fopen(scanFileReco,'rt');
    if strcmp(temp,'')
        FileString=fscanf(fidReco, '%c');
        
        RecoType=get_Param(FileString, 'RECO_image_type=');
        
        RecoFOV =get_Param(FileString, 'RECO_fov='); %the field of view
        SplitSpot=strfind(RecoFOV, ')');
        RecoFOV=RecoFOV(SplitSpot+1:end);
        RecoLineSplit = strsplit(RecoFOV,' ');
        RECO_fov = [str2num(char(RecoLineSplit))]'; %convert string into 3x1 matrix
        RECO_fov = RECO_fov.*10; % convert FOV from cm to mm
        
        RecoSize =get_Param(FileString, 'RECO_size='); %the matrix size
        SplitSpot=strfind(RecoFOV, ')');
        RecoFOV=RecoSize(SplitSpot+1:end);
        RecoLineSplit = strsplit(RecoFOV,' ');
        RECO_size = [str2num(char(RecoLineSplit))]'; %convert string into 3x1 matrix
        
        voxel_size=RECO_fov./RECO_size;
    end
    
    % Search VISU_PARS FILE...
    scanFileVisuPars=strcat(pathname,'\',scanNum,'\pdata\',recoNum,'\visu_pars');
    [fidVisuPars temp] = fopen(scanFileVisuPars,'rt');
    
    if strcmp(temp,'')
        FileString=fscanf(fidReco, '%c');
        
        Slope=get_Param(FileString, 'VisuCoreDataSlope=');
        SplitSpot=strfind(Slope, ')');
        Slope=Slope(SplitSpot+1:end);
        SlopeLineSplit = strsplit(Slope,' ');
        SLOPE1 = str2num(char(slopelineSplit(1)));
        ss=size(slopelineSplit);
        if ss(1)>1
            SLOPE2=str2num(char(slopelineSplit(2))); % Handles recos with 2 slopes too
        end
        
        DataType=get_Param(FileString, 'VisuCoreWordType=');
        if strcmp('16BIT_SGN_INT',DataType)==1
            DataType = 'bit16';
        elseif strcmp('32BIT_SGN_INT',DataType)==1
            DataType = 'bit32';
        end 
    end
    
end
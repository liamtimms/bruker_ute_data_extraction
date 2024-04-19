% AUTOMATING PROCESSING : LIAM MRI EXTRACTION
% This code was written using by Liam Timms while at
% Northeastern to handle additional scan types, extract more
% information for an overview file, and work with PV6 data.
% An older PV5 code written by Codi A. Gharagouzloo and Praveen Kulkarni 
% while at Northeastern University was used an initial 
% reference implementation.

clear
clc
close all

% INITATE PROGRAM ->

% USER SELECTS PATIENT FILE
if  ~exist('pathname','var')
    [pathname] = uigetdir('*','select patient raw data folder');
end
dirData = dir(pathname);     % Get the data for the current directory
dirIndex = [dirData.isdir];  % Find the index for directories
s=size(dirData);

if s(1) <1
    msg = sprintf('...no scans found in selected folder ');
    error(msg);
    return;
end

FileSubject=strcat(pathname,'\subject');
[FileID, temp] = fopen(FileSubject,'rt');
FileString=fscanf(FileID, '%c');

ParamName='SUBJECT_name=';
SubjectName = get_Param(FileString, ParamName);
Less=strfind(SubjectName, '<');
More=strfind(SubjectName, '>');
SubjectName=SubjectName(Less(1)+1:More(1)-1);
SubjectName
ScanOverview.Subject(1,1)={SubjectName};

ParamName='SUBJECT_weight=';
SubjectWeight = get_Param(FileString, ParamName);
ScanOverview.Subject(1,2)={SubjectWeight};

ParamName='SUBJECT_sex_animal=';
SubjectSex = get_Param(FileString, ParamName);
ScanOverview.Subject(1,3)={SubjectSex};

ParamName='SUBJECT_position=';
SubjectPos = get_Param(FileString, ParamName);
ScanOverview.Subject(1,4)={SubjectPos};


% EXTRACT TIMES AND FOLDER NUMBERS
% We want to order and number the scans according to the time they were
% taken rather than the "folder numer" which varies according to bullshit
for i=1:s(1)
    i
    currName=dirData(i).name;
    [folderNumber, status]=str2num(currName);
    if status==1
        % SEARCH METHOD FILE FOR TIMES
        scanNum=currName;
        scanFileMethod=strcat(pathname,'\',scanNum,'\method');
        [FileID, temp] = fopen(scanFileMethod,'rt');
        FileString=fscanf(FileID, '%c');
        Folders(i,1)={currName};
        BrukerTime = get_BrukerTime(FileString);
        Dates(i,1)=BrukerTime;
    end
end

[~, order]=sort(Dates(:,1));
SortedFolders=Folders(order,:);
SortedFolders=SortedFolders(~cellfun('isempty',SortedFolders)); %remove empties

%%%%%% BEGIN EXTRACTION %%%%%%
% SAVENAME TYPE: 1 for scan number, 2 for progression, 3 to not save (just
% make ScanOverview
savenametype=1;
UTE_Count=0;
AFI_Count=0;
RARE_Count=0;
MSME_Count=0;
n=0;

s=size(SortedFolders);

for i=1:s(1)
    currName=SortedFolders{i,1};
    status=0;
    
    % SEARCH METHOD FILE...
    scanNum=currName;
    scanFileMethod=strcat(pathname,'\',scanNum,'\method');
    [FileID, temp] = fopen(scanFileMethod,'rt');
    FileString=fscanf(FileID, '%c');
    
    ScanOverview.FolderNumber(i,1)={currName};
    
    ParamName='Method=';
    MethodType = get_Param(FileString, ParamName);
    ScanOverview.ScanTypes(i,1)={MethodType};
    
    ParamName='PVM_EchoTime=';
    TimeToEcho = get_Param(FileString, ParamName);
    ScanOverview.TEandTR(i,1)={TimeToEcho};
    
    ParamName='PVM_RepetitionTime=';
    TimeToRepetition = get_Param(FileString, ParamName);
    ScanOverview.TEandTR(i,2)={TimeToRepetition};
    
    ParamName='ExcPulse1=';
    ExcPulseParams = get_Param(FileString, ParamName);
    ScanOverview.Pulse(i,2)={ExcPulseParams};
    Commas=strfind(ExcPulseParams, ',');
    FlipAngle=str2num(ExcPulseParams(Commas(2)+1:Commas(3)-1));
    ScanOverview.Pulse(i,1)={FlipAngle};
    
    ParamName='PVM_SpatResol=';
    SpatRes = get_Param(FileString, ParamName);
    Spaces=strfind(SpatRes, ' ');
    if numel(Spaces)>3
        Dim=str2num(SpatRes(Spaces(1)+1:Spaces(2)-1));
        if Dim==3
            X_Res=str2num(SpatRes(Spaces(2)+3:Spaces(3)-1));
            Y_Res=str2num(SpatRes(Spaces(3)+1:Spaces(4)-1));
            Z_Res=str2num(SpatRes(Spaces(4)+1:end));
            ScanOverview.Resolution(i,1)={X_Res};
            ScanOverview.Resolution(i,2)={Y_Res};
            ScanOverview.Resolution(i,3)={Z_Res};
        elseif Dim==2
            X_Res=str2num(SpatRes(Spaces(2)+3:Spaces(3)-1));
            Y_Res=str2num(SpatRes(Spaces(3)+1:Spaces(4)-1));
            ScanOverview.Resolution(i,1)={X_Res};
            ScanOverview.Resolution(i,2)={Y_Res};
        end
    end
    
    ParamName='PVM_SPackArrNSlices=';
    NumSlices = get_Param(FileString, ParamName);
    Spaces=strfind(NumSlices, ' ');
    NSlices=str2num(NumSlices(Spaces(2)+3:end));
    ScanOverview.NSlices(i,1)={NSlices};
    
    ParamName='PVM_NAverages=';
    NAve = get_Param(FileString, ParamName);
    ScanOverview.NAve(i,1)={NAve};
    
    ParamName='PVM_ScanTime=';
    ScanLength = get_Param(FileString, ParamName);
    ScanOverview.ScanLength(i,1)={ScanLength};
    
    ParamName='PVM_SPackArrSliceOrient=';
    SPackArrSliceOrient = get_Param(FileString, ParamName);
    Spaces=strfind(SPackArrSliceOrient, ' ');
    Direction=str2num(SPackArrSliceOrient(Spaces(1)+1:Spaces(2)-1));
    Slice_orientation=SPackArrSliceOrient(Spaces(2)+2:end);
    
    ParamName='PVM_SPackArrSliceDistance=';
    SliceDistStr = get_Param(FileString, ParamName);
    Spaces=strfind(SliceDistStr, ' ');
    SliceDistance=str2num(SliceDistStr(Spaces(2)+3:end));
    
    BrukerTime = get_BrukerTime(FileString);
    ScanOverview.DateAndTime(i,1)={BrukerTime};
    
    [RG, RecoType, voxel_size, DataType, SLOPE1, SLOPE2] = get_RECO_Params(pathname, scanNum);
    
    if strcmp(MethodType, '<Bruker:UTE3D>')
        UTE_Count=UTE_Count+1; % Used for savenametype2
        
        [UTENii, slopelineSplit]=ThreeDUTEextract_PV6_v2(pathname,scanNum,savenametype,UTE_Count,n,SubjectName);
        ScanOverview.UTE.Niis(i,1)={UTENii};
    elseif strcmp(MethodType,'<Bruker:RARE>')
        if Direction==1
            RARE_Count=RARE_Count+1; % Used for savenametype2
            RARENii=RAREextract_PV6_v2(pathname,scanNum,savenametype,RARE_Count,n,SubjectName,NSlices);
            ScanOverview.RARE.Niis(i,1)={RARENii};
        end
    elseif strcmp(MethodType, '<Bruker:MSME>')
        if Direction==1
            MSME_Count=MSME_Count+1; % Used for savenametype2
            MSMENii=MSMEextract_PV6_v2(pathname,scanNum,savenametype,MSME_Count,n,SubjectName,NSlices);
            ScanOverview.MSME.Niis(i,1)={MSMENii};
        end
    end
    
    if status==1
        fprintf('scan %d extracted', i);
        n=n+1;
    end
    
    clearvars UTENii AFINii MSMENii RARENii
    fclose all;
    
end

% Get StudyName for SaveName
pathsplit = strsplit(pathname,'\');
[temp nameNum]=size(pathsplit);
pathsplit=pathsplit(nameNum);
pathsplit = strsplit(char(pathsplit),'.');
StudyName=pathsplit(1);

% Save ScanOverview
SaveName=char(strcat(StudyName,'_ScanOverview'));
save(SaveName,'ScanOverview');

%Old version from Codi
%     while feof(fid) == 0
%         tline = fgetl(fid);
%         % CASE 1: UTE3D
%         scanMethodCheck=  regexpi(tline,'Method=UTE3D');
%         if ~isempty(scanMethodCheck)
%             pathname
%             scanNum
%             ThreeD_UTEextract(pathname,scanNum,savenametype,UTE_Count)
%             UTE_Count=UTE_Count+1; % Used for savenametype2
%         end
%         clearvars scanMethod;
%         % CASE 2: AFI
%         scanMethodCheck=  regexpi(tline,'Method=afi');
%         if ~isempty(scanMethodCheck)
%             pathname
%             scanNum
%             AFIextract(pathname,scanNum,savenametype,AFI_Count)
%             AFI_Count=AFI_Count+1; % Used for savenametype2
%         end
%         % CASE 3: T1
%         scanMethodCheck=  regexpi(tline,'Method=RAREVTR');
%         if ~isempty(scanMethodCheck)
%             pathname
%             scanNum
%             T1extract(pathname,scanNum,savenametype,T1_Count)
%             T1_Count=T1_Count+1; % Used for savenametype2
%         end
%         % IF ANATOMICAL: EXTRACT
%         % IF 3DUTE: EXTRACT, FIND RG AND SLOPE, SCALE, EXPORT ORIG AND SCALED
%         % -> ASSUME PDATA 1 IS INTENSITY, PDATA 2 IS COMPLEX (REAL AND IMAGINARY)
%
%         % IF T1: TAKE T1 STUFF LIKE WHAT PRAVEENS PROGRAM DOES
%         % extractTvalNC program from praveen
%         % IF T2: TAKE T2 STUFF LIKE WHAT PRAVEENS PROGRAM DOES
%     end

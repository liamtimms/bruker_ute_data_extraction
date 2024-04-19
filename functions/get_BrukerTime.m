function [BrukerTime] = get_BrukerTime(FileString)

%Extract the relevant part of the string
ParamName='OWNER';
ParamPos=strfind(FileString, ParamName);
if ~isempty(ParamPos)
    Grab=FileString(ParamPos(1):end);
    EndPos=strfind(Grab, '##');
    Grab(EndPos(1)-1:end)=[]; % the -1 deals with the hidden "return" character
    StartPos=strfind(Grab, '$$');
    Grab(1:StartPos(1))=[];
    
    Spaces=strfind(Grab, ' ');
    DateString=Grab(Spaces(1)+1:Spaces(2)-1);
    DateDash=strfind(DateString, '-');
    TimeString=Grab(Spaces(2)+1:Spaces(3)-1);
    TimeSep=strfind(TimeString, ':');
    TimeSep(3)=strfind(TimeString, '.');
    
    Y=str2num(DateString(1:DateDash(1)-1));
    M=str2num(DateString(DateDash(1)+1:DateDash(2)-1));
    D=str2num(DateString(DateDash(2)+1:end));
    
    H=str2num(TimeString(1:TimeSep(1)-1));
    MI=str2num(TimeString(TimeSep(1)+1:TimeSep(2)-1));
    S=str2num(TimeString(TimeSep(2)+1:TimeSep(3)-1));
    %     MS=str2num(TimeString(TimeSep(3)+1:end));
    %     BrukerTime = datetime(Y,M,D,H,MI,S,MS); new version of matlab can
    %     do ms but not this version
    
    BrukerTime = datetime(Y,M,D,H,MI,S);
else
    BrukerTime=NaN;
end

end
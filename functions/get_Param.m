function [Parameter] = get_Param(FileString, ParamName)

%Extract the relevant part of the string
ParamPos=strfind(FileString, ParamName);
if ~isempty(ParamPos)
    Grab=FileString(ParamPos(1):end);
    EndPos=strfind(Grab, '##');
    Grab(EndPos(1)-1:end)=[]; % the -1 deals with the hidden "return" character
    StartPos=strfind(Grab, '=');
    Grab(1:StartPos(1))=[];
    
    EndPos=strfind(Grab, '$$');
    if ~isempty(EndPos)
        Grab(EndPos(1)-1:end)=[];
    end
    Parameter=Grab;
    
    %test to see if this is a number and convert if so
    DigitTest= isstrprop(Parameter,'digit');
    Decimal = strfind(Parameter,'.');
    DigitTest(Decimal)=1;
    s1=size(nonzeros(DigitTest));
    s2=size(DigitTest);
    if s1(1)==s2(2)
        Parameter=str2num(Parameter);
    end  
else
    Parameter=NaN;
end

end
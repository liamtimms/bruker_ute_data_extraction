function [AFIimg] = get_AFIcalculation(IntensityImg_TR1,IntensityImg_TR2, TR1, TR2)

n=TR2/TR1;
s=size(IntensityImg_TR1);
if s~=size(IntensityImg_TR2)
    fprintf('Bro, wtf. The TR images are different sizes')
else
    rImg=IntensityImg_TR2./IntensityImg_TR1;
    for i=1:1:s(1)
        for j=1:1:s(2)
            for k=1:1:s(3)
                x=(n*rImg(i,j,k)-1)/(n-rImg(i,j,k));
                if abs(x)<1
                    AFIimg(i,j,k)=rad2deg(acos(x));
                else
                    AFIimg(i,j,k)=NaN;
                end
            end
        end
    end
end

end
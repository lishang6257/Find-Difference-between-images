% 同态滤波器
% ImageIn   -   需要进行滤波的灰度图像
% High      -   高频增益,需要大于1
% Low       -   低频增益,取值在0和1之间
% C         -   锐化系数
% Sigma     -   截止频率，越大图像越亮
% 输出为进行滤波之后的灰度图像
function [ImageOut] = HomoFilter(ImageIn, High, Low, C, Sigma)
    Img = double(ImageIn);      %   转换图像矩阵为双精度型,不会改变数据本身
    [Height, Width] = size(ImageIn);      % 返回的行数和列数
    CenterX = floor(Width / 2);     %   中心点坐标
    CenterY = floor(Height / 2);
    LogImg = log(Img + 1);      %   图像对数数据
    Log_FFT = fft2(LogImg);     %   傅里叶变换
    for Y = 1 : Height 
        for X = 1 : Width 
            Dist= (X - CenterX) * (X - CenterX) + (Y - CenterY) * (Y - CenterY);            %   点（X,Y）到频率平面原点的距离
            H(Y, X)=(High - Low) * (1 - exp(-C * (Dist / (2 * Sigma * Sigma)))) + Low;      %   同态滤波器函数
        end 
    end
    H = ifftshift(H);                   %   对H做反中心化                             
    Log_FFT = H.* Log_FFT;              %   滤波，矩阵点乘                                                
    Log_FFT = ifft2(Log_FFT);           %   反傅立叶变换 
    Out = exp(Log_FFT)-1;               %   取指数 

    % 指数处理ge = exp(g)-1;% 归一化到[0, L-1]
    Max = max(Out(:));
    Min = min(Out(:));
    Range = Max - Min;
    for Y = 1 : Height 
        for X = 1 : Width    
            ImageOut(Y, X) = uint8(255 * (Out(Y, X) - Min) / Range); 
        end
    end
end
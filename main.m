close all;
% clear all;

[a b] = imread('fly.gif');
if ~isempty(b)
    fly = ind2rgb(a,b);
end

img1_gray=rgb2gray(fly);
imwrite(img1_gray,'fly.jpg','quality',100);
fly_point = match('fly.jpg','fly.jpg');


imgOutVert = 'imgOutVert.png';
imgOutHori = 'imgOutHori.png';
matchings = 'matchings.txt';
keys1 = 'keys1.txt';
keys2 = 'keys2.txt';
file_img1 = 'fly.jpg';
file_img2 = 'region.jpg';
flag_resize = 0;




Start = 1;
End = 4;
for name = Start : End
    %%%%%%%%%% step 0 读取图片 %%%%%%%%%%%%%%%%%%%%
    imgName1 = ['test' num2str(name) '-1.jpg'];
    imgName2 = ['test' num2str(name) '-2.jpg'];
    img1 = imread(imgName1);
    img2 = imread(imgName2);
    img1_gray = rgb2gray(img1);
    img2_gray = rgb2gray(img2);

    moving = img2_gray;
    fixed = img1_gray;
    
%     figure(name+30)
%     figure, imshowpair(img1, img2, 'montage');
    
    %%%%%%%%%% step 1 初始配准（粗配准） %%%%%%%%%%%%%%%%%%%%
    [optimizer, metric] = imregconfig('multimodal');
%     % 有两种选择‘monomodal’, 'multimodal'两种，分别质量两幅图像是单一模态还是多模态，根据需要自己选择。
    
    %%%%%%%%%%% step 2  优化 %%%%%%%%%%%%%%
    % disp(optimizer);
    % disp(metric);
    
    %%%% 2.1 改变优化器的步长已达到对更加精细的变换。
    optimizer.InitialRadius = optimizer.InitialRadius/3.5;
    
    %%%% 2.2 （2.1）基础上改变最大迭代次数
    optimizer.MaximumIterations = 400;
    
    %%%% 3 初始条件
    tformSimilarity = imregtform(moving,fixed,'similarity',optimizer,metric);
    
%     %%%% 4 精对准
    tform = imregtform(moving,fixed,'affine',optimizer,metric,...
        'InitialTransformation',tformSimilarity);
    
    Rmoving = imref2d(size(moving));
    Rfixed = imref2d(size(fixed));
    [movingReg_rgb] = imwarp(img2,Rmoving,tform,'OutputView',Rfixed, 'SmoothEdges', true);
    movingRegX = rgb2gray(movingReg_rgb);

    
    % 转lab空间
    img_sub1_lab = rgb2lab(movingReg_rgb);
    img_sub2_lab = rgb2lab(img1);
    % 计算色差（色彩相似度）
    img_deltaE = sqrt( sum( (img_sub1_lab - img_sub2_lab).^2 , 3));

    
    threshold = 25;
    img_deltaE_bw = ones(size(img_deltaE));
    img_deltaE_bw(img_deltaE <= threshold) = 0;
    img_deltaE_bw(img_deltaE > threshold) = 1;
    img_deltaE_bw = logical(img_deltaE_bw);

   
    
    img_sub1 = imsubtract(fixed,movingRegX);
    img_sub2 = imsubtract(movingRegX,fixed);
    img_sub = img_sub1 + img_sub2;
    res1 = medfilt2(img_sub);
    if name == 1
        res = HomoFilter(res1, 0.05 ,0.1 , 500 , 200 );
    elseif name == 2
        res = HomoFilter(res1, 0.4 ,0.9 , 500 , 200 );
    elseif name == 3
        res = HomoFilter(res1, 0.7 ,0.4 , 500 , 200 );
    elseif name == 4
        res = HomoFilter(res1, 1 ,1, 500 , 200 );
    end
    t = graythresh(res);
    res = im2bw( real(res),t);
    

    res2 = imreconstruct(res,img_deltaE_bw);
    res = max(res2,res);
    
    res = imclose(res,strel('square',3));
    res = bwfill(res,'holes');
    res = imopen(res,strel('square',3));
    res = medfilt2(res);
    
    dis = 4;
    res(1:1+dis,:) = 0;
    res(:,1:1+dis) = 0;
    res(size(img1,1)-dis:size(img1,1),:) = 0;
    res(:,size(img1,2)-dis:size(img1,2)) = 0;
    
    
%     figure(name+10)
%     subplot(221),imshowpair(fixed, movingRegX);
%     subplot(222),imshow(img_sub1+img_sub2);
%     subplot(223),imshow(res);
%     subplot(224),imshow(img_deltaE_bw);
    
    figure(name+20)
    subplot(121)
    imshow(img1);
    hold on;
    region = regionprops(bwlabel(res),'Boundingbox');
    for i = 1 : size(region,1)
        rec = floor(region(i).BoundingBox);
        imwrite(  res( rec(2):rec(2)+rec(4) , rec(1):rec(1)+rec(3) ) ,'region.jpg','quality',100);
        p = match('region.jpg','fly.jpg');
        if p > 0
            rectangle('Position',region(i).BoundingBox,'Curvature',[0,0],'LineWidth',3 ,'LineStyle','-','EdgeColor','b');
        else
            rectangle('Position',region(i).BoundingBox,'Curvature',[0,0],'LineWidth',3 ,'LineStyle','-','EdgeColor','r');
        end
    end
    hold off
    subplot(122),imshow(img2);
    
end
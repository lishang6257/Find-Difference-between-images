% num = match(image1, image2)
%
% This function reads two images, finds their SIFT features, and
%   displays lines connecting the matched keypoints.  A match is accepted
%   only if its distance is less than distRatio times the distance to the
%   second closest match.
% It returns the number of matches displayed.
%
% Example: match('scene.pgm','book.pgm');

function [in,out,p] = matchDir(image1, image2)% image1 = 'v1.jpg';image2 = 'v2.jpg';
% Find SIFT keypoints for each image
[im1, des1, loc1] = sift(image1);
[im2, des2, loc2] = sift(image2);

% For efficiency in Matlab, it is cheaper to compute dot products between
%  unit vectors rather than Euclidean distances.  Note that the ratio of 
%  angles (acos of dot products of unit vectors) is a close approximation
%  to the ratio of Euclidean distances for small angles.
%
% distRatio: Only keep matches in which the ratio of vector angles from the
%   nearest to second nearest neighbor is less than distRatio.
distRatio = 0.6;   

% For each descriptor in the first image, select its match to second image.
des2t = des2';                          % Precompute matrix transpose

match = ones(1, max( size(loc1,1) , size(loc2,1) ) );
for i = 1 : size(des1,1)
   dotprods = des1(i,:) * des2t;        % Computes vector of dot products
   [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results

   % Check if nearest neighbor has angle less than distRatio times 2nd.
   if (vals(1) < distRatio * vals(2))
      match(i) = indx(1);
   else
      match(i) = 0;
   end
end

% Create a new image showing the two images side by side.
% Show a figure with lines joining the accepted matches.


% im3 = appendimages(im1,im2);
% figure('Position', [100 100 size(im3,2) size(im3,1)]);
% colormap('gray');
% imagesc(im3);
% hold on;
% cols1 = size(im1,2);


num = 0;
sum_a = 0;
sum_b = 0;
res = zeros(sum(match > 0),5);
for i = 1: size(des1,1)
  if (match(i) > 0)
      
%     line([loc1(i,2) loc2(match(i),2)+cols1], ...
%          [loc1(i,1) loc2(match(i),1)], 'Color', 'c');
     
    a = loc2(match(i),2) - loc1(i,2);
    b = loc2(match(i),1) - loc1(i,1);
    scale = loc2(match(i),3) - loc1(i,3);
    dir = loc2(match(i),4) - loc1(i,4);
    num = num +1;
    res(num,:) = [i a  b  scale dir];
  end
end

hold off;

% 奇异点剔除[正态分布]
in = zeros(num,2);
out = zeros(num,2);
res_std = std( res(:,2:5) );
res_mean = mean( res(:,2:5) );
res_new = res;
num = 0;
for i = 1 : size(res,1)
    if  match(i) > 0 ...
        && ( res(i,2) < (res_mean(1) + 3*res_std(1)) && res(i,2) > (res_mean(1) - 3*res_std(1)) ) ...
        && ( res(i,3) < (res_mean(2) + 3*res_std(2)) && res(i,3) > (res_mean(2) - 3*res_std(2)) ) ...
        && ( res(i,4) < (res_mean(3) + 3*res_std(3)) && res(i,4) > (res_mean(3) - 3*res_std(3)) ) ...
        && ( res(i,5) < (res_mean(4) + 3*res_std(4)) && res(i,5) > (res_mean(4) - 3*res_std(4)) ) %...
%         && abs(res(i,2)) < 30 && abs(res(i,3)) < 30 && abs(res(i,4)) < 0.2
    % 满足3-sigma法则
        num = num + 1;
        res_new(num,:) = res(i,:);
        in(num,:) = [loc1(i,2) loc1(i,1)];
        out(num,:) = [loc2(match(i),2) loc2(match(i),1)];
%         select = res(i,:)
        
    end
end
res_new = res_new(1 : num,:);

% for j = 1: size(res_new,1)
%     i = res_new(j,1);
%     if i < 1
%         continue;
%     end
% %     [loc1(i,2) loc2(match(i),2)+cols1]
% %          [loc1(i,1) loc2(match(i),1)]
%     line([loc1(i,2) loc2(match(i),2)+cols1], ...
%          [loc1(i,1) loc2(match(i),1)], 'Color', 'c');
% end
% hold off;
res_new_mean = mean(res_new(:,2:4));
p = res_new_mean(1:2);
fprintf('Found (%d , %d )\n', p);
% num = round(num);

% fprintf('Found %d matches.\n', num);





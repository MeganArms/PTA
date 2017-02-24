function [molPixelSubs,BW2,altBW] = RoughScan(obj,RawImage,k)

Option = obj.Option;
R = Option.spotR;
img = RawImage;

[M,N] = size(img);
img = imcomplement(img);
if strcmp(Option.illumination,'on')
    % Modified tophat filtering to remove uneven background
    img1 = imopen(img,strel('diamond',15));
    img2 = imgaussfilt(img1,14);
    img3 = img - img2;
    img4 = img3 + abs(min(img3(:)));
    img5 = img4/max(img4(:));
    img_2 = imadjust(img5);
else
    img_2 = imadjust(img);
end
if Option.exclude ~= 0
    x1 = Option.exclude(1,1);
    y1 = Option.exclude(1,2);
    x2 = Option.exclude(2,1);
    y2 = Option.exclude(2,2);
    img_2 = imadjust(img_2);
    img_2(x1:x2,y1:y2) = min(img_2(:));
end
if Option.include
    minval = Option.include(1,1);
    maxval = Option.include(1,2)-1;
    I = img_2(minval:maxval, minval:maxval);
    topad = (M - (maxval-minval+1))/2;
    img_2 = padarray(I,[topad topad],background,'both');
end
img2 = double(img_2);
% th2 = mean(img2(:)) + 2*std(img2(:));
% th = mean(img2(:)) + 2.5*std(img2(:));
[~,e] = histcounts(img2,100);
th = e(100); th2 = e(85);
obj.Frame(k).Threshold = th;
BW = imbinarize(img2,th);
altBW = imbinarize(img2,th2);

BWpad = padarray(BW,[1 1],0,'both');                                                                                                                                                                                                                                                                                                      
testmat1 = zeros(M+2,N+2); testmat1(1:end-2,1:end-2) = BW;
testmat2 = zeros(M+2,N+2); testmat2(3:end,3:end) = BW;
testmat3 = zeros(M+2,N+2); testmat3(3:end,1:end-2) = BW;
testmat4 = zeros(M+2,N+2); testmat4(1:end-2,3:end) = BW;

testmat5 = zeros(M+2,N+2); testmat5(3:end,2:end-1) = BW;
testmat6 = zeros(M+2,N+2); testmat6(1:end-2,2:end-1) = BW;
testmat7 = zeros(M+2,N+2); testmat7(2:end-1,3:end) = BW;
testmat8 = zeros(M+2,N+2); testmat8(2:end-1,1:end-2) = BW;

tallymat = BWpad + testmat1 + testmat2 + testmat3 + testmat4 + testmat5 ...                                                                                          
    + testmat6 + testmat7 + testmat8;

BW2 = logical(BW.*(tallymat(2:end-1,2:end-1) > 2)); % if >1 -> there's a particle

% BW2 = imclose(BW,strel('diamond',3));
CC = regionprops(BW2,'PixelIdxList');
% j = 1;
% for i = 1:CC.NumObjects
%     if length(CC.PixelIdxList{i}) >= 50
%         PixelIdxList{j} = CC.PixelIdxList{i};
%         j = j + 1;
%     end
% end

molPixelSubs = cell(1);
l = 1;
for k = 1:length(CC)
    % if numel(CC(k).PixelIdxList) < 50
    [i,j] = getcentroid(CC(k).PixelIdxList); %center of image
    if ge(i,R+1) && ge(M-R,i) && ge(j,R+1) && ge(N-R,j) %not in a corner/edge
        molPixelSubs{l} = [i,j];
        l = l+1;
        %         else
        %             CC(k) = []; % delete the objects that are on the corner or edge of image
    end
    % end
end

    function [c_row,c_col] = getcentroid(pixelIdxList)
        [rows,cols] = ind2sub([M,N],pixelIdxList);
        weight = img2(pixelIdxList);
        c_row = dot(rows,weight)/sum(weight);
        c_col = dot(cols,weight)/sum(weight);
        c_row = round(c_row);
        c_col = round(c_col);
    end

end

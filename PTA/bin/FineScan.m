function FineScan(obj,RawImage,k)
% FINESCAN(MOLPIXELIDX,RAWIMAGE) gets the detailed information for molecules
% identified in RAWIMAGE

Option = obj.Option;
R = Option.spotR;   % radius (pixel) of diffraction limited spot
% bg = min(RawImage(:));
[molPixelIdx,BW,altBW] = RoughScan(obj,RawImage,k);
th = obj.Frame(k).Threshold;
img = convert2double(RawImage);
NumMolecule = length(obj.Molecule);
nm = 1; altInd = [];

for k = 1:length(molPixelIdx)
    if isempty(molPixelIdx{k})
        break
    end
    s = molPixelIdx{k}(1); s1 = s - R; s2 = s + R;
    t = molPixelIdx{k}(2); t1 = t - R; t2 = t + R;
    subImage = img(s1:s2,t1:t2);
    BW_sub = BW(s1:s2,t1:t2);
    altBW_sub = altBW(s1:s2,t1:t2);
    CC_sub = regionprops(BW_sub,'Centroid','PixelIdxList','FilledArea','MajorAxisLength','MinorAxisLength','EulerNumber','Eccentricity'); %output struct 
    altCC = regionprops(altBW_sub,'PixelIdxList','EulerNumber','Centroid');
    
    % Get indices relative to entire image
%     [sg,tg] = meshgrid(s1:s2,t1:t2);
%     subs = reshape(cat(3,sg,tg),[numel(subImage),2]);
%     full_subim_inds = sub2ind(size(RawImage),subs(:,1),subs(:,2));
    
    % Deal with the above threshold pixels in the peripheral of subimage
    N = length(CC_sub); % number of objects in a subimage - should be 1
    if N > 1
        obj_idx = 2*R^2+2*R+1;  
        for l = 1:N
            pixIdxList = CC_sub(l).PixelIdxList;
            if ~ismember(obj_idx,pixIdxList)
                subImage(pixIdxList) = min(min(subImage));
            end
        end
    end
    % Perform centroid fitting
    edgeThreshold = th;
    x1 = ceil(R/2); x2 = x1 + R;
    edgeImage = subImage; edgeImage(x1:x2, x1:x2) = 0;
    subImage(edgeImage > edgeThreshold) = min(min(subImage));
    centroid = regionprops(true(size(subImage)),subImage,'WeightedCentroid');
    % Center the centroid in the center of the center pixel
    u = centroid.WeightedCentroid(2)-R-0.5;
    v = centroid.WeightedCentroid(1)-R-0.5;
    
    % Get the largest object in the subimage
    if N > 1
        centroid_idx = sub2ind(size(subImage),...
            round(centroid.WeightedCentroid(2)),round(centroid.WeightedCentroid(1)));
        obj_idx = [];
        % Find the most central object
        for l = 1:N
            if sum(ismember(CC_sub(l).PixelIdxList,centroid_idx)) == 1
                obj_idx = l;
                EN = CC_sub(obj_idx).EulerNumber; break
            end
        end
        % If no object is in the center, it's probably a giant ring around
        % the middle. Find the object closest to the center
        if isempty(obj_idx)
            EN = 0;
            x = round(centroid.WeightedCentroid(2));
            y = round(centroid.WeightedCentroid(1));
            d = squareform(pdist([x,y;cat(1,CC_sub.Centroid)]));
            inds = 1:length(d) - 1;
            obj_idx = inds(d(1,2:end) == min(d(1,2:end)));
        end
        % If this object is part of an alternate object, don't fit,continue
        [p,q] = ind2sub([2*R+1,2*R+1],CC_sub(obj_idx).PixelIdxList);
        subs = [p,q]; spot = [s,t];
        full_subs = spot - (R + 1 - subs);
        fullim_objinds = sub2ind(size(RawImage),full_subs(:,1),full_subs(:,2)); % Object indices in the full image
        if sum(ismember(fullim_objinds,altInd)) > 0
            continue
        end
        % Or find the largest object if none is on the centroid
%         if isempty(obj_idx)
%             lengths = zeros(1,N); idxs = 1:N;
%             for l = 1:N
%                 lengths(l) = length(CC_sub(l).PixelIdxList);
%             end
%             obj_idx = idxs(lengths == max(lengths));
%         end
    else
        obj_idx = 1;
        % If this object is part of an alternate object, don't fit,
        % continue
        [p,q] = ind2sub([2*R+1,2*R+1],CC_sub(obj_idx).PixelIdxList);
        subs = [p,q]; spot = [s,t];
        full_subs = spot - (R + 1 - subs);
        fullim_objinds = sub2ind(size(RawImage),full_subs(:,1),full_subs(:,2)); % Object indices in the full image
        if sum(ismember(fullim_objinds,altInd)) > 0
            continue
        end
        EN = CC_sub(obj_idx).EulerNumber;
    end
    % Find the object in the alt image that's closest to the main object in
    % the original sub image if the EulerNumber has changed
    if sum(cat(1,altCC.EulerNumber) < 1) > 0 && sum(cat(1,altCC.EulerNumber) < 1) > sum(cat(1,CC_sub.EulerNumber) < 1)
        centers = cat(1,altCC.Centroid);
        distances = squareform(pdist(cat(1,CC_sub(obj_idx).Centroid,centers)));
        distfromcenter = distances(1,2:end); label_inds = 1:length(centers);
        objofinterest = min(label_inds(distfromcenter == min(distfromcenter))); % Sometimes they're equidistant
        obj.Molecule(NumMolecule+nm).altEN = altCC(objofinterest).EulerNumber;
        % Add indices of this subimage's object to the alternate indices list
        % isfull_subim_inds(fullim_objinds == subim_inds)
        [p,q] = ind2sub([2*R+1,2*R+1],altCC(objofinterest).PixelIdxList);
        subs = [p,q]; spot = [s,t];
        full_subs = spot - (R + 1 - subs);
        fullim_objinds = sub2ind(size(RawImage),full_subs(:,1),full_subs(:,2));
        altInd = cat(1,altInd,fullim_objinds);
    end
    Maj = CC_sub(obj_idx).MajorAxisLength;
    Min = CC_sub(obj_idx).MinorAxisLength;
    Elongation = (Maj - Min)/Maj; % closer to 1, closer to a line. closer to 0, closer to a circle
    obj.Molecule(NumMolecule+nm).Area = CC_sub(obj_idx).FilledArea;
    obj.Molecule(NumMolecule+nm).Elongation = Elongation;
    obj.Molecule(NumMolecule+nm).Eccentricity = CC_sub(obj_idx).Eccentricity;
    obj.Molecule(NumMolecule+nm).EulerNumber = EN;
    obj.Molecule(NumMolecule+nm).coordinate = [s t];
    % Store the centroid in nm
    obj.Molecule(NumMolecule+nm).centroid = [u,v]*obj.Option.pixelSize;
    nm = nm + 1;
end

end
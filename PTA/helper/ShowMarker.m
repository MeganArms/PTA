  function ShowMarker(obj,index)
% Show where the identified single molecules are and number them according
% to which trajectory that they're a part of
% Input
%   IMAGE - image file name
%   INDEX - the index of the image
I = obj.filename;
Option = obj.Option;
Molecule = obj.Molecule;
Frame = obj.Frame;

Result = struct;
if isfield(Molecule,'connectedResult')
    for i = 1:length(Molecule)
        if ~isempty(Molecule(i).connectedResult)
            Result(i).trajectory = Molecule(i).connectedResult;
        else
            continue
        end
    end
else
    Result = obj.Result;
end
[~,~,ext] = fileparts(I);
if strcmp(ext, '.tiff') || strcmp(ext,'.tif')
    img = imread(I,index);
elseif strcmp(ext,'.nd2')
    data = bfopen(I);
    img = data{1}{index,1};
end
mol_ind = Frame(index).MoleculeIndex; % Vector of the molecule indices in this frame
N = length(mol_ind);
pts = zeros(N,2);
for i = 1:N
    pts(i,1) = Molecule(mol_ind(i)).coordinate(1);
    pts(i,2) = Molecule(mol_ind(i)).coordinate(2);
end
pts = int32(pts);


% Search each trajectory individually for the molecule index of interest
currentTraj = zeros(N,1);
for j = 1:N 
    for i = 1:length(Result)
        findValue = Result(i).trajectory == mol_ind(j);
        if sum(findValue) == 1
            currentTraj(j,1) = i;
            break;
%         else
%             currentTraj(j,1) = 0;
        end
        clear findValue
    end
end
          
%[M,~] = size(img);
img = convert2double(img);
if strcmp(Option.illumination,'on')
    img1 = imopen(img,strel('diamond',7));
    img2 = imgaussfilt(img1,14);
    img3 = img - img2;
    img4 = img3 + abs(min(img3(:)));
    img5 = img4/max(img4(:));
    img_1 = imadjust(img5);
%     % High pass filtering to remove uneven background
%     mid = floor(M/2)+1;
%     Img = fft2(img);
%     Img1 = fftshift(Img);
%     Img2 = Img1;
%     Img2(mid-3:mid+3,mid-3:mid+3) = min(min(Img1));
%     Img2(257,257) = Img1(257,257);
%     img1 = ifft2(ifftshift(Img2));
%     img12 = abs(img1);
%     img13 = img12-min(min(img12));
%     img14 = img13/max(max(img13));
%     % Mulitply pixels by the sum of their 8-connected neighbors to increase
%     % intensities of particles
%     img_1 = colfilt(img14,[3 3],'sliding',@colsp);
else
    img_1 = imadjust(img);
end

figure,imshow(imcomplement(img_1))
hold on, plot(pts(:,2),pts(:,1),'go'),hold off
for i = 1:length(pts)
    if currentTraj(i) ~= 0
        text((double(pts(i,2))+4),double(pts(i,1)),sprintf('%d %s, %s',...
            currentTraj(i),Result(currentTraj(i)).Motility_flag,...
            Result(currentTraj(i)).Aggregate_flag),'FontSize',14,'Color','g');
    end
end

title(['Frame number ',num2str(index)],'Fontsize',14);

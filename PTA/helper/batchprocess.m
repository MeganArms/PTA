function batchprocess
% load('/Users/MeganArmstrong 1/Documents/Hess Lab/BSA Project/Corina Data/20140514/80C/out/2014-05-14_BSA_flowcell_pluronic_80C_exp200ms_l10_EM200.mat');
% clearvars('-except','allDisplacements','allCounts','allLogCounts','allCoordinates','allSizes','allTraj');
clear
folderIn = uigetdir('/Users/MeganArmstrong 1/Documents/MATLAB');
folderOut = [folderIn, '/out']; %filtered images + filtered
dirListing = dir(folderIn); %struct - nb of files+names
numFiles = length(dirListing);
exptime = 0.2; % seconds
movietime = 60; % seconds
frames = movietime/exptime;

allDisplacements = cell(frames,3); % longest possible movie frames for preallocation
allCounts = zeros(1,frames); %histogram data, adds vertically
allLogCounts = zeros(1,frames);
allCoordinates = []; %all the individual trajectories 
allSizes = []; %gaussian or disk 
allTraj = []; %all molecular indices

initialvars = who; %all variables in the workspace

for i = 4:numFiles
    if dirListing(i).isdir %is directory -> next file
        break
    end
    clear obj
    % Initialize object
    obj = FSMIA([folderIn, '/', dirListing(i).name]);
    
    % Get threshold
    img = imread([folderIn, '/', dirListing(i).name],'Index',1);
    %img = I(513:1536,513:1536);
    obj.Option.illumination = 'off';
    if strcmp(obj.Option.illumination,'on')
        % High pass filtering to remove uneven background
        [M,~] = size(img); % tilda to say don't output
        mid = floor(M/2)+1; % midpoint
        Img = fft2(img);
        Img1 = fftshift(Img); %low frequency in center
        Img2 = Img1;
        if M == 2048
            Img2(mid-1:mid+1,:) = min(Img1(:)); % For 2048x2048 images, remove the center 3 lines
            Img2(:,mid-1:mid+1) = min(Img1(:));
        else
            Img2(mid-9:mid+9,mid) = min(min(Img1)); % 18x18 cross around the midpoint
            Img2(mid,mid-9:mid+9) = min(min(Img1));
        end
        Img2(mid,mid) = Img1(mid,mid); 
        img1 = ifft2(ifftshift(Img2)); % shifting back
        img12 = abs(img1); % magnitude
        img13 = img12-min(min(img12));
        img14 = img13/max(max(img13)); % btw 0 and 1
        % Mulitply pixels by the sum of their 8-connected neighbors to increase
        % intensities of particles
        outImage = imadjust(colfilt(img14,[3 3],'sliding',@colsp));
        %outImage = imadjust(img14);
    else
        outImage = imadjust(img);
    end
    imwrite(uint16(outImage),[folderOut,'/',dirListing(i).name]);
    th = .9999; bg = 0.9;
    % figure(1),imshow(outImage), th = input('Set threshold: '); 
    % bg = input('Set background: ');
    % close(figure(1))
    
    % Set options
    obj.Option.threshold = th;
    obj.Option.spotR = 10; %radius in pixels
    obj.Option.pixelSize = 70;
    obj.Option.include = 0;
    obj.Option.exclude = 0;
    obj.Option.connectDistance = 5*70; %
    obj.Option.ds = 1;
    obj.Option.fitting = 'slow'; %gaussian or dist
    obj.Option.isolation = 'fast'; %always
    obj.Option.bg = bg;
    obj.Option.wavelength = 647;
    obj.Option.na = 1.49;
    
    % Begin stack analysis
    analyzestack(obj, obj.filename);
    % Save output
    save([folderOut,'/',dirListing(i).name(1:end-3),'mat']);
    
    % Perform analysis
    createTrajectories(obj);
    longTraj = connectShortTraj(obj,exptime);
    coords = getCoordinates(obj,'yes');
    sizes = particleSize(obj);
    [~,Displacement,~] = findSteps(coords,1);
    logcount = logResidenceTimeStat(coords,'ExposureTime',exptime);
    logcount = padarray(logcount,[0 frames-length(logcount)],0,'post');
    count = ResidenceTimeStat(coords,'ExposureTime',exptime);
    count = padarray(count,[0 frames-length(count)],0,'post');
    % [msd,D] = Dcoeff(Displacement,exptime);
    % [logT,logSF] = logSurvivalFunction(298,exptime,logcount);%not used
    % [T,SF] = survivalFunction(ceil(movietime/exptime),exptime,count);
    % dSF = diffSurvival(Displacement,exptime,1,1,100);
    
    % Save output
    save([folderOut,'/',dirListing(i).name(1:end-3),'mat']);
    
    % Collect results
    allTraj = [allTraj, longTraj];
    allCoordinates = [allCoordinates; coords];
    allSizes = [allSizes; sizes];
    for k = 1:length(Displacement)
        allDisplacements{k} = [allDisplacements{k}; Displacement{k}];
    end
    allCounts = sum([allCounts; count],1);
    allLogCounts = sum([allLogCounts; logcount],1);
    clearvars('-except',initialvars{:},'i','initialvars')
    fprintf('%d videos analyzed!\n',i-3);
end
save('CollectedResults test.mat','allTraj','allCoordinates','allSizes','allDisplacements','allCounts','allLogCounts');

end


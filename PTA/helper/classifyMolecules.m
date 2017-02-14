function classifyMolecules(obj,traj,molInds)

% Threshold test for broken circles
% mInitial = obj.Result(trajNum).trajectory(1);
% centroid = obj.Molecule(mInitial).centroid;
% i = centroid(1); j = centroid(2);
% R = obj.Option.spotR;
% i1 = i - R; i2 = i + R; j1 = j - R; j2 = j + R;
% subIm = rawImage(i1:i2,j1:j2);
% subim = imcomplement(subim); 
% subim2 = imadjust(double(subim));
% th = mean(subim2(:)) + 2*std(subim2(:));
% bw = imbinarize(subim2,th);


trajlength = size(traj,1); samat = zeros(1,trajlength);
inds = 1:trajlength;
for i = 1:trajlength
    mIndex = molInds(i);
    % If EN never >= 1, bulk diffusion only, otherwise Lévy or Surface.
    % Only test for the "never" condition on the first iteration
    if sum(traj(:,1) < 1) == trajlength
        obj.Molecule(mIndex).Motility_flag = 'Bulk';
        obj.Molecule(mIndex).Aggregate_flag = 'Unknown';
    % elseif traj(i,1) < 1 && sum(traj(:,1)) >= 1
    elseif sum(diff(inds(traj(:,1) < 1)) > 1) >= 1 && traj(i,1) < 1 || sum(diff(inds(traj(:,1) >= 1)) > 1) >=1 && traj(i,1) < 1 || sum(traj(:,1) < 1) == 1 && traj(i,1) < 1
        obj.Molecule(mIndex).Motility_flag = 'Levy'; 
        obj.Molecule(mIndex).Aggregate_flag = [];
    elseif sum(diff(inds(traj(:,1) < 1)) > 1) == 0 && traj(1,1) < 1 && traj(i,1) < 1
        obj.Molecule(mIndex).Motility_flag = 'Landing';
        obj.Molecule(mIndex).Aggregate_flag = [];
    elseif sum(diff(inds(traj(:,1) < 1)) > 1) == 0 && traj(end,1) < 1 && traj(i,1) < 1
        obj.Molecule(mIndex).Motility_flag = 'Desorbed';
        obj.Molecule(mIndex).Aggregate_flag = [];
    elseif sum(diff(inds(traj(:,1) < 1)) > 1) >= 1 && traj(i,1) >= 1 && traj(i,2) <= 25 || sum(diff(inds(traj(:,1) >= 1)) > 1) >=1 && traj(i,1) >= 1 && traj(i,2) <= 25
        obj.Molecule(mIndex).Motility_flag = 'Surface'; 
        obj.Molecule(mIndex).Aggregate_flag = 'Single';
        samat(i) = false;
    elseif sum(diff(inds(traj(:,1) < 1)) > 1) >= 1 && traj(i,1) >= 1 && traj(i,2) > 25 || sum(diff(inds(traj(:,1) >= 1)) > 1) >=1 && traj(i,1) >= 1 && traj(i,2) > 25
        obj.Molecule(mIndex).Motility_flag = 'Surface'; 
        obj.Molecule(mIndex).Aggregate_flag = 'Aggregate';
        samat(i) = true;
    % Desorbed -> adsorbed
    elseif sum(diff(inds(traj(:,1) < 1)) > 1) == 0 && traj(1,1) < 1 && traj(i,1) >= 1 && traj(i,2) <= 25
        obj.Molecule(mIndex).Motility_flag = 'Adsorbed';
        obj.Molecule(mIndex).Aggregate_flag = 'Single';
        samat(i) = false;
    elseif sum(diff(inds(traj(:,1) < 1)) > 1) == 0 && traj(1,1) < 1 && traj(i,1) >= 1 && traj(i,2) > 25
        obj.Molecule(mIndex).Motility_flag = 'Adsorbed';
        obj.Molecule(mIndex).Aggregate_flag = 'Aggregate';
        samat(i) = true;
    % Adsorbed -> desorbed
    elseif sum(diff(inds(traj(:,1) < 1)) > 1) == 0 && traj(end,1) < 1 && traj(i,1) >= 1 && traj(i,2) <= 25
        obj.Molecule(mIndex).Motility_flag = 'Adsorbed';
        obj.Molecule(mIndex).Aggregate_flag = 'Single';
        samat(i) = false;
    elseif sum(diff(inds(traj(:,1) < 1)) > 1) == 0 && traj(end,1) < 1 && traj(i,1) >= 1 && traj(i,2) > 25
        obj.Molecule(mIndex).Motility_flag = 'Adsorbed';
        obj.Molecule(mIndex).Aggregate_flag = 'Aggregate';
        samat(i) = true;
    elseif sum(traj(:,1) < 1) == 0 && traj(i,2) <= 25
        obj.Molecule(mIndex).Motility_flag = 'Surface';
        obj.Molecule(mIndex).Aggregate_flag = 'Aggregate';
    elseif sum(traj(:,1) < 1) == 0 && traj(i,2) > 25
        obj.Molecule(mIndex).Motility_flag = 'Surface';
        obj.Molecule(mIndex).Aggregate_flag = 'Aggregate';
        
%     elseif traj(i,1) >= 1 && traj(i,2) <= 25
%         obj.Molecule(mIndex).Motility_flag = 'Surface';
%         obj.Molecule(mIndex).Aggregate_flag = 'Single';
%         samat(i) = false;
%     elseif traj(i,1) >= 1 && traj(i,2) > 25
%         obj.Molecule(mIndex).Motility_flag = 'Surface';
%         obj.Molecule(mIndex).Aggregate_flag = 'Aggregate';
%         samat(i) = true;
    end
end

% Assign aggregate status for levy flights
if sum(traj(:,1) < 1) > 0 && sum(traj(:,1) < 1) < trajlength
% if sum(traj(:,1)) >= 1
    for i = 1:trajlength
        mIndex = molInds(i);
        if isempty(obj.Molecule(mIndex).Aggregate_flag) && sum(samat) == length(samat)
            obj.Molecule(mIndex).Aggregate_flag = 'Aggregate';
        elseif isempty(obj.Molecule(mIndex).Aggregate_flag) && sum(samat) == 0
            obj.Molecule(mIndex).Aggregate_flag = 'Single';
        elseif sum(samat) > 0 && sum(samat) < length(samat)
            obj.Molecule(mIndex).Aggregate_flag = 'Aggregate';
        end
    end
end
    
    
    
    
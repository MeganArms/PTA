function classifyTrajectories(obj)

Molecule = obj.Molecule;
Result = obj.Result;
Levy = struct('trajectory',{}); Aggregate = struct('trajectory',{}); Unknown = struct('trajectory',{});
Surface = struct('trajectory',{}); Single = struct('trajectory',{}); Bulk = struct('trajectory',{});
Adsorption = struct('trajectory',{}); Desorption = struct('trajectory',{});
k = 1; l = 1; m = 1; n = 1; o = 1; p = 1; q = 1; r = 1;
% LevyTraj = 0; SurfaceTraj = 0;
for i = 1:length(Result)
    molInds = Result(i).trajectory;
    inds = 1:length(molInds);
    % Get descriptors in trajectory form
    traj = zeros(length(molInds),3);
    for j = 1:length(molInds)
        mIndex = molInds(j);
        if ~isempty(Molecule(mIndex).altEN)
            traj(j,1) = Molecule(mIndex).altEN;
        else
            traj(j,1) = Molecule(mIndex).EulerNumber;
        end
        traj(j,2) = Molecule(mIndex).Area;
        traj(j,3) = Molecule(mIndex).Eccentricity;
    end
    % Classify molecules
    if size(traj,1) > 2
        classifyMolecules(obj,traj,molInds);

        % Classify current trajectory
        % Euler number < 1 indicates a hole in the object, which indicates bulk 
        % diffusion
        if sum(traj(:,1) < 1) == size(traj,1)
            obj.Result(i).Motility_flag = 'Bulk';
            Bulk(k).trajectory = molInds; k = k + 1;
            obj.Result(i).Aggregate_flag = 'Unknown';
            Unknown(l).trajectory = molInds; l = l + 1;
        % Euler number < 1 indicates a hole in the object, if the whole
        % trajectory does NOT have EN < 1, then at some point it is on the
        % surface, making it a Levy flight.
        elseif sum(diff(inds(traj(:,1) < 1)) > 1) >= 1 || sum(diff(inds(traj(:,1) >= 1)) > 1) >=1
            obj.Result(i).Motility_flag = 'Levy';
            Levy(m).trajectory = molInds; m = m + 1;
        elseif sum(diff(inds(traj(:,1) < 1)) > 1) == 0 && traj(1,1) < 1  
            obj.Result(i).Motility_flag = 'Adsorption';
            Adsorption(q).trajectory = molInds; q = q + 1;
        elseif sum(diff(inds(traj(:,1) < 1)) > 1) == 0 && traj(end,1) < 1
            obj.Result(i).Motility_flag = 'Desorption';
            Desorption(r).trajectory = molInds; r = r + 1;
    %     elseif sum(traj(:,1)) >= 1 && sum(traj(:,1) < 1) >= 1 %If at some point it is on the surface and at some point it is off the surface
    %         obj.Result(i).Motility_flag = 'Levy';
    %         Levy(m).trajectory = molInds; m = m + 1;
        % If the Euler is never < 1, the molecule remains on the surface the
        % entire time
        elseif sum(traj(:,1) < 1) == 0 % If it is never off the surface 
            obj.Result(i).Motility_flag = 'Surface';
            Surface(n).trajectory = molInds; n = n + 1;
        end
        % Beads with an area <= 25 are single beads, larger are aggregates. If
        % there is a mix within one trajectory, assume aggregate.
        if sum(traj(:,2) > 25) >= 1 && sum(traj(:,2) <= 25) == 0 && ~strcmp(obj.Result(i).Motility_flag,'Bulk')
            obj.Result(i).Aggregate_flag = 'Aggregate';
            Aggregate(o).trajectory = molInds; o = o + 1;
        elseif sum(traj(:,2) <= 25) >= 1 && sum(traj(:,2) > 25) == 0 && ~strcmp(obj.Result(i).Motility_flag,'Bulk')
            obj.Result(i).Aggregate_flag = 'Single';
            Single(p).trajectory = molInds; p = p + 1;
        elseif sum(traj(:,2) <= 25) >= 1 && sum(traj(:,2) > 25) >= 1 && ~strcmp(obj.Result(i).Motility_flag,'Bulk')
            obj.Result(i).Aggregate_flag = 'Aggregate';
            Aggregate(o).trajectory = molInds; o = o + 1;
        end
    end
%     numHoles = sum(traj(:,1) < 1);
%     if numHoles >= 1
%         holes_at = molInds(traj(:,1) < 1);
%         if numHoles > 1
%             for j = 1:numHoles
%                 obj.Molecule(holes_at(j)).Motility_flag = 'Levy';
%             end
%             obj.Result(i).Motility_flag = 'Levy';
%             if molInds(end) == holes_at(end)
%                 Levy(k).trajectory = molInds(1:end-1);
%                 k = k + 1;
%                 LevyTraj = i;
%             else
%                 Levy(k).trajectory = molInds;
%                 k = k + 1;
%                 LevyTraj = i;
%             end
%         elseif numHoles == 1 && molInds(end) ~= holes_at
%             obj.Molecule(holes_at).Motility_flag = 'Levy';
%             obj.Result(i).Motility_flag = 'Levy';
%             Levy(k).trajectory = molInds;
%             k = k + 1;
%             LevyTraj = i;
%         elseif numHoles == 1 && molInds(end) == holes_at
%             obj.Molecule(holes_at).Motility_flag = 'Surface';
%             obj.Result(i).Motility_flag = 'Surface';
%             Surface(m).trajectory = molInds;
%             m = m + 1;
%             SurfaceTraj = i;
%         end
%     end
%     
%     % Area > mean + std of area of trajectory indicates out of focus bulk
%     % diffusion
% %     mean(traj(:,2))
% %     std(traj(:,2))
% %     threshold = mean(traj(:,2)) + std(traj(:,2))
% %     numOOF = sum(traj(:,2) > threshold);
% %     if numOOF >= 1
% %         flights_at = molInds(traj(:,2) > threshold);
% %         if numOOF > 1
% %             for j = 1:numOOF
% %                 obj.Molecule(flights_at(j)).Motility_flag = 'Levy';
% %             end
% %             obj.Result(i).Motility_flag = 'Levy';
% %             if molInds(end) == flights_at(end) && LevyTraj ~= i
% %                 Levy(k).trajectory = molInds(1:end-1);
% %                 k = k + 1;
% %                 LevyTraj = i;
% %             elseif LevyTraj ~= i
% %                 Levy(k).trajectory = molInds;
% %                 k = k + 1;
% %                 LevyTraj = i;
% %             end
% %         elseif numOOF == 1 && molInds(end) ~= flights_at
% %             obj.Molecule(flights_at).Motility_flag = 'Levy';
% %             obj.Result(i).Motility_flag = 'Levy';
% %             if LevyTraj ~= i
% %                 Levy(k).trajectory = molInds;
% %                 k = k + 1;
% %                 LevyTraj = i;
% %             end
% %         elseif numOOF == 1 && molInds(end) == flights_at
% %             obj.Molecule(flights_at).Motility_flag = 'Surface';
% %             obj.Result(i).Motility_flag = 'Surface';
% %             if SurfaceTraj ~= i
% %                 Surface(m).trajectory = molInds;
% %                 m = m + 1;
% %                 SurfaceTraj = i;
% %             end
% %         end
% %     end
%     
%     % If was set as neither Surface or Lévy, it's a surface trajectory
%     if SurfaceTraj ~= i && LevyTraj ~= i
%         Surface(m).trajectory = molInds;
%         m = m + 1;
%         SurfaceTraj = i;
%         for j = 1:length(molInds)
%             obj.Molecule(molInds(j)).Motility_flag = 'Surface';
%         end
%         obj.Result(i).Motility_flag = 'Surface';
%     end
%     
%     % Eccentricity > 0.8 or Area > 25 indicates aggregate if on surface
%     if strcmp(obj.Result(i).Motility_flag,'Surface') && sum(traj(:,3) > 0.8) > 1 || sum(traj(:,2) > 25) > 1
%         obj.Molecule(molInds(1)).Aggregate_flag = 'Aggregate';
%         obj.Result(i).Aggregate_flag = 'Aggregate';
%         Aggregate(l).trajectory = molInds;
%     elseif ~strcmp(obj.Result(i).Motility_flag,'Surface') && sum(traj(:,2) > 25) > 1
%         obj.Molecule(molInds(1)).Aggregate_flag = 'Unknown';
%         
% %         l = l + 1;
% %     elseif sum(traj(:,3) > 0.8) > 1 && sum(traj(:,2) > 25) <= 1
% %         obj.Molecule(molInds(1)).Aggregate_flag = 'Wonky particle';
% %         obj.Result(i).Aggregate_flag = 'Wonky particle';
%     else        
%         obj.Molecule(molInds(1)).Aggregate_flag = 'Single';
%         obj.Result(i).Aggregate_flag = 'Single';
%         Single(n).trajectory = molInds;
%         n = n + 1;
%     end
end
obj.Surface = Surface;
obj.Levy = Levy;
obj.Bulk = Bulk;
obj.Adsorption = Adsorption;
obj.Desorption = Desorption;
obj.Aggregate = Aggregate;
obj.Single = Single;
obj.Unknown = Unknown;
end
    
    
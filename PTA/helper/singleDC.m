function [D,D_er,eps,eps_er] = singleDC(data,exptime)
% individual diffusion coefficients
D = zeros(size(data,1),1); D_er = D; eps = D; eps_er = D;
for i = 1:size(data,1)
    if length(data{i,1}) >= 4
        [~,dis,~] = findSteps(data(i,1),1,0.1);
        [D(i),D_er(i),eps(i),eps_er(i),~] = msdcell(dis,exptime,'off');
    end
end
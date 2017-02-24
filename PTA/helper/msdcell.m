function [D,D_er,rho,rho_er,f] = msdcell(displacements,exptime,plotting)

SD = cellfun(@(x) x.^2,displacements,'UniformOutput',0);
MSD = cellfun(@mean,SD);
sd = cellfun(@std,displacements);
SEM = sd./sqrt(cellfun(@length,displacements));

maxstep = length(MSD);
if maxstep > 4
    maxstep = 4;
end

timevec = exptime:exptime:exptime*length(MSD);


fitopt = fitoptions('poly1');
fitopt.Weights = SEM(1:maxstep);
f = fit(timevec(1:maxstep)',MSD(1:maxstep),'poly1',fitopt);
ci = confint(f);
D = f.p1/4; D_er = (ci(3) - ci(1))/2;
rho = MSD(1); rho_er = SEM(1);
if strcmp(plotting,'on')
    figure,errorbar(timevec,MSD,SEM,'ko');
    hold on, plot(f), hold off
    xlabel('Time (s)'), ylabel('<\Delta x^2> (µm^2)');
end
end
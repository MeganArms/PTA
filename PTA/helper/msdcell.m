function [D,f] = msdcell(displacements,exp)

SD = cellfun(@(x) x.^2,displacements,'UniformOutput',0);
MSD = cellfun(@mean,SD);
sd = cellfun(@std,displacements);
SEM = sd./sqrt(cellfun(@length,displacements));

timevec = exp:exp:exp*length(MSD);

figure,errorbar(timevec,MSD,SEM,'ko');

fitopt = fitoptions('poly1');
fitopt.Weights = SEM;
f = fit(timevec',MSD,'poly1',fitopt);
D = f.p1;
hold on, plot(f), hold off
xlabel('Time (s)'), ylabel('<\Delta x^2> (µm^2)');

end
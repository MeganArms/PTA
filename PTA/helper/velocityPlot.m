function g = velocityPlot(data,exptime,L)

g = zeros(round(L));
r = zeros(size(data,1),1); ini = [r,r]; x = r; y = r; theta = r; m = r;
for i = 1:size(data,1)
    t = size(data{i,1},1)*exptime;
    ini(i,:) = data{i,1}(1,:);
    fin = data{i,1}(end,:);
    r(i) = pdist([ini(i,:);fin]);
    x(i) = (ini(i,2) - fin(2))/t;
    y(i) = (ini(i,1) - fin(1))/t;
    m(i) = sqrt(x(i)^2 + y(i)^2);
    theta(i) = atan360(x(i),y(i));
    g(round(ini(i,1)),round(ini(i,2))) = theta(i);
end

meanV = [mean(x), mean(y)]; mag = sqrt(meanV(1)^2 + meanV(2)^2);
meanTheta = atan360(meanV(1),meanV(2));
sdV = [std(x),std(y)];
sdTheta = atan360(sdV(1),sdV(2));
meanTheta = round(meanTheta*180/pi); sdTheta = round(sdTheta*180/pi);

% Visualize
figure,quiver(ini(:,1),ini(:,2),x,y);
xlabel('X (µm)','FontSize',14), ylabel('Y (µm)','FontSize',14)
title(['Avg. drift vel. magnitude: ',num2str(mag),' µm/s'],'FontSize',14);
h = gca; h.XLim = [0 L]; h.YLim = [0 L];
% figure,histogram(theta,0:pi/12:2*pi,'FaceColor',[0.5 0.5 0.5])
figure,histogram(theta*180/pi,0:180/12:360,'FaceColor',[0.5 0.5 0.5])
xlabel('Drift direction, \theta (º)','FontSize',14)
ylabel('f(\theta)','FontSize',14)
title(['Avg. drift vel. direction: ',num2str(meanTheta),'±',num2str(sdTheta),'º'],'FontSize',14);
h = gca; h.XTick = 0:45:360;
h.XTickLabel = {'0','45','90','135','180','225','270','315','360'};
% h = gca; h.XTick = 0:pi/4:2*pi;
% h.XTickLabel = {'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'};
end

function ang = atan360(X,Y)
if X < 0
    ang = pi - atan(Y/X);
elseif Y < 0
    ang = 2*pi + atan(Y/X);
elseif X < 0 && Y < 0
    ang = pi + atan(Y/X);
else
    ang = atan(Y/X);
end
end
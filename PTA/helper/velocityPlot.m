function velocityPlot(data,exptime)

r = zeros(size(data,1),1); ini = [r,r]; x = r; y = r; theta = r;
for i = 1:size(data,1)
    t = size(data,1)*exptime;
    ini(i,:) = data{i,1}(1,:);
    fin = data{i,1}(end,:);
    r(i) = pdist([ini;fin]);
    x(i) = (ini(2) - fin(2))/t;
    y(i) = (ini(1) - fin(1))/t;
    theta(i) = tan(y/x);
end

figure,quiver(ini(:,1),ini(:,2),x,y);
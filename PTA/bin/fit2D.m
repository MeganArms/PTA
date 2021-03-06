function [f,gof] = fit2D(obj,img)
% Fit PSF to ROI
Option = obj.Option;
R = Option.spotR;
ps = Option.pixelSize;
lambda = Option.wavelength;
NA = Option.na;
[y,x] = meshgrid(-R:R);
y = y(~isnan(img))*ps;
x = x(~isnan(img))*ps;
z = img(~isnan(img));

% PSF model
ft = fittype('A*exp((-(x-x0)^2-(y-y0)^2)/(2*sigma^2))+z0',...
    'independent',{'x','y'},'dependent','z');

% fit options
opts = fitoptions(ft);
opts.Display = 'off';
z_0 = Option.bg;
A_0 = img(R+1,R+1) - z_0;
x_0 = 0;
y_0 = 0;
sigma_0 = lambda/2/NA;
opts.StartPoint = [A_0 sigma_0 x_0 y_0 z_0];
opts.Lower = [0 0.1*R*ps -0.5*R*ps -0.5*R*ps z_0-200];
opts.Upper = [A_0+200 0.6*R*ps 0.5*R*ps 0.5*R*ps z_0+200];
for i=1:5
    if opts.Lower(i)>opts.Upper(i)
        f=1;
    end
end

[f,gof] = fit([x,y],z,ft,opts);

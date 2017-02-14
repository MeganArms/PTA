function circle_im = makecircle(R,r)

L = 2*R + 1; % length of a side
S = L^2; % Number of elements
im = zeros(L);
midp = [R+1,R+1];

lg = 1:L;
[XG,YG] = meshgrid(lg,lg);

subgrid = cat(3,XG,YG);
subs = reshape(subgrid,S,2);
distvect = cat(1,midp,subs);
a = squareform(pdist(distvect));
radii = a(1,:);
circlepoints = subs(radii > r - 3 & radii < r,:);
circlepoints(:,2) = circlepoints(:,2) - 1;
circleinds = sub2ind([L,L],circlepoints(:,1),circlepoints(:,2));
circle_im = zeros(size(im));
circle_im(circleinds) = 1;
end
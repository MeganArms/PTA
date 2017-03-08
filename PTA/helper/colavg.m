function bp = colavg(subimg)

% COLAVG (COLumn Average) is the function designed to work with COLFILT
% and a 'sliding' type 41x41 neighborhood. The purpose is to average the
% the kernel.

% Generalizable with padding of the lut's to the size of nhood and changing
% the "9"s to be m*n product of size of nhood.

% subimg(subimg > 0) = 1;
[~,n] = size(subimg);
% dispm = floor(m/3):floor(2*m/3);
% dispn = 1:floor(n/3);
% padm = heaviside(m - 9)*(m - 9);
% padn = heaviside(n - 1)*(n - 1);
nonzero = length(nonzeros(subimg));
if nonzero == 0
    nonzero = 1;
end
lut = 1/(nonzero^2)*ones(41);
bp = zeros(1,n);

lutcurrent = reshape(lut,[41^2,1]);
lutmatrix = repmat(lutcurrent,1,n);
outs = sum(subimg.*lutmatrix,1);
bp(1,:) = subimg(ceil(41^2/2),:).*outs; % Five is the middle pixel in column

end

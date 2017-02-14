function db = convert2double(img)

maxint = double(intmax(class(img)));
img = double(img);
db = img./maxint;

end
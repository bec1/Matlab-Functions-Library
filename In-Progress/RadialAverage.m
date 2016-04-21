function result=RadialAverage(img)
%% RADIALAVERAGE radially averages an input image

%% Apply a square centered crop
r0=round(centerfinder(img));
matrix = crop(img,r0,180);

%% Find distances to each pixel
[m,~] = size(matrix);
center = round(m/2); 
[i,j] = meshgrid(1:m,1:m);
i = i(:); j = j(:);
dist = round(sqrt((i-center).^2 + (j-center).^2));
[dist,y] = sort(dist);

%% Bin the distances
hh = hist(dist,max(dist)+1);

vec = matrix(:);
vec = vec(y); % sort the same way as dist

ini = 2;

result(1:max(dist))=0;

for k = 1:max(dist) %maybe this loop still can be optmized

   index = ini:ini+hh(k+1)-1; % this indexing is working perfectly
   result(k) = mean(vec(index));
   ini = max(index)+1;

end


end

function imgout=crop(img,r,d)
%% Crop a square subimage
w=floor(d/2);
x=floor(r(2));
y=floor(r(1));
imgout = img((x-w):(x+w),(y-w):(y+w));
end
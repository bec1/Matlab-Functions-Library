function avg = imgAvg(images)
for i=1:length(images)
    img(:,:,i) = loadrawimage(images{i});
end
avg = mean(img,3);
end

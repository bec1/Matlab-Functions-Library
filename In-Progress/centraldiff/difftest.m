function difftest
points = -10:.1:10;
sindata = sin(points);
cosdata = cos(points);
noises = 0:.1:2;

for i=1:length(noises)
    noiseamp = noises(i);
    sindatanoisy = sin(points) -noiseamp/2 + noiseamp*rand([1 length(points)]);
    residuals = cosdata - polydiff(points,sindatanoisy,10,2);
    

end

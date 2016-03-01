function difftest
points = -10:.1:10;
sindata = sin(points);
cosdata = cos(points);
noises = [linspace(0,.1, 50) logspace(0.1,1, 50) logspace(1,4,100)];


noiseamp = 1;
sindatanoisy = sin(points) -noiseamp/2 + noiseamp*rand([1 length(points)]);
cdiffdata = centraldiff(points,sindatanoisy,1,4);
pdiffdata = polydiff(points,sindatanoisy,2,2);

figure(1);
subplot(2,1,1)
plot(points,sindata)
hold all
plot(points,sindatanoisy)
subplot(2,1,2)
plot(points,cosdata)
hold all
plot(points,cdiffdata,'Color',[0.7 .7 .7])
plot(points,pdiffdata)


for i=1:length(noises)
    noiseamp = noises(i);
    sindatanoisy = sin(points) -noiseamp/2 + noiseamp*rand([1 length(points)]);
    
    residuals_p = cosdata - polydiff(points,sindatanoisy,10,2);
    msq_err_p(i) = sum(residuals_p.^2)/length(residuals_p);
    
    residuals_c = cosdata - centraldiff(points,sindatanoisy,1,4);
    msq_err_c(i) = sum(residuals_c.^2)/length(residuals_c);
    
end

SNR = 2./(noises);
figure(2); 
loglog(SNR,msq_err_p)
hold all
loglog(SNR,msq_err_c)

end

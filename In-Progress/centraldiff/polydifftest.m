function difftest
points = [linspace(-10,-10,200)];
sindata = sin(points);
cosdata = cos(points);
noises = [linspace(0,1, 10) logspace(0.1,1,10) logspace(1,4,20)];

ranges = 1:30;

for i=1:length(noises)
    noiseamp = noises(i);
    sindatanoisy = sin(points) -noiseamp/2 + noiseamp*rand([1 length(points)]);
    
    residuals_p = cosdata - polydiff(points,sindatanoisy,1,2);
    msq_err_p(i) = sum(residuals_p.^2)/length(residuals_p);
    
    
end

SNR = 2./(noises);
figure(2); 
loglog(SNR,msq_err_p)
hold all
loglog(SNR,msq_err_c)
grid on
xlabel('SNR')
ylabel('MSE')

end

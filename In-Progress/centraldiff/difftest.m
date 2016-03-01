function difftest
points = [linspace(-10,0,200) linspace(0, 10, 100)];
sindata = sin(points);
cosdata = cos(points);
noises = [linspace(0,1, 10) logspace(0.1,1,10) logspace(1,4,20)];


noiseamp = 0.5;
sindatanoisy = sin(points) -noiseamp/2 + noiseamp*rand([1 length(points)]);
cdiffdata = centraldiff_2(points,sindatanoisy);
pdiffdata = polydiff(points,sindatanoisy,10,2);

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


% for i=1:length(noises)
%     noiseamp = noises(i);
%     sindatanoisy = sin(points) -noiseamp/2 + noiseamp*rand([1 length(points)]);
%     
%     residuals_p = cosdata - polydiff(points,sindatanoisy,1,2);
%     msq_err_p(i) = sum(residuals_p.^2)/length(residuals_p);
%     
%     residuals_c = cosdata - centraldiff_2(points,sindatanoisy);
%     msq_err_c(i) = sum(residuals_c.^2)/length(residuals_c);
%     
% end
% 
% SNR = 2./(noises);
% figure(2); 
% loglog(SNR,msq_err_p)
% hold all
% loglog(SNR,msq_err_c)
% grid on
% xlabel('SNR')
% ylabel('MSE')

end

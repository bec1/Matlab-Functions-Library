function centers= sloshing(data)

for i=1:length(data)
    profiles(:,i) = sum(data(i).od,2);
    [xData, yData] = prepareCurveData( [], profiles(:,i) );

    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    opts.StartPoint = [228.524340783294 202 58.7049532423034];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    centers(i) = fitresult.b1;
end

end
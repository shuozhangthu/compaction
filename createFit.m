function [fitresult, gof] = createFit(depth, porosity)
%CREATEFIT(DEPTH,POROSITY)
%  Create a fit.
%
%  Data for 'exponential porosity' fit:
%      X Input : depth
%      Y Output: porosity
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 23-May-2018 09:45:20


%% Fit: 'exponential porosity'.
[xData, yData] = prepareCurveData( depth, porosity );

% Set up fittype and options.
ft = fittype( 'a+b*exp(-x/c)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.StartPoint = [0.525 0.125 175];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'exponential porosity' );
h = plot( fitresult, xData, yData );
legend( h, 'porosity vs. depth', 'exponential porosity', 'Location', 'NorthEast' );
% Label axes
xlabel depth
ylabel porosity
grid on


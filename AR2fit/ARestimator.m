function [parameters]=ARestimator(theta)
% ARestimator.m: function to estimate parameters for an
% autoregressive model with Gaussian disturbance, for an input angle series
% theta. Outputs parameters to a structure 'parameters' with fields for the
% fitted coefficients and p-values.
% Created April 10, 2020 by Ajay B. Limaye (ajay@virginia.edu).
% Last edited July 26, 2021 by Ajay B. Limaye (ajay@virginia.edu).

% For more see https://www.mathworks.com/help/econ/arma-models.html
% or execute graphically using: econometricModeler

% theta: coordinate series
mdl = arima('ARLags',[1 2],'constant',0); % AR-2 model with no constant
try
    mdlEstimate = estimate(mdl,theta,'display','off');
catch
    parameters.b1 = NaN; % if estimated variance model is invalid
    parameters.b2 = NaN;
    parameters.sigma = NaN;
    return
end

% save the output table
results = summarize(mdlEstimate);

% extract values and p-values for fitted AR coefficients
parameters.b1 = table2array(results.Table('AR{1}','Value'));
parameters.b1_pValue = table2array(results.Table('AR{1}','PValue'));

parameters.b2 = table2array(results.Table('AR{2}','Value'));
parameters.b2_pValue = table2array(results.Table('AR{2}','PValue'));

variance = table2array(results.Table('Variance','Value'));
parameters.sigma = sqrt(variance);
parameters.sigma_pValue = table2array(results.Table('Variance','PValue'));

end
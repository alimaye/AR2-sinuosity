function [parameters]=ARestimator(theta)
% ARestimator.m: function to estimate parameters for an
% autoregressive model with Gaussian disturbance, for an input angle series
% theta. Outputs parameters to a structure 'parameters' with fields for the
% fitted coefficients and p-values.
% Created April 10, 2020 by Ajay B. Limaye (ajay@virginia.edu).
% Last edited April 8, 2021 by Ajay B. Limaye (ajay@virginia.edu).

%https://www.mathworks.com/help/econ/arma-models.html
% or graphically using: econometricModeler

% theta: coordinate series
mdl = arima('ARLags',[1 2],'constant',0); % AR-2 model with no constant
try
    mdlEstimate = estimate(mdl,theta,'display','off');
catch
    parameters.AR1_coeff = NaN; % if estimated variance model is invalid
    parameters.AR2_coeff = NaN;
    parameters.sigma = NaN;
    return
end

% save the output table like this
results = summarize(mdlEstimate);

% extract values and p-values for AR coefficients
parameters.AR1_coeff = table2array(results.Table('AR{1}','Value'));
parameters.AR1_pValue = table2array(results.Table('AR{1}','PValue'));

parameters.AR2_coeff = table2array(results.Table('AR{2}','Value'));
parameters.AR2_pValue = table2array(results.Table('AR{2}','PValue'));

variance = table2array(results.Table('Variance','Value'));
parameters.sigma = sqrt(variance);
parameters.sigma_pValue = table2array(results.Table('Variance','PValue'));

end
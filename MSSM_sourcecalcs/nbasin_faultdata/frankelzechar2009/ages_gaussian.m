function [t, p_t] = ages_gaussian(ages, uncertainties, confidence, ...
    doPlot, minT, maxT)
% Syntax: 
% 
% age_gaussian(ages, uncertainties, confidence, doPlot, minT, maxT)
% 
% Given the sample ages and uncertainties and assuming a Gaussian 
% uncertainty distribution for each age, compute the summed, normalized 
% age probability density function, the corresponding c.d.f. and
% associated common statistics.  We return the pdf.  We may optionally plot
% the pdf and cdf, and we may show the plot over the (optionally) specified 
% range.
% 
% Input parameters:
%   ages                       	measured ages (in a)
%   uncertainty             	age uncertainties (standard  deviation)
%   confidence                  level of confidence in which you're
%                               interested (i.e., 0.95 will yield 95% 
%                               confidence)
%   doPlot                      Should we plot the resulting distributions?
%   minT                        What should be the minimum age be
%                               (for the plot)?
%   maxT                        What should be the maximum age be
%                               (for the plot)?
% 
% Output parameters:
%   t                           array of age values
%   p_t                         array representing p.d.f. of values in t
% 
% J. Douglas Zechar (zechar@ldeo.columbia.edu)

%%%%%%%%%%%%%%%%%%%%
% BEGIN AGE ANALYSIS
%%%%%%%%%%%%%%%%%%%%
assert(nargin == 6, 'Usage: age_gaussian(ages, uncertainties, confidence, doPlot, minT, maxT)');
assert(length(ages) == length(uncertainties), 'Age vector must be the same length as the uncertainties vector');
assert(min(ages) > 0, 'All sample ages must be greater than zero');
assert(min(uncertainties) > 0, 'All uncertainties must be greater than zero');
assert(confidence > 0 && confidence <= 1, 'Confidence must be between 0 and 1');

% We'll compute the age pdf from (the minimum age minus 5 of its 
% uncertainties) to (the maximum age plus 5 of its uncertainties)
t1 = min(ages); 
uncertaintyOfMinAge = uncertainties(find(ages == t1));
t1 = t1 - 5 * uncertaintyOfMinAge;
% if the minimum age would be less than zero (due to large uncertainties, or
% just b/c of young surfaces) we (yes, somewhat arbitrarily) set the min.
% age to be 1/10 the youngest age
if (t1 < 0)
    t1 = 0.1 * min(ages);
end
t2 = max(ages);
uncertaintyOfMaxAge = uncertainties(find(ages == t2));
t2 = t2 + 5 * uncertaintyOfMaxAge;

ageRange = t2 - t1;

numberOfTimeSteps = 1000;

% We'll finely discretize displacement
t = linspace(t1, t2, numberOfTimeSteps);
p_t = zeros(size(t)); % reserve storage for the summed normalized age pdf
D_t = zeros(size(t)); % reserve storage for the summed normalized age cdf

% Compute the composite pdf and cdf for all ages
numberOfAges = length(ages);
for i = 1:numberOfAges
    pdf = normpdf(t, ages(i), uncertainties(i)) / numberOfAges;
    cdf = normcdf(t, ages(i), uncertainties(i)) / numberOfAges;
    p_t = p_t + pdf;
    D_t = D_t + cdf;
%     plot(t, pdf, 'r'); % Plot each age's p.d.f.
end

h = max(p_t);

% Compute mean, median, mode, and bounds abt the median
t_mean = t * transpose(p_t ./ sum(p_t));
t_lowerBoundCDFValue = (1 - confidence) / 2;
t_upperBoundCDFValue = (1 + confidence) / 2;
t_lowerBound = t(find(D_t >= t_lowerBoundCDFValue, 1, 'first'));
t_upperBound = t(find(D_t >= t_upperBoundCDFValue, 1, 'first'));
t_median = t(find(D_t >= 0.5, 1, 'first'));
t_mode = t(find(p_t==max(p_t), 1, 'first'));

s1 = sprintf('==========================\nAGE STATISTICS (ka)\n==========================');
disp(s1);
s2 = sprintf('Mean: %.1f', t_mean / 1000);
disp(s2);
s6 = sprintf('Mode: %.2f', t_mode / 1000);
disp(s6);
s3 = sprintf('Median and %.2f%% bounds: %.1f +%.1f/-%.1f', (confidence * 100), t_median / 1000, (t_upperBound - t_median) / 1000, (t_median - t_lowerBound) / 1000);
disp(s3);

if (doPlot)
    % plot age cdf and pdf
    figure;
    hold on;
    title('Age density functions', 'FontSize', 16);
    [AX,H1,H2] = plotyy(t, p_t, t, D_t);
    set(get(AX(1),'Ylabel'),'String','Probability density');
    set(get(AX(1),'Ylabel'),'FontSize', 14);
    set(get(AX(2),'Ylabel'),'String','Cumulative density');
    set(get(AX(1),'Ylabel'),'FontSize', 14);
    set(get(AX(2),'Ylabel'),'FontSize', 14);
    set(H1, 'LineWidth', 3, 'Color', 'blue');
    set(H2, 'LineWidth', 3, 'LineStyle', '--', 'Color', [0 0.5 0]);
    axis(AX, 'square');
    ylim(AX(1), [0 1.5 * h]);
    ylim(AX(2), [0 1]);
    if (isempty(minT))
    else
        xlim(AX(1), [minT maxT]);
        xlim(AX(2), [minT maxT]);        
    end
    
    xlabel('Age (a)', 'FontSize', 14);
end
%%%%%%%%%%%%%%%%%%
% END AGE ANALYSIS
%%%%%%%%%%%%%%%%%%
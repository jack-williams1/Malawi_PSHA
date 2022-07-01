function [d, p_d] = displacements_gaussian(displacements, uncertainties, ...
    confidence, doPlot, minD, maxD)
% Syntax: 
% 
% displacements_gaussian(displacements, uncertainties, confidence, doPlot, minD, maxD)
% 
% Given the sample displacements and uncertainties and assuming a Gaussian 
% uncertainty distribution for each displacements, compute the summed, 
% normalized displacement probability density function, the corresponding c.d.f. and
% associated common statistics.  We return the pdf.  We may optionally plot 
% the pdf and cdf, and we may show the plot over the (optionally) specified range.
% 
% Input parameters:
%   displacements              	measured displacements, in m
%   uncertainty             	displacement uncertainties (standard  deviation, in m)
%   confidence                  level of confidence in which you're
%                               interested (i.e., 0.95 will yield 95% 
%                               confidence)
%   doPlot                      Should we plot the resulting distributions?
%   minD                        What should be the minimum offset be
%                               (for the plot)?
%   maxD                        What should be the maximum offset be
%                               (for the plot)?
% 
% Output parameters:
%   d                           array of displacement values
%   p_d                         array representing p.d.f. of values in d
% 
% J. Douglas Zechar (zechar@ldeo.columbia.edu)

%%%%%%%%%%%%%%%%%%%%
% BEGIN DISPLACEMENT ANALYSIS
%%%%%%%%%%%%%%%%%%%%
assert(nargin == 6, 'Usage: displacements_gaussian(displacements, uncertainties, confidence, doPlot, minD, maxD)');
assert(length(displacements) == length(uncertainties), 'Displacement vector must be the same length as the uncertainties vector');
assert(min(uncertainties) > 0, 'All uncertainties must be greater than zero');
assert(confidence > 0 && confidence <= 1, 'Confidence must be between 0 and 1');

% We'll compute the displacement pdf from (the minimum displacement minus 5 of its 
% uncertainties) to (the maximum displacement plus 5 of its uncertainties)
d1 = min(displacements - 5 * uncertainties);
%EDIT BY JW 28/04/21
% if the minimum displacement calculated above is <0, 
% set it to be 1/10 of min displacement (as is done in ages_gaussian)
% avoids -ve slip rates
if (d1 < 0)
    d1 = 0.1 * min(displacements);
end
d2 = max(displacements + 5 * uncertainties);

dispRange = d2 - d1;

numberOfDispSteps = 1000;

% We'll finely discretize displacement
d = linspace(d1, d2, numberOfDispSteps);
p_d = zeros(size(d)); % reserve storage for the summed normalized displacement pdf
D_d = zeros(size(d)); % reserve storage for the summed normalized displacement cdf

% Compute the composite pdf and cdf for all displacements
numberOfDisplacements = length(displacements);
for i = 1:numberOfDisplacements
    pdf = normpdf(d, displacements(i), uncertainties(i)) / numberOfDisplacements;
    cdf = normcdf(d, displacements(i), uncertainties(i)) / numberOfDisplacements;
    p_d = p_d + pdf;
    D_d = D_d + cdf;
%     plot(d, pdf, 'r'); % Plot each displacement's p.d.f.
end

h = max(p_d);

% Compute mean, median, mode, and bounds abt the median
d_mean = d * transpose(p_d ./ sum(p_d));
d_lowerBoundCDFValue = (1 - confidence) / 2;
d_upperBoundCDFValue = (1 + confidence) / 2;
d_lowerBound = d(find(D_d >= d_lowerBoundCDFValue, 1, 'first'));
d_upperBound = d(find(D_d >= d_upperBoundCDFValue, 1, 'first'));
d_median = d(find(D_d >= 0.5, 1, 'first'));
d_mode = d(find(p_d==max(p_d), 1, 'first'));

s1 = sprintf('==========================\nDISPLACEMENT STATISTICS (m)\n==========================');
disp(s1);
s2 = sprintf('Mean: %.1f', d_mean);
disp(s2);
s6 = sprintf('Mode: %.1f', d_mode);
disp(s6);
s3 = sprintf('Median and %.2f%% bounds: %.1f +%.1f/-%.1f', (confidence * 100), d_median, (d_upperBound - d_median), (d_median - d_lowerBound));
disp(s3);

if (doPlot)
    % plot displacement cdf and pdf
    figure;
    hold on;
    title('Displacement density functions', 'FontSize', 16);
    [AX,H1,H2] = plotyy(d, p_d, d, D_d);
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
    if (isempty(minD))
    else
        xlim(AX(1), [minD maxD]);
        xlim(AX(2), [minD maxD]);        
    end
    xlabel('Displacement (m)', 'FontSize', 14);
end
%%%%%%%%%%%%%%%%%%
% END DISPLACEMENT ANALYSIS
%%%%%%%%%%%%%%%%%%
function [d, p_d] = displacement_arbitrary(displacement, density, ...
    confidence, doPlot, minD, maxD)
% Syntax: 
% 
% displacement_abitrary(displacement, density, confidence, doPlot, minD, maxD)
% 
% For the specified displacment probability density, we compute the most
% the common statistics for this distribution and return the pdf.
% We may optionally plot the pdf and cdf, and we may show the plot 
% over the (optionally) specified range.
% 
% Input parameters:
%   displacement                array of displacement values of interest, in m
%   density                 	probability density values corresponding to
%                               the specified displacements
%   confidence                  level of confidence in which you're
%                               interested (i.e., 0.95 will yield 95% 
%                               confidence)
%   doPlot                      Should we plot the resulting distributions?
%   minT                        What should be the minimum offset be
%                               (for the plot)?
%   maxT                        What should be the maximum offset be
%                               (for the plot)?
%
% Output parameters:
%   d                           array of displacement values
%   p_d                         array representing p.d.f. of values in d
% 
% J. Douglas Zechar (zechar@ldeo.columbia.edu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN DISPLACEMENT ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(nargin == 6, 'Usage: displacement_abitrary(displacement, density, confidence, doPlot, minD, maxD)');
assert(length(displacement) == length(density), 'Displacement vector must be the same length as the density vector');
assert(confidence > 0 && confidence <= 1, 'Confidence must be between 0 and 1');

d = displacement;
p_d = density;

d1 = min(displacement);
d2 = max(displacement);
displacementRange = d2 - d1;

% Compute displacment cdf
D_d = zeros(length(p_d), 1);
D_d(1) = p_d(1);
for i=2:length(d)
    D_d(i) = D_d(i - 1) + p_d(i);
end
D_d = D_d / (sum(p_d));

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
    xlabel('Displacement (m)', 'FontSize', 14);
    if (isempty(minD))
    else
        xlim(AX(1), [minD maxD]);
        xlim(AX(2), [minD maxD]);        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END DISPLACEMENT ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
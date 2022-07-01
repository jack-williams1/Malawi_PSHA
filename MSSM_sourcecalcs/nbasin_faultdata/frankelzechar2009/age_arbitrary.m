function [t, p_t] = age_arbitrary(age, density, confidence, ...
    doPlot, minT, maxT)
% Syntax: 
% 
% age_abitrary(age, density, confidence, doPlot, minT, maxT)
% 
% For the specified age probability density, we compute the most
% the common statistics for this distribution and return the pdf.
% We may optionally plot the pdf and cdf, and we may show the plot 
% over the (optionally) specified range.
% 
% Input parameters:
%   age                         array of age values of interest, in a
%   density                 	probability density values corresponding to
%                               the specified ages
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
assert(nargin == 6, 'Usage: age_arbitrary(age, density, confidence, doPlot, minT, maxT)');
assert(length(age) == length(density), 'Age vector must be the same length as the density vector');
assert(confidence > 0 && confidence <= 1, 'Confidence must be between 0 and 1');
assert(min(age) > 0, 'All ages must be greater than zero');

t = age;
p_t = density;

t1 = min(age);
t2 = max(age);
ageRange = t2 - t1;

% Compute age cdf
D_t = zeros(length(p_t), 1);
D_t(1) = p_t(1);
for i=2:length(t)
    D_t(i) = D_t(i - 1) + p_t(i);
end
D_t = D_t / (sum(p_t));

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
    plot([t1 t1 t2 t2], [0 h h 0], 'Color', 'blue', 'LineWidth', 3);
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
    xlabel('Age (a)', 'FontSize', 14);
    if (isempty(minT))
    else
        xlim(AX(1), [minT maxT]);
        xlim(AX(2), [minT maxT]);        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END AGE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
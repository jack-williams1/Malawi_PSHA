function [v, p_v, v_mean,v_lowerBound ,v_upperBound] = slip_rate(ages, ageDensity, displacements, ...
    displacementDensity, confidence, doPlot, minV, maxV)

%Noted amended by jw 25/04/21 so that function also returns v_mean,
%v_lowerbound and v_upperbound

% Syntax: 
% 
% slip_rate(ages, ageDensity, displacements, displacementDensity, ...
% confidence, doPlot, minV, maxV)
% 
% Using the specified displacement distribution and the age distribution, 
% we compute the corresponding slip rate pdf and cdf and the associated 
% common statistics.  We may optionally plot the pdf and cdf, and we may
% show the plot over the (optionally) specified slip rate range.
% 
% Note: Slip rate pdf computation is based on Eq'n 29 from Peter Bird's 
% 2007 Geosphere paper.
% 
% Input parameters:
%   ages                       	ages of interest (in a)
%   ageDensity                  probability density values corresponding to
%                               the specified ages
%   displacements             	displacements of interest (in m)
%   displacementDensity         probability density values corresponding to
%                               the specified ages
%   confidence                  confidence region in which you're
%                               interested (i.e., 0.95 will yield 95% 
%                               confidence region about the median)
%   doPlot                      Should we plot the resulting distributions?
%   minV                        What should be the minimum slip rate be
%                               (for the plot)?
%   maxV                        What should be the maximum slip rate be
%                               (for the plot)?
% 
% Output parameters:
%   v                           array of slip rate values
%   p_v                         array representing p.d.f. of values in v
% 
% J. Douglas Zechar (zechar@ldeo.columbia.edu)

%%%%%%%%%%%%%%%%%%%%
% BEGIN RATE ANALYSIS
%%%%%%%%%%%%%%%%%%%%
assert(nargin == 8, 'Usage: slip_rate(ages, ageDensity, displacements, displacementDensity, confidence, doPlot, minV, maxV)');
assert(length(ages) == length(ageDensity), 'Age vector must be the same length as the age pdf vector');
assert(length(displacements) == length(displacementDensity), 'Displacement vector must be the same length as the displacement pdf vector');
assert(min(ages) > 0, 'All ages must be greater than zero');
assert(confidence > 0 && confidence <= 1, 'Confidence must be between 0 and 1');

% Minimum possible slip rate is given by the minimum displacement divided 
% by the maximum age
v1 = min(displacements) / max(ages) * 1000; % multiply by 1000 to get to mm/a

% Maximum possible slip rate is given by the maximum displacement divided 
% by the minimum age.  Here, we limit the maximum slip rate to be 100
% mm/a; otherwise the slip rate range can be unrealistically large and we
% may undersample the distribution
v2 = min(max(displacements) / min(ages) * 1000, 100); % multiply by 1000 to get to mm/a
slipRateRange = v2 - v1;
numberOfSlipRateSteps = 1000;

assert(v2 > v1, 'Maximum slip rate in range must be greater than minimum slip rate.');
t = ages;
p_t = ageDensity;
d = displacements;
p_d = displacementDensity;

% We'll finely discretize slip rate space and make room to store the p.d.f.
% and cumulative probability
v = linspace(v1, v2, numberOfSlipRateSteps);
p_v = zeros(size(v));
D_v = zeros(size(v));

% Compute the slip rate pdf and cdf for all values of slip rate
for i=1:length(v)
    rate = v(i);
    for j=1:length(t)
        age = t(j);
        displacement = rate * age / 1000; % divide by 1000 to get to m
        displacementIndex = find(d >= displacement, 1);
        if isempty(displacementIndex) == 0
                p_v(i) = p_v(i) + p_t(j) * age * p_d(displacementIndex);
        end
    end
    if (i > 1)
        D_v(i) = D_v(i - 1) + p_v(i);
    else
        D_v(i) = p_v(i);
    end
end
D_v = D_v / sum(p_v);
p_v = p_v / sum(p_v);

% Find the maximum probability density (used for plotting)
h = max(p_v);

% Compute mean, median, mode, and bounds abt the median
v_mean = v * transpose(p_v ./ sum(p_v));
v_lowerBoundCDFValue = (1 - confidence) / 2;
v_upperBoundCDFValue = (1 + confidence) / 2;
v_lowerBound = v(find(D_v >= v_lowerBoundCDFValue, 1, 'first'));
v_upperBound = v(find(D_v >= v_upperBoundCDFValue, 1, 'first'));
v_median = v(find(D_v >= 0.5, 1, 'first'));
v_mode = v(find(p_v == max(p_v), 1, 'first'));

s1 = sprintf('========================\nSLIP RATE STATISTICS (mm/a)\n========================');
disp(s1);
s2 = sprintf('Mean: %.1f', v_mean);
disp(s2);
s6 = sprintf('Mode: %.1f', v_mode);
disp(s6);
s3 = sprintf('Median and %.2f%% bounds: %.1f +%.1f/-%.1f', (confidence * 100), v_median, (v_upperBound - v_median), (v_median - v_lowerBound));
disp(s3);

if (doPlot)
    % plot slip rate cdf and pdf
    figure;
    hold on;
    title('Slip rate density functions', 'FontSize', 16);
    [AX,H1,H2] = plotyy(v, p_v, v, D_v);
    set(get(AX(1),'Ylabel'),'String','Probability density');
    set(get(AX(1),'Ylabel'),'FontSize', 14);
    set(get(AX(2),'Ylabel'),'String','Cumulative probability');
    set(get(AX(1),'Ylabel'),'FontSize', 14);
    set(get(AX(2),'Ylabel'),'FontSize', 14);
    set(H1, 'LineWidth', 3, 'Color', 'blue');
    set(H2, 'LineWidth', 3, 'LineStyle', '--', 'Color', [0 0.5 0]);
    axis(AX, 'square');
    ylim(AX(1), [0 1.5 * h]);
    ylim(AX(2), [0 1]);

    % If the user specified some limits for the horizontal axis, use them
    if (isempty(minV))
    else
        xlim(AX(1), [minV maxV]);
        xlim(AX(2), [minV maxV]);        
    end
    xlabel('Slip rate (mm/a)', 'FontSize', 14);
end

% Check to see if the slip rate is well-constrained according to the
% criterion of Bird (2007): compute the width of the 2-sigma confidence
% interval and see if it is less than the median
confidence = 0.9545;
v_lowerBoundCDFValue = (1 - confidence) / 2;
v_upperBoundCDFValue = (1 + confidence) / 2;
v_lowerBound_1 = v(find(D_v >= v_lowerBoundCDFValue, 1, 'first'));
v_upperBound_1 = v(find(D_v >= v_upperBoundCDFValue, 1, 'first'));
v_median = v(find(D_v >= 0.5, 1, 'first'));
v_widthOfConfidenceInterval = v_upperBound_1 - v_lowerBound_1;
if (v_widthOfConfidenceInterval > v_median)
    s4 = sprintf('NOTE: Slip rate is not well-constrained (Median = %.1f, width of 95.45 percent confidence interval is %.1f)', v_median, v_widthOfConfidenceInterval);
    disp(s4);
end

%%%%%%%%%%%%%%%%%%
% END AGE ANALYSIS
%%%%%%%%%%%%%%%%%%
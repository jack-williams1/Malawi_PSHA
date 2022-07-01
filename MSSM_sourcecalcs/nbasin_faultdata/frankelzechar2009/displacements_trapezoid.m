function [d, p_d] = displacements_trapezoid(minDisplacements, ...
    minPreferredDisplacements, maxPreferredDisplacements, maxDisplacements, ...
    confidence, doPlot, minD, maxD)
% Syntax: 
% 
% displacements_trapezoid(minDisplacements, minPreferredDisplacements, ...
%   maxPreferredDisplacements, maxDisplacements, confidence, doPlot, ...
%   minD, maxD)
% 
% Given multiple trapezoid displacement pdfs, compute the summed, 
% normalized displacement probability density function, the corresponding c.d.f. and
% associated common statistics.  We return the pdf.  We may optionally plot 
% the pdf and cdf, and we may show the plot over the (optionally) specified range.
% 
% Input parameters:
%   minDisplacements            array of absolute min displacements you
%                               feel are reasonable, in m
%   minPreferredDisplacements	preferred minimum displacements, in m
%   maxPreferredDisplacements	preferred maximum displacements, in m
%   maxDisplacements            absolute max displacements you feel are 
%                               reasonable
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
%   p_d                         array representing pdf of values in d
% 
% J. Douglas Zechar (zechar@ldeo.columbia.edu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN DISPLACEMENT ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(nargin == 8, 'Usage: displacements_trapezoid(minDisplacements, minPreferredDisplacements, maxPreferredDisplacements, maxDisplacements, confidence, doPlot, minD, maxD)');
assert(length(minDisplacements) == length(minPreferredDisplacements), ...
    'All parameter vectors must be of the same length');
assert(length(minPreferredDisplacements) == length(maxPreferredDisplacements), ...
    'All parameter vectors must be of the same length');
assert(length(maxPreferredDisplacements) == length(maxDisplacements), ...
    'All parameter vectors must be of the same length');
assert(sum(minDisplacements >= minPreferredDisplacements) == 0, ...
       'minDisplacements must be less than minPreferredDisplacements');
assert(sum(minPreferredDisplacements >= maxPreferredDisplacements) == 0, ...
    'minPreferredDisplacements must be less than maxPreferredDisplacements');
assert(sum(maxPreferredDisplacements >= maxDisplacements) == 0, ...
    'maxPreferredDisplacements must be less than maxDisplacements');
assert(confidence > 0 && confidence <= 1, 'Confidence must be between 0 and 1');

minDisp = min(minDisplacements);
maxDisp = max(maxDisplacements);

displacementRange = maxDisp - minDisp;
numberOfDisplacementSteps = 1000;

% We'll finely discretize displacement
d = linspace(minDisp, maxDisp, numberOfDisplacementSteps);
p_d = zeros(size(d));
D_d = zeros(size(d));
maxH = 0;

for j = 1:length(minDisplacements)
    % To compute each probability density, we use some simple geometry and
    % the constraint that the area under each pdf must sum to unity.
    d1 = minDisplacements(j);
    d2 = minPreferredDisplacements(j);
    d3 = maxPreferredDisplacements(j);
    d4 = maxDisplacements(j);

    l1 = d2 - d1;
    l2 = d3 - d2;
    l3 = d4 - d3;
    h = 1 / (0.5 * l1 + l2 + 0.5 * l3); % maximum displacement pdf value

    % Compute the displacement pdf and cdf for all values of displacement
    for i=1:length(d)
        dValue = minDisp + i * (displacementRange / numberOfDisplacementSteps);
        if (dValue > d1 && dValue <= d2)
            p_d(i) = p_d(i) + (h / l1) * (dValue - d1);
        end
        if (dValue > d2 && dValue < d3)
            p_d(i) = p_d(i) + h;
        end
        if (dValue >= d3 && dValue < d4)
            p_d(i) = p_d(i) + (h / l3) * (d4 - dValue);
        end
    end
    if (h > maxH)
        maxH = h;
    end
end

% Normalize the pdf
p_d = p_d / length(minDisplacements);

for i = 1:length(d)
    if (i > 1)
        D_d(i) = D_d(i - 1) + p_d(i);
    else
        D_d(i) = p_d(i);
    end
end
D_d = D_d / sum(p_d);

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
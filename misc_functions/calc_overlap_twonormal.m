% numerical integral of the overlapping area of two normal distributions:
% s1,s2...sigma of the normal distributions 1 and 2
% mu1,mu2...center of the normal distributions 1 and 2
% xstart,xend,xinterval...defines start, end and interval width
% example: [overlap] = calc_overlap_twonormal(2,2,0,1,-10,10,0.01)

function [overlap2, h1, h2, h3] = calc_overlap_twonormal(s1,s2,mu1,mu2,xstart,xend,xinterval,pval,fault_name,std_sr)
%clf
x_range=xstart:xinterval:xend;
h1=plot(x_range,normpdf(x_range,mu1,s1)','k-','linewidth',1.2);hold on

hold on

h3=area(x_range,min([normpdf(x_range,mu1,s1)' normpdf(x_range,mu2,s2)']'),'FaceColor',[0.75 0.75 0.75]);
h2=plot(x_range,normpdf(x_range,mu2,s2)','r-','linewidth',1.2);
overlap=cumtrapz(x_range,min([normpdf(x_range,mu1,s1)' normpdf(x_range,mu2,s2)']'));
overlap2 = overlap(end);
h4=plot([mu2-2*std_sr ,mu2-2*std_sr],[0 4],'k--');
h5=plot([mu2+2*std_sr ,mu2+2*std_sr],[0 4],'k--');
xlabel('Slip Rate mm/yr'); ylabel('f(x)'); xlim([0 1.2]); ylim([0 4]);
set(gca,'FontSize',13); axis square


if pval>0.05
    h_test=('Accepted');
    else
    h_test=('Rejected');
end

if pval>0.05
    pval_str=(num2str(pval,2));
    else
    pval_str=('<0.05');
end

title(fault_name,'fontweight','normal');
subtitle(['Overlap: ',num2str(overlap2,2), ', p = ', pval_str]);
end
%legend([h1 h2 h3],{'simulation slip rate','seismic reflector slip rate',['overlap = ',num2str(overlap2,2)]});

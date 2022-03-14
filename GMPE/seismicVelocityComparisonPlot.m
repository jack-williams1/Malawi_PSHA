% PLOT VARIOUS VELOCITY PROFILES (VS) AGAINST EACH OTHER FOR GMPE SELECTION

%read and sort data

ebinger2019=readtable('MalawiSeismicVelocityComparison.xlsx','Sheet','Ebinger2019_clean');
maguire1994=readtable('MalawiSeismicVelocityComparison.xlsx','Sheet','Maguire1994');
boore2016=readtable('MalawiSeismicVelocityComparison.xlsx','Sheet','Boore2016');
stevens2021=readtable('MalawiSeismicVelocityComparison.xlsx','Sheet','Stevens2021');
%% PLOT DATA
figure (1)

%subplot(1,2,1)
plot(ebinger2019.vs,ebinger2019.depth,'r-','LineWidth',1.5); hold on
plot(stevens2021.vs,stevens2021.depth,'b-','LineWidth',1.5); hold on
%plot(maguire1994.vs,maguire1994.depth,'m-','LineWidth',1.5); hold on
plot(boore2016.vs,boore2016.depth,'k-','LineWidth',1.5); hold on
set(gca, 'YDir','reverse'); pbaspect([1 2 1]); set(gca,'fontsize',13); hold on;
ylabel('Depth (km)'); xlabel('vs (km/s)'); xlim([0 5]); ylim([0 40]); legend({'Ebinger2019','Stevens2021','Boore2016'},'Location','southwest');


%Option to plot seismic velocities in 0-10 km depth interval
%{

t1 = title('(a)','fontweight','normal');
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]);

subplot(1,2,2)
plot(ebinger2019.vs,ebinger2019.depth,'r-','LineWidth',1.5); hold on
plot(stevens2021.vs,stevens2021.depth,'b-','LineWidth',1.5); hold on
%plot(maguire1994.vs,maguire1994.depth,'m-','LineWidth',1.5); hold on
plot(boore2016.vs,boore2016.depth,'k-','LineWidth',1.5); hold on
 set(gca, 'YDir','reverse'); pbaspect([1 2 1]); set(gca,'fontsize',13); hold on;
ylabel('Depth (km)'); xlabel('vs (km/s)'); xlim([0 5]); ylim([0 10]); legend({'Ebinger2019','Stevens2021','Boore2016'},'Location','southwest')

t2 = title('(b)','fontweight','normal');
set(t2, 'horizontalAlignment', 'left','units','normalized');
h2 = get(t1, 'position');
set(t2, 'position', [-0.2 h1(2) h1(3)]);

%}
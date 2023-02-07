%% Plot Length- and Layer-Limited source aspect ratios

close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/MSSM_sourcecalcs']);

[num_raw1,txt_raw1,raw_raw1]=xlsread('MSSM.xlsx','MSSM_AdaptedSources');
[num_raw2,txt_raw2,raw_raw2]=xlsread('MSSM.xlsx','FaultGeometry');
[num_raw3,txt_raw3,raw_raw3]=xlsread('MSSM.xlsx','MultiFaultGeometry');

source_dim=zeros(length(num_raw1(:,1)),3);
count=0;
for ii=1:length(num_raw1(:,1))
    
    if num_raw1(ii,1)<600
       tmp_indx=find(num_raw1(ii,1)==num_raw2(:,1));
       source_dim(ii,1)=num_raw2(tmp_indx,6); 
    else
       tmp_indx=find(num_raw1(ii,1)==num_raw3(:,1));
       source_dim(ii,1)=num_raw3(tmp_indx,4); 
       
    end
       source_dim(ii,2)=num_raw1(ii,9)/source_dim(ii,1);
       source_dim(ii,3)=num_raw1(ii,11)/source_dim(ii,1);
       
       if (source_dim(ii,2))>source_dim(ii,3)
           count=count+1;
       end
end

%Length-width scaling from eq1 in PSHA manuscript
len_tmp=[5:1:200]'; width_tmp=zeros(length(len_tmp),1);
tmp2=35/sin(deg2rad(53));

for jj=1:length(len_tmp)

    tmp1=(17.5.*(len_tmp(jj)*1000)^0.667)/1000;
    
    if tmp1<tmp2
        width_tmp(jj)=tmp1;
    else
        width_tmp(jj)=tmp2;
    end
end

%% Plot length vs width for layer and length limited sources

figure(101);

plot(len_tmp,width_tmp,'k--');hold on
plot(source_dim(:,1),source_dim(:,2),'bv'); hold on
plot(source_dim(:,1),source_dim(:,3),'r^');

xlim([0,160]);ylim([0,50]);
xlabel('Length (km)'),ylabel('Width (km)'); set(gca,'fontsize',11);
legend('Scaling from Eq. 1','Length-limited sources','Layer-limited sources','Location','southeast');
axis square
%% Code for indexing Malawi-Mozambique border in Lake Malawi
%Ensures border is visible if lake background is different 
%from rest of map (e.g. PSHA maps)

load map_data_EastAfrica

LakeMalawi = shaperead('malawi_lake.shp');
LakeMalawiCoord = [LakeMalawi.Y(1,1:end-1)' LakeMalawi.X(1,1:end-1)'];

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
load('syncat_PSHA_MSSD_input','Region');
%index border coordindates in Lake

indx=1;

for ii=1:length(MapData2)

    if inpolygon(MapData2(ii,2),MapData2(ii,1),LakeMalawiCoord(:,2),LakeMalawiCoord(:,1))~=0
        LakeBorder_indx(indx)=ii; indx=indx+1;
    end
end

%Group lake border into the individual segments that are indexed
kk=1; ll=1;

for jj=1:length(LakeBorder_indx)-1
    
   if LakeBorder_indx(jj+1)-LakeBorder_indx(jj)==1 
      LakeBorder{kk}(ll,:)=MapData2(LakeBorder_indx(jj),:); ll=ll+1;
      
   else    
      kk=kk+1; ll=1; 
      LakeBorder{kk}(ll,:)=MapData2(LakeBorder_indx(jj),:);
      
   end
    
end

%% Index only parts of Lake Border that are of interest

LakeMalawiBorder=vertcat(LakeBorder{10},LakeBorder{2},LakeBorder{1},[-13.4788 34.8611]);

%Plot figure to check
figure(14);

plot(MapData2(:,2),MapData2(:,1),'k-','LineWidth',1); hold on
for kk=1:length(LakeBorder)    
plot(LakeBorder{kk}(:,2),LakeBorder{kk}(:,1),'b-');hold on
end

plot(LakeMalawiBorder(:,2),LakeMalawiBorder(:,1),'r-','LineWidth',3)
axis equal; axis(Region);

%%
save('LakeMalawiBorder','LakeMalawiBorder');
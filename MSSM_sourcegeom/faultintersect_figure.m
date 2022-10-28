%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   %% FAULT INTERSECTION FIGURE   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Same as fault geom_1 to create grid points
% Then checks fault intersection and creates figures

close all

%Load functions needed for fault gemoetry calculations
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/misc_functions']);
addpath([mydir(1:idcs(end)-1) '/gis_files']);

load MSSM_sources

ellipsoid = almanac('earth','wgs84','meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);

%% Group sec info by fault and create polygons

fault_end=zeros(length(num_fault),11);

for i=1:length(num_fault)
   
    if num_fault(i,13)>0 %index fault by sections 
        
    sec_index = find(num_fault(i,1)==num_sec(:,15));
    
    %index sec 1) lat 2)long 3)length 4) sec strike 5) dip 6) h_level
    %7) h_top 8)% fault_id  %8) sec_id 9) 10) fault strike 11) width
    
        if num_fault(i,23)>0
            %select multi_fault width
            fault_geom_{i} =  vertcat(num_sec(sec_index,[6:10,12:13,15,1,16,19]));
        else
            %select fault width
            fault_geom_{i} =  vertcat(num_sec(sec_index,[6:10,12:13,15,1,16,17]));
        end
        
    clear sec_index
    
    else %single section fault geom info taken from num_fault
    
         if num_fault(i,23)>0
            %select multi_fault width
            fault_geom_{i} =  num_fault(i,[4:8,10:11,1,1,7,24]);
         else
            fault_geom_{i} =  num_fault(i,[4:8,10:11,1,1,7,12]); 
         end
    end
 
    %Sort fault sections to find end point
   
    if fault_geom_{i}(1,10) >90 && fault_geom_{i}(1,10) <270
        
    % For right hand rule, faults with strike 090-270 will have end
    % point at southern end 
         fault_geom_{i} = sortrows(fault_geom_{i},1);
    else    
        
    % For right hand rule, other strikes will have end
    % point at northern end     
        fault_geom_{i} = sortrows(fault_geom_{i},1,'descend');
    end   
    
    %Manual correction for Mwanza as it is not unidirectional N to S
    if i == 2
        fault_geom_{i}=fault_geom_{i}([1 2 4 3 5 6 7],:);
    end
    
    %Find end point of fault (info not included in num_sec)
    fault_end(i,1:2) = reckon(fault_geom_{i}(1,1),fault_geom_{i}(1,2),km2deg(fault_geom_{i}(1,3)),fault_geom_{i}(1,4));
    %Assign fault end properties of closest section end point
    fault_end(i,:)=horzcat(fault_end(i,1:2),fault_geom_{i}(1,3:11));
    %Combine fault end and rest of geometry into one cell array
    geometry_tmp{i} = vertcat(fault_end(i,:),fault_geom_{i});   
    
    %Project fault through fault or multi_fault width and strike to get geometry bottom
    for jj=1:height(geometry_tmp{i})
        
        if geometry_tmp{i}(jj,7)==0
        %if fault propagates to surface
        geometry_top{i}(jj,:)=geometry_tmp{i}(jj,[1,2,7]);
        geometry_bottom{i}(jj,1:2) = reckon(geometry_top{i}(jj,1),geometry_top{i}(jj,2),km2deg(geometry_tmp{i}(jj,11)*cos(geometry_tmp{i}(jj,5)*pi/180)),geometry_tmp{i}(jj,10)+90);
        geometry_bottom{i}(jj,3) = geometry_tmp{i}(jj,11)*sin(geometry_tmp{i}(jj,5)*pi/180);
      
        else
        %if fault mapped below surface, project back to surface up-dip
        geometry_top{i}(jj,1:2) = reckon(geometry_tmp{i}(jj,1),geometry_tmp{i}(jj,2),km2deg((geometry_tmp{i}(jj,6)-geometry_tmp{i}(jj,7))/tan(geometry_tmp{i}(jj,5)*pi/180)),geometry_tmp{i}(jj,10)-90);
        geometry_top{i}(jj,3) = geometry_tmp{i}(jj,7);
        %project base, where lateral extent is difference between top of
        %fault and base of crust
        geometry_bottom{i}(jj,1:2) = reckon(geometry_top{i}(jj,1),geometry_top{i}(jj,2),km2deg((geometry_tmp{i}(jj,11)*cos(geometry_tmp{i}(jj,5)*pi/180))-(geometry_top{i}(jj,3)/tan(geometry_tmp{i}(jj,5)*pi/180))),geometry_tmp{i}(jj,10)+90);
        geometry_bottom{i}(jj,3) = (geometry_tmp{i}(jj,11)- geometry_top{i}(jj,3))*sin(geometry_tmp{i}(jj,5)*pi/180);    

        end
    end
    
    fault_geometry_points{i}=[geometry_top{i}; flipud(geometry_bottom{i}); geometry_top{i}(1,:)];
    fault_geometry_points{i}=horzcat(fault_geometry_points{i},ones(height(fault_geometry_points{i}),1)*num_fault(i,1));
   clear  geometry_tmp{i}
end

%% Create fault grids

dist={};%Min distance between grid point and surface fault trace
idx_c={};%Fault segment point is associated with
x_d_min={};
y_d_min={};
is_vertex={};

simu_count = 1;
clock_old = cputime;

for i=1:length(num_fault)
   
    if i== 20*simu_count
        clock_new = cputime;
        disp(['Simulation cycle is: ',num2str(i),' out of ',num2str(length(num_fault)),' & Required time is: ',num2str(clock_new-clock_old)]);
        clock_old = clock_new;
        simu_count = simu_count + 1;
    end
    
 dip=num_fault(i,8);  
  
%Fault verticies in utm (so in same metre units like depth)
 [faultx_{i}, faulty_{i}] = deg2utm(fault_geometry_points{i}(:,1),fault_geometry_points{i}(:,2));
%Find verticies in utm of fault surface trace 
 [faultxs_{i}, faultys_{i}] = deg2utm(geometry_top{i}(:,1),geometry_top{i}(:,2));
 
%First grid fault in 2D, but with number of points needed in 3D 
% in a 1x1 km grid (where the scale is in metres)
 fault_grid{i} = polygrid(faultx_{i}, faulty_{i},1e-6);
 
    for jj=1:length(fault_grid{i})
    %Find horizontal distance between each grid point and fault surface
    %trace
    [dist{i}(jj), x_d_min{i}(jj), y_d_min{i}(jj), is_vertex{i}(jj),idx_c{i}(jj)]=p_poly_dist1(fault_grid{i}(jj,1),fault_grid{i}(jj,2),faultxs_{i}, faultys_{i},false);
    dist{i}(jj)=dist{i}(jj)/1000;
    %Convert grid points back to lat and long
    [faultgrid_wgs84{i}(jj,1), faultgrid_wgs84{i}(jj,2)]=utm2deg(fault_grid{i}(jj,1),fault_grid{i}(jj,2),'36 L');
    %Find depth of point on fault surface
    faultgrid_wgs84{i}(jj,3)=num_fault(i,11)+(dist{i}(jj)*tan(dip*pi/180));
    
    %Assign fault grid coordinates section ids
    %Only works if order of fault sections in geom cell array 
    % same as searched through in p_poly_dist1

    %Note some short (<1 km) fault sections don't have points 
    %assigned to them
    if fault_geom_{i}(1,10) >90 && fault_geom_{i}(1,10) <270
    
        for kk=1:height(fault_geom_{i})
            if idx_c{i}(jj)==kk
               faultgrid_wgs84{i}(jj,4) = fault_geom_{i}(kk,8);%index fault_id
               faultgrid_wgs84{i}(jj,5) = fault_geom_{i}(kk,9);%index sec_id
               faultgrid_wgs84{i}(jj,6) = idx_c{i}(jj);%index number of sections
            end
            
            if idx_c{i}(jj)==0 %bug where some sections near fault ends not assigned to sections
               faultgrid_wgs84{i}(jj,4) = fault_geom_{i}(1,8);
               faultgrid_wgs84{i}(jj,5) = fault_geom_{i}(1,9);
               faultgrid_wgs84{i}(jj,6) = 1;
            end
            
        end
    
   else
        
        for kk=1:height(fault_geom_{i})
            if idx_c{i}(jj)==(height(fault_geom_{i})+1-kk)
                faultgrid_wgs84{i}(jj,4) = fault_geom_{i}((height(fault_geom_{i})+1-kk),8);
                faultgrid_wgs84{i}(jj,5) = fault_geom_{i}((height(fault_geom_{i})+1-kk),9);
                faultgrid_wgs84{i}(jj,6) = idx_c{i}(jj);
            end
            
            if idx_c{i}(jj)==0
                 faultgrid_wgs84{i}(jj,4) = fault_geom_{i}(1,8);
                 faultgrid_wgs84{i}(jj,5) = fault_geom_{i}(1,9);
                 faultgrid_wgs84{i}(jj,6) = 1;
            end
              
        end %end kk loop
      
   end %end if/else loop
  
    end %end jj loop
clear dip
end %end i loop

%% Create figure of intersecting faults
%Figure created out of faultintersection function
%Figure will not generate if faults don't intersect
%Figure currently set for intersection of Mlungusi and Liwawadzi faults with Chingale Step

flt_1 = 314; flt_2 = [315 316];

flt2a={}; cutoffpoints={};
fig_option=0; col={'c','r'};
% Finds area of flt2
flt_indx1=find(num_fault(:,1)==flt_1); %indx of fault geometry not changed


figure(408);
for i=1:length(flt_2)
    
    %Run function, script in 'misc functions'
    flt_indx2(i)=find(num_fault(:,1)==flt_2(i)); %indx of fault that is cut off by flt 1
   [area2,flt2a{i},newpoints]=faultintersect(flt_indx1,flt_indx2(i),geometry_top,geometry_bottom,fault_geometry_points,faultgrid_wgs84,fault_geom_,num_fault,fig_option);

   hold on
        patch('XData',fault_geometry_points{flt_indx1}(:,2),'YData',fault_geometry_points{flt_indx1}(:,1),'ZData',fault_geometry_points{flt_indx1}(:,3),'FaceColor','g');hold on
        plot3(fault_geometry_points{flt_indx1}(:,2),fault_geometry_points{flt_indx1}(:,1),fault_geometry_points{flt_indx1}(:,3),'k-'); 
        plot3(fault_geometry_points{flt_indx2(i)}(:,2),fault_geometry_points{flt_indx2(i)}(:,1),fault_geometry_points{flt_indx2(i)}(:,3),'k-');   
        patch('XData', newpoints(:,2),'YData', newpoints(:,1),'ZData', newpoints(:,3),'FaceColor',col{i});hold on 
    ax=gca; ax.ZDir = 'reverse'; daspect([0.01 0.01 1]); view(3); view(125,50); grid on %
    xlabel('Longitude'); ylabel('Latitude'); zlabel('Depth (km)');
    
   set(get(gca,'ylabel'),'rotation',330); set(get(gca,'xlabel'),'rotation',50);
    
    end

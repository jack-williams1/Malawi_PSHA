%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   %% SIMPLYFING THE MSSD FAULT GEOMETRY-2   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FAULT GEOMETRY DEFINED BY GRID WITHIN FAULT POLYGON %%%%

%%% FOR FAULT WIDTHS LIMITED BY 35 KM CRUST%%%%%

%%%% SEE FAULT GEOM_1 FOR WIDTHS LIMITED BY LEONARD (2010)  %%%%%%

%%%% THIS REPACES THE 'faultsource_SMAFD.m' FILE THAT PRODUCED FAULT GRIDS
%%%% FROM DIRECT EXTRAPOLATION OF SURFACE FAULT GEOMETRY %%%%%

%%% THIS FILE AND ITS OUTPUTS CAN BE FOUND IN ARCHIVE FOLDER %%%%

close all

%Load functions needed for fault gemoetry calculations
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/misc_functions']);
addpath([mydir(1:idcs(end)-1) '/gis_files']);

load map_data_EastAfrica
load MSSD_sources

ellipsoid = almanac('earth','wgs84','meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);

s_th = 35; % Thickness of seismogenic crust

%% Group sec info by fault and create polygons

fault_end=zeros(length(num_fault),11);

for i=1:length(num_fault)
   
    if num_fault(i,13)>0 %index fault by sections 
        
    sec_index = find(num_fault(i,1)==num_sec(:,15));
    
    %index sec 1) lat 2)long 3)length 4) sec strike 5) dip 6) h_level
    %7) h_top 8)% fault_id  %8) sec_id 9) 10) fault strike 11) width = s_th
    
    tmp_geom_{i} =  vertcat(num_sec(sec_index,[6:10,12:13,15,1,16]));
    
    %append width based on thickness of seismogenic crust
    fault_geom_{i} = horzcat(tmp_geom_{i},s_th/sin(tmp_geom_{i}(1,5)*pi/180)*ones(height(tmp_geom_{i}),1));
        
    clear sec_index tmp_geom_
    
    else %single section fault geom info taken from num_fault
   
    tmp_geom_{i} =  vertcat(num_fault(i,[4:8,10:11,1,1,7]));   
    fault_geom_{i} = horzcat(tmp_geom_{i},s_th/sin(tmp_geom_{i}(1,5)*pi/180)*ones(height(tmp_geom_{i}),1));
    
    clear tmp_geom_
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

%% Find if two faults intersect and correct geometry
% If they do intersect, removes grid points of shorter faults
% Flt 1 is fault geometry maintaned; Flt 2 is cut off fault
% Faults to correct:

flt_num = [ 308 303;... % Elephant Marsh (308) & Panga-1 (303)
            308 304;... % Elephant Marsh (308) & Panga-2 (304)
            302 305;... % Mwanza(302) & Panga-3 (305)  
            302 306;... % Mwanza(302) % Panga-4(306)
            314 315;... % Chingale Step (314) & Liwawadzi (315)
            314 316;... % Chingale Step (314) & Mlungusi (316)
            314 317;... % Chingale Step (314) & Mtsimukwe (317)
            314 318;... % Chingale Step (314) & Lisungwe-1 (318)
            314 319;... % Chingale Step (314) & Lisungwe-2 (319)
            314 320;... % Chingale Step (314) & Wamkurumadzi-1 (320)
            314 321;... % Chingale Step (314) & Wamkurumadzi-2 (321)
            322 325;... % Tsikulamowa (322) & Bwanwie (325)
            333 323;... % Malombe (333) & Liwonde National Park (323)
            333 324;... % Malombe (333) & Chiloli (324)
            334 333;... % Makanjira (334) & Malombe (333)
            355 329;... % South Basin 13a (355) & BMF-1 (329)
            355 330;... % South Basin 13a (355) & BMF-2 (330)
            330 326;... % BMF-2 (330) & Bwanjwe (326)
            328 355;... % Chirobwe-Ncheu-1 (328) & South Basin 13a (355)
            328 356;... % Chirobwe-Ncheu-1 (328) & South Basin 13b (356)
            337 335;... % Metangula-1 (337) & Cassimo (335)
            355 364;... % South Basin 13a (355) & South Basin 10 (364)
            358 363;... % South Basin 6 (358) & South Basin 8 (363)
            358 351;... % South Basin 6 (358) & South Basin 3 (351)
            352 349;... % South Basin 5a (352) & South Basin 11a (349)
            352 350;... % South Basin 5a (352) & South Basin 11b (350)
            361 339;... % South Basin 12 (361) & Liwaladzi (339)
            361 344;... % South Basin 12 (361) & Phirilanyama-1 (344)
            361 345;... % South Basin 12 (361) & Phirilanyama-2 (345)
            359 366;... % South Basin 7a (359) & South Basin 15 (366)
            359 358;... % South Basin 7a (359) & South Basin 6 (358)
            369 384;... % Usisya Main (369) & Central Basin 19 (384)
            384 377;... % Central Basin 19 (384) & Central Basin 2 (377)
            384 385;... % Central Basin 19 (384) & Central Basin 20 (385)
            384 383;... % Central Basin 19 (384) & Central Basin 11 (383)
            369 386;... % Usisya Main (369) & Lipichili  (386)
            370 386;... % Usisya North (370) & Lipichili  (386)
            370 387;... % Usisya North (370) & Lipichili North  (387)
            369 388;... % Usisya Main (369) & Lipichili South  (388)
            406 372;... % Livingstone (406) & Usisya Tip-2 (372) 
            406 373;... % Livingstone (406) & Usisya Tip-3 (373) 
            406 374;... % Livingstone (406) & Usisya Tip-4 (374)
            398 393;... % Mbiri (398) & Wovwe-1 (393)
            398 391;... % Mbiri (398) & South Karonga-1 (391)
            398 392;... % Mbiri (398) & Sabi-1 (392)
            398 396;... % Mbiri (398) & Wovwe-2 (396)
            398 394;... % Mbiri (398) & South Karonga-2 (394)
            398 395;... % Mbiri (398) & Sabi-2 (395)
            407 405;... % St Mary (407) & Karonga-1 (405)
            407 404;... % St Mary (407) & Karonga-2 (404)
            ]; 


flt_indx=zeros(height(flt_num),2);  
fig_option=0; %0 = no figures from fault intersect function, 1 for figures

for ii=1:height(flt_num)

% Finds area of flt2

    flt_indx(ii,1)=find(num_fault(:,1)==flt_num(ii,1)); %indx of fault geometry not changed
    flt_indx(ii,2)=find(num_fault(:,1)==flt_num(ii,2)); %indx of fault that is cut off by flt 1

    %Run function, script in 'misc functions'

   [area2(ii),flt2a,newpoints{ii}]=faultintersect(flt_indx(ii,1),flt_indx(ii,2),geometry_top,geometry_bottom,fault_geometry_points,faultgrid_wgs84,fault_geom_,num_fault,fig_option);
   
   %Replace fault geom
   faultgrid_wgs84{flt_indx(ii,2)}=[];
   faultgrid_wgs84{flt_indx(ii,2)}=flt2a;
   fault_geometry_points{flt_indx(ii,2)}=[];
   fault_geometry_points{flt_indx(ii,2)}=newpoints{ii};
   clear flt2a
end


%% Find if part of multiple fault system and combine geometry

for m=1:height(num_multi_fault)
    %find faults that are part of multifault system
    multi_fault_index{m} = find(num_multi_fault(m,1)==num_fault(:,23));

    %for multi faults
    if isempty(multi_fault_index{m})==0
     tmp1={};   
        for mm=1:length(multi_fault_index{m})
            tmp1{mm}=faultgrid_wgs84{multi_fault_index{m}(mm)};
        end
            tmp2=tmp1';
            %Combine different faults into multi fault system
            faultgrid_wgs84{multi_fault_index{m}(1)}=cell2mat(tmp2);
            %Append multi_fault ID
            tmp3=ones(length(faultgrid_wgs84{multi_fault_index{m}(1)}),1)*num_multi_fault(m,1);
            faultgrid_wgs84{multi_fault_index{m}(1)}=horzcat(faultgrid_wgs84{multi_fault_index{m}(1)},tmp3);
             clear tmp1 tmp2 tmp3
          % 
            %clear geometry from combined cell arrays
        for mm=2:length(multi_fault_index{m})
            faultgrid_wgs84{multi_fault_index{m}(mm)}=zeros(1,width(faultgrid_wgs84{multi_fault_index{m}(1)}));  
        end
            
            
    end
  end

 for m=1:height(num_fault)
    %find faults that AREN'T part of multifault system
    non_multi_fault_index{m} = find(num_fault(m,23)==num_multi_fault(:,1));   
    
      if isempty(non_multi_fault_index{m})==1
           tmp3=zeros(length(faultgrid_wgs84{m}),1);
           %Append 6th column for these faults where value is 0
           faultgrid_wgs84{m}=horzcat(faultgrid_wgs84{m},tmp3);
           clear tmp3
      end

 end
    
% Create variable to index faultgrid by fault id or multifault id
grid_indx2 = zeros(length(faultgrid_wgs84),1);

 for i=1:length(faultgrid_wgs84)
     
     %Index by multi fault id
     if faultgrid_wgs84{i}(1,7)>0
     grid_indx2(i)=faultgrid_wgs84{i}(1,7);
     
     %Index by fault id
     elseif faultgrid_wgs84{i}(1,4)>0
     grid_indx2(i)=faultgrid_wgs84{i}(1,4);
     
     %Empty cell
     else
     grid_indx2(i)=0;
     end
 end
%% Save fault section geometrical model

tmp=faultgrid_wgs84';
faultgrid_wgs84_num=cell2mat(tmp);

%Turn num array back into cell array, but this time indexed by sec
for i=1:length(num_sec)
  secgrid_wgs84{i}=faultgrid_wgs84_num(find(num_sec(i,1)==faultgrid_wgs84_num(:,5)),:);   
end

secgrid_wgs84_2=secgrid_wgs84;

save('sec_geom_2','secgrid_wgs84_2');

%% Save fault grid cell array and create section grid array

faultgrid_wgs84_2=faultgrid_wgs84; %append 1 to mat file so to differentiate 
%it from other fault_geom file
fault_geometry_points_2=fault_geometry_points;

save('fault_geom_2','faultgrid_wgs84_2','grid_indx2','fault_geometry_points_2');

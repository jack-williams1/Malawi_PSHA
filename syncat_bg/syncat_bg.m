%% Create background syncat for multiple source zones
% Creates location (lat and long) and mag of event
% Number of events for each source zone determined by its G-R relationship

%Events run in one big loop, will need to run for a few days!!!

%PSHA set-up
load new_SZ %with updated source zones for Mozambique

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%load in PSHA options. To change, update excel spreadsheet and then update
%related Matlab variable
load('syncat_PSHA_MSSD_input','NumSimu','bg_mmin','PSHA_Zone');

addpath([mydir(1:idcs(end)-1) '/gis_files']); addpath([mydir(1:idcs(end)-1) '/misc_functions']);

load map_data_EastAfrica

load Malawi_GDEM2 %A list of altitudes for a grid of points covering Malawi

ellipsoid = almanac('earth','wgs84','meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);

%% Create background syn cat for multiple source zones
% Creates location (lat and long) and mag of event
% Number of events for each source zone determined by its G-R relationship

malawi_bg_sz = [new_SZ(7) new_SZ(8) new_SZ(9) new_SZ(14) new_SZ(16) new_SZ(17) new_SZ(19) new_SZ(20) new_SZ(23) new_SZ(26) new_SZ(27)];

%Convert each source zone into a polygon for event location selection
%Manaully update with source zones used
SZ_poly7=polyshape(new_SZ(7).ZONE(:,1),(new_SZ(7).ZONE(:,2))); %Lake Victoria
SZ_poly8=polyshape(new_SZ(8).ZONE(:,1),(new_SZ(8).ZONE(:,2))); %Western Rift-Tanganyika
SZ_poly9=polyshape(new_SZ(9).ZONE(:,1),(new_SZ(9).ZONE(:,2))); % Malawi Rift
%Polygon 14 needs hole in as area also includes polygon 15
SZ_poly14=polyshape({new_SZ(14).ZONE(:,1),new_SZ(15).ZONE(:,1)},{new_SZ(14).ZONE(:,2),new_SZ(15).ZONE(:,2)}); %Mweru-1
%Polygon 16 needs hole in as area also includes polygon 17
SZ_poly16=polyshape({new_SZ(16).ZONE(:,1),new_SZ(17).ZONE(:,1)},{new_SZ(16).ZONE(:,2),new_SZ(17).ZONE(:,2)}); %Kariba-Okavango-1
SZ_poly17=polyshape(new_SZ(17).ZONE(:,1),(new_SZ(17).ZONE(:,2))); %Kariba-Okavango-2
SZ_poly19=polyshape(new_SZ(19).ZONE(:,1),(new_SZ(19).ZONE(:,2))); %Eastern Rift
SZ_poly20=polyshape(new_SZ(20).ZONE(:,1),(new_SZ(20).ZONE(:,2))); %Davie Rift
SZ_poly23=polyshape(new_SZ(23).ZONE(:,1),(new_SZ(23).ZONE(:,2))); %Rovuma Basin
SZ_poly26=polyshape(new_SZ(26).ZONE(:,1),(new_SZ(26).ZONE(:,2))); %Nyanga
SZ_poly27=polyshape(new_SZ(27).ZONE(:,1),(new_SZ(27).ZONE(:,2))); %Northeast Mozambique

%Combine polygons
SZ_poly=[SZ_poly7; SZ_poly8; SZ_poly9; SZ_poly14 ;SZ_poly16; SZ_poly17; SZ_poly19; SZ_poly20; SZ_poly23; SZ_poly26; SZ_poly27];

%Distance in km of largest event from Region 1 to consider in source zone
max_dist = 200;

%If use 300 km as max dist need to include SZ15 (Mweru-2) SZ22 (Mozambique), 
%SZ_poly15=polyshape(new_SZ(15).ZONE(:,1),(new_SZ(15).ZONE(:,2))); %Mweru-2
%SZ_poly22=polyshape(new_SZ(22).ZONE(:,1),(new_SZ(22).ZONE(:,2))); %Mozambique

%% CREATE RECORD

for i=1:length(malawi_bg_sz)

    bg_a    = malawi_bg_sz(i).a;  
    bg_b    = malawi_bg_sz(i).b;
    bg_mmax    = malawi_bg_sz(i).Mwmax;
    
    % Annual number of events between Mmin and Mmax within the background zone 
    Npa_bg = 10^(bg_a-bg_b*bg_mmin) - 10^(bg_a-bg_b*bg_mmax);
    
    %Number of events each year in the simulation of length NumSimu
    % Follows Poisson distribution, where the average is Npa_bg
    Num_event_bg = poissrnd(Npa_bg,NumSimu,1);
       
    source_count_bg = sum(Num_event_bg);%Number of events 
    
    count_gm_bg = 0;
    tmp_count_bg = 1;
    %Randomly distribute numbers between 0-1 that is length of number of events
    t_bg   = rand(source_count_bg,1); %Note stored in cell array
  
    %Derives mag of each event, such that magnitudes have a log normal
    %distribution between Mmax and Mmin. As rand function -> 0, so mags approach mmax
    %For every order of magnitude increase in rand, mg increases by 1/b
   
    mag_bg   = - log10(10^(-bg_b*bg_mmax)+rand(source_count_bg,1)*(10^(-bg_b*bg_mmin)-10^(-bg_b*bg_mmax)))/bg_b;
    
    %Randomly distribute events within each source zone polygon

    numEventsIn=source_count_bg;
    dist=zeros(source_count_bg,1); bg_source{i}=zeros(source_count_bg,6);
    lat=zeros(source_count_bg,1); long=zeros(source_count_bg,1);
  
      for j =1:numEventsIn
            flagIsIn = 0; %Ignore events outside polygon
    
            while ~flagIsIn
                lat(j) =min(SZ_poly(i).Vertices(:,2))+(max(SZ_poly(i).Vertices(:,2)) - min(SZ_poly(i).Vertices(:,2)))*rand(1); %randomly distribute lat
                long(j)  =min(SZ_poly(i).Vertices(:,1))+(max(SZ_poly(i).Vertices(:,1)) - min(SZ_poly(i).Vertices(:,1)))*rand(1); %randomly distribute long
                flagIsIn = inpolygon(long(j),lat(j),SZ_poly(i).Vertices(:,1),SZ_poly(i).Vertices(:,2));  
            end
            
            %Determine distance from event to PSHA area
            %Distance is negative if event within PSHA area
                dist(j)=deg2km(p_poly_dist(long(j),lat(j),PSHA_Zone(:,2),PSHA_Zone(:,1)));
            %Create cell array, one for each source zone
            %Rows are
            %1) distance from PHSA source zones
            %2) long,
            %3) lat
            %4) mag
            %5) simulation cycle
            %6) source zone
               
             bg_source{i}(j,:) =horzcat(dist(j),long(j),lat(j),mag_bg(j),t_bg(j), i); %
            
             count_gm_bg = count_gm_bg + 1;
            
            if count_gm_bg == tmp_count_bg*20000                
               tmp_count_bg = tmp_count_bg + 1;
            disp(['Current bg event number is : ',num2str(count_gm_bg),' out of ',num2str(source_count_bg), ' for source ',num2str(i)]);
             end
        
      end %close j loop
   
      %remove events >max dist from assessed region
      rem_events=find(bg_source{i}(:,1)>max_dist);
      bg_source{i}(rem_events,:)=[];
      
        clear count_gm_bg tmp_count_bg
end


%% clean and save catalog

%Convert seperate catalog from each source zone into one array
bg_catalog=cell2mat(bg_source');
save('bg_catalog','bg_source');

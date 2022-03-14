% RANDOMLY FLOAT RUPTURES ON FAULT PLANE

% Associate ruptures with grid points from faultgrid_wgs84
% Ruptures from Adapted MSSD event catalog 
% Rupture position constrain by fault length and section width
% For events with Mw>5.4

function [gp]= source2site_v3(syncat_ijk,faultgrid_wgs84_1,faultgrid_wgs84_2,num_fault,num_multi_fault,num_MSSD,grid_indx1,grid_indx2)
 
clear gp gp_tmp gp_tmp1

rigidity = 33*10^9; 

%emperical constants from Leonard 2010
c1=17.5; c2 = 3.8*10^-5; a=3;  b=6.1;

    %clear gp_tmp gp_tmp1 source_id m_flt_indx flt_all_indx flt_indx mssd_indx dip
    %clear l w r_area z grd_indx source_pts source_geom bound bound_top l_yy l_yy_deg 
    %clear nst rst rst_sp M I gp_check

    %Derive indexes for num_fault
    source_id = syncat_ijk(4);

    if source_id>600
    m_flt_indx = find(num_multi_fault(:,1)==source_id);%for indexing num_multi_fault
    flt_all_indx = find(source_id==num_fault(:,23));
        
    dip=num_fault(flt_all_indx(1),8);  ttd = num_fault(flt_all_indx(1),11); %toptop depth
    
    else
    flt_indx = find(num_fault(:,1)==source_id);%for indexing num_fault
    dip = num_fault(flt_indx,8); ttd = num_fault(flt_indx,11);
    end
    
    mssd_indx =  find(num_MSSD(:,1)==source_id);%for indexing MSSD
   
   %Define rupture geomety using inverse of Leonard 2010 equations
   l=10^((log10(10^(1.5*syncat_ijk(3)+9.05))-1.5*log10(c1)-log10(c2*rigidity))*2/5)/1000;
   w=(c1*(l*1000)^0.66)/1000;%Leonard (2010) equation 6
   r_area = 10^(2/3*(1.5*syncat_ijk(3)+9.05-log10(c2*rigidity)))/10^6;%Leonard (2010) equation 7

   z=w*sin(dip*pi/180);%Depth extent of rupture
   
   %Index source that rupture is associated with
    if syncat_ijk(5)==1  %Partial fault width, use fault_geom_1 model
        grd_indx = find(syncat_ijk(4)==grid_indx1(:));
        source_pts = faultgrid_wgs84_1{grd_indx};
        s_area =  num_MSSD(mssd_indx,9);

    else %Full crust fault width, use fault_geom_2 model
        grd_indx = find(syncat_ijk(4)==grid_indx2(:));
        source_pts = faultgrid_wgs84_2{grd_indx};
        s_area =  num_MSSD(mssd_indx,11);
    end
     
    indx=1;

    while length(indx)<4 %floating depth must result in indexing at least 4 fault points 
    
     %bottom top depth. Difference between fault and rupture depth
     btd = max(source_pts(:,3)-z);
        
     %float top depth randomly between ttd and btd
     ftd = (btd-ttd)*rand(1,1);
        
     %New bottom depth of rupture given random assignment of ftd
     std_r = ftd+z; clear indx
     
     %Index cells in fault grid between ftd and std
     indx=find(source_pts(:,3) > ftd & source_pts(:,3) <std_r);
    end
    
    gp_tmp(:,:) = source_pts(indx,1:4);
    
    %gp(:,3) = gp(:,3) + EQCAT(ijk,6); %Still vertically projecting, but can change
   
    %if  magnitude is close to or bigger than Mmax due to +0.15 variation in 
    %Stochastic Event Catalog method, take entire source geometry
    
        if (syncat_ijk(3)+0.05)>=num_MSSD(mssd_indx,12) && syncat_ijk(5)==1
            gp=source_pts(:,1:3);  gp_check = size(source_pts,1)*(1/cos(dip*pi/180))/r_area;
            return %function can end here
            
        elseif (syncat_ijk(3)>+0.05)>=num_MSSD(mssd_indx,12) && syncat_ijk(5)==2
            gp=gp_tmp(:,1:3); gp_check = size(gp_tmp,1)*(1/cos(dip*pi/180))/r_area;
            return  %function can end here
            
        end
      
 
    %Laterally vary position of rupture tip
    
    if source_id>600 %Find boundary around each fault polygon
        tmp_id = unique(gp_tmp(:,4));
        bound_geom_tmp={};bound_geom_top_tmp={};
       for pp=1:length(unique(gp_tmp(:,4)))
           tmp_indx = find(tmp_id(pp)==gp_tmp(:,4));
           tmp_geom = gp_tmp(tmp_indx,:);%Index all grid points for each fault
           bound_tmp=boundary(tmp_geom(:,2),tmp_geom(:,1)); %Find indicies of fault boundary
           bound_geom_tmp{pp} = [tmp_geom(bound_tmp,1),tmp_geom(bound_tmp,2)];%Cooridinates of each point around fault
           bound_top_tmp=find(tmp_geom(bound_tmp,3)<(0.7+min(tmp_geom(bound_tmp,3)))); %Find top of polygon
           bound_geom_top_tmp{pp} = [tmp_geom(bound_tmp(bound_top_tmp),1),tmp_geom(bound_tmp(bound_top_tmp),2)];%Cooridinates of each point at top of fault
       
       end
       bound_geom = cell2mat(bound_geom_tmp');  bound_top_geom = cell2mat(bound_geom_top_tmp');
       %find indicies of fault boundaries in gp_tmp
       [~, bound] = ismember(bound_geom,gp_tmp(:,1:2)); bound=bound(:,1);
       [~, bound_top] = ismember(bound_top_geom,gp_tmp(bound,1:2)); bound_top=bound_top(:,1);
       
    else %for single fault system, single boundary to define
     bound=boundary(gp_tmp(:,2),gp_tmp(:,1));
     bound_top=find(gp_tmp(bound,3)<(0.7+min(gp_tmp(bound,3)))); %Find top of polygon
    end
   
    l_yy=abs(l*cos(num_MSSD(mssd_indx,4)*pi/180)); %N-S extent of rupture (in km)
    l_yy_deg=km2deg(l_yy);%N-S extent of rupture (in degrees)
    nst = max(gp_tmp(bound(bound_top),1))-l_yy_deg;%northern most south tip of rupture
       
    %float southern rupture tip somewhere between southern tip of fault and nst
    rst = ((nst-min(gp_tmp(bound(bound_top),1)))*rand(1,1))+min(gp_tmp(bound(bound_top),1));
    
    % find closest point in grid to rupture southern tip
    [M, I]=min(abs(rst-gp_tmp(bound(bound_top),1)),[],1);
    
    % find coordinate and fault for rupture southern tip
    rst_sp = [gp_tmp(bound(bound_top(I)),1),gp_tmp(bound(bound_top(I)),2),gp_tmp(bound(bound_top(I)),4)];
    
    try %first tries to float rupture using dimensions from Leonard 2010 eqs.
    
    if source_id>600
        
        flt1_indx_gp=find(gp_tmp(bound(bound_top),4)==rst_sp(3));%index fault at top of geom
        [tmp, fnt_indx] = max(gp_tmp(bound(bound_top(flt1_indx_gp)),1));
        fnt = gp_tmp(bound(bound_top(flt1_indx_gp(fnt_indx))),[1:2]);%find northern most point fault at top of geom
        l1_yy = deg2km(fnt(:,1)-rst_sp(1));%N-S extent of fault that has ruptured
        flt_indx = find(num_fault(:,1)==rst_sp(3));%index fault in num_fault
        
    end
    
    % for single fault sources, or ruptures short enough to rupture on one fault
    if source_id<600 || l_yy<l1_yy
    
    %Use coordinate to make polygon around rupture with
    %defined length and width, and fault strike
    
    %Aziumth is +180 if fault RHR strike is south bearing
    if num_fault(flt_indx,7) >90 && num_fault(flt_indx,7)<270
    fnt_sp = reckon(rst_sp(1),rst_sp(2),km2deg(l),num_fault(flt_indx,7)+180);
    else
    fnt_sp = reckon(rst_sp(1),rst_sp(2),km2deg(l),num_fault(flt_indx,7));    
    end
    
    %Find bottom parts of rupture polygon using fault dip_dir
    fst_bp = reckon(rst_sp(1),rst_sp(2),km2deg(w*cos(dip*pi/180)),num_fault(flt_indx,7)+90);
    fnt_bp = reckon(fnt_sp(1),fnt_sp(2),km2deg(w*cos(dip*pi/180)),num_fault(flt_indx,7)+90);
    rupture_poly = polyshape([rst_sp(2) fnt_sp(2) fnt_bp(2) fst_bp(2)],[rst_sp(1) fnt_sp(1) fnt_bp(1) fst_bp(1)]);
    %Find and index grid points in r
    [in_r,on_r] =inpolygon(gp_tmp(:,2),gp_tmp(:,1),rupture_poly.Vertices(:,1),rupture_poly.Vertices(:,2));
    gp_tmp1=[gp_tmp(in_r,1),gp_tmp(in_r,2),gp_tmp(in_r,3),gp_tmp(in_r,4)];
    gp=gp_tmp1(:,1:3);
    gp_check = size(gp_tmp1,1)*(1/cos(dip*pi/180))/r_area;
    clear in_r on_r fnt_sp fst_sp fst_bp fnt_bp rupture_poly 
    
    else %if multi_fault rupture
    
    %if random rupture tip at northern fault tip, don't float geometry on fault
    if rst_sp(1)~=fnt(1) && rst_sp(2) ~=fnt(2)
    
        %find bottom parts of geometry of fault-1 
        fst_bp = reckon(rst_sp(1),rst_sp(2),km2deg(w*cos(dip*pi/180)),num_fault(flt_indx,7)+90);
        fnt_bp = reckon(fnt(1),fnt(2),km2deg(w*cos(dip*pi/180)),num_fault(flt_indx,7)+90);
        rupture_poly = polyshape([rst_sp(2) fnt(2) fnt_bp(2) fst_bp(2)],[rst_sp(1) fnt(1) fnt_bp(1) fst_bp(1)]);
    
        %find grid points that ruptured in fault-1
        [in_r,on_r] =inpolygon(gp_tmp(:,2),gp_tmp(:,1),rupture_poly.Vertices(:,1),rupture_poly.Vertices(:,2));
        gp_tmp1=[gp_tmp(in_r,1),gp_tmp(in_r,2),gp_tmp(in_r,3),gp_tmp(in_r,4)];
    
        clear in_r on_r fnt_sp fst_sp fst_bp fnt_bp rupture_poly
    
    else    
        gp_tmp1=[];
    end
    
    %remaining length of rupture geometry to define
    l_need = l-l1_yy/abs(cos(num_fault(flt_indx,7)*pi/180));
    clear flt_indx
    
    %index next fault to north
    flt_sp = num_fault(flt_all_indx,[1,4]);
    flt_sp = sortrows(flt_sp,2); 
    
    %manual correction to multifault systems SB 7-12 and CB19-20
    %converging flt dip mean faults not correctly sorted N-S by sp
    if source_id == 611 || source_id == 616 
       flt_sp=flip(flt_sp);   
    end
        
    flt_1_indx=find(rst_sp(3)==flt_sp(:,1));
    flt_n_id = flt_sp(flt_1_indx+1,1);
    flt_n_indx=flt_1_indx+1;
    
    %if length of next fault <remaining length of rupture
    %take full fault geometry. Repeat until length of next fault>remaining
    %length of geometry
  
    while  l_need > num_fault(find(num_fault(:,1)==flt_n_id),6)
        
        %index fault geom from geom_tmp
        flt_gp_indx=find(gp_tmp(:,4)==flt_n_id);
        
        if isempty(gp_tmp1)==1
        gp_tmp1 = [gp_tmp(flt_gp_indx,1),gp_tmp(flt_gp_indx,2),gp_tmp(flt_gp_indx,3),gp_tmp(flt_gp_indx,4)];  
        else
        gp_tmp1 = vertcat(gp_tmp1,[gp_tmp(flt_gp_indx,1),gp_tmp(flt_gp_indx,2),gp_tmp(flt_gp_indx,3),gp_tmp(flt_gp_indx,4)]);
        end
        l_need = l_need-num_fault(find(num_fault(:,1)==flt_n_id),6);
        
        %if reach northern most part of multi_fault system, break
        if max(gp_tmp1(:,1)) ==  max(gp_tmp(:,1))
            gp=gp_tmp1(:,1:3);
            gp_check = size(gp_tmp1,1)*(1/cos(dip*pi/180))/r_area;
            clear flt_gp_indx
            break 
        end
        
        %add indicies by one
        flt_n_id = flt_sp(flt_n_indx+1,1);
        flt_n_indx=flt_n_indx+1; 
        clear flt_gp_indx
    end %end while loop
    
        %index next fault to rupture in num_fault
        flt_indx = find(flt_sp(flt_n_indx,1)==num_fault(:,1));
    
    if l_need<2 %if only small part of fault left to index, don't find next participating fault
        gp=gp_tmp1(:,1:3);
        gp_check = size(gp_tmp1,1)*(1/cos(dip*pi/180))/r_area;
    else
    
    %For remaining l_need, find geometry of northernmost participating fault
        flt1_indx_gp=find(gp_tmp(bound(bound_top),4)==flt_n_id);%index northern most fault at top of geom
        [~, fst_indx] = min(gp_tmp(bound(bound_top(flt1_indx_gp)),1));%index its southern most point
        fst = gp_tmp(bound(bound_top(flt1_indx_gp(fst_indx))),[1:2]);%find northern most point fault at top of geom
    
        %Find northermnmost part of rupture
        %Aziumth is +180 if fault RHR strike is south bearing
        if num_fault(flt_indx,7) >90 && num_fault(flt_indx,7)<270
            fnt_sp = reckon(fst(1),fst(2),km2deg(l_need),num_fault(flt_indx,7)+180);
        else
            fnt_sp = reckon(fst(1),fst(2),km2deg(l_need),num_fault(flt_indx,7));    
        end
    
        fst_bp = reckon(fst(1),fst(2),km2deg(w*cos(dip*pi/180)),num_fault(flt_indx,7)+90);
        fnt_bp = reckon(fnt_sp(1),fnt_sp(2),km2deg(w*cos(dip*pi/180)),num_fault(flt_indx,7)+90);
    
        rupture_poly = polyshape([fst(2) fnt_sp(2) fnt_bp(2) fst_bp(2)],[fst(1) fnt_sp(1) fnt_bp(1) fst_bp(1)]);
    
        %find grid points that ruptured in fault-1
        [in_r,on_r] =inpolygon(gp_tmp(:,2),gp_tmp(:,1),rupture_poly.Vertices(:,1),rupture_poly.Vertices(:,2));
        gp_tmp1=vertcat(gp_tmp1,[gp_tmp(in_r,1),gp_tmp(in_r,2),gp_tmp(in_r,3),gp_tmp(in_r,4)]);
    
        clear in_r on_r fnt_sp fst_sp fst_bp fnt_bp rupture_poly 
    
        gp=gp_tmp1(:,1:3);
        gp_check = size(gp_tmp1,1)*(1/cos(dip*pi/180))/r_area;
        
    end %end if l_need<2 statement
    
    end %end fault or multi-fault if statement

    %if floating rupture area results in area less than 67% of area estimated by Leonard 2010 equations
    %force error so that it runs through catch, and area is iteratively defined instead
    assert(gp_check(ijk)>0.67);

    %if cannot float rupture, iteratively define it, first based on gp_tmp,
    %then source_pts starting from gp_tmp1, and then source_pts starting from rst_sp
    catch ME
        
        %if try completely fails to find any points, create dummy gm_tmp1
        %and randomly find new source points
    if  exist('gp_tmp1','var')==0
        gp_tmp1=[];   
        r_indx = randi([1 length(source_pts)],1,1);
        rst_sp=[source_pts(r_indx,1),source_pts(r_indx,2)];
        clear r_indx
    end
    
        if isempty(gp_tmp1)==0 
        %number of grid points needed to match r_area
        %account that grid spacing is 1 km along strike, and 1.6 km down dip
        
        g_need = round((r_area-(length(gp_tmp1)*(1/cos(dip*pi/180))))/(1/cos(dip*pi/180)));
       
        %First try to extend rupture along strike across gp_tmp        
        %find points in gp_tmp not in gp_tmp1, and distance between them
         
        mindiff=zeros(size(gp_tmp,1),2);
          for gg=1:size(gp_tmp,1)
              if ismember(gp_tmp(gg,:),gp_tmp1,'rows')==0
                mindiff(gg,1)=gg;
                mindiff(gg,2) = min((gp_tmp(gg,1)-gp_tmp1(:,1)).^2 + (gp_tmp(gg,2)-gp_tmp1(:,2)).^2);
              end
          end
           
        mindiff=sortrows(mindiff,2);%sort rows into increasing distance
        mindiff(find(mindiff(:,1)==0),:)=[];%remove zeros
        
          %if enough points in mindiff to satisfy g_need take points from gp_2
         if size(mindiff,1)>g_need
            gp_take=mindiff(1:g_need,1);
        
            gp_tmp1=vertcat(gp_tmp1,gp_tmp(gp_take,1:4)); clear g_need
            gp=gp_tmp1(:,1:3);
            gp_check = size(gp_tmp1,1)*(1/cos(dip*pi/180))/r_area;
            
        else %if insufficient points along gmp_2, use source_pts to indx geom_pts 
            mindiff=zeros(size(source_pts,1),2);
        
            for gg=1:size(source_pts,1)
                if ismember(source_pts(gg,1:4),gp_tmp1,'rows')==0
                mindiff(gg,1)=gg;
                mindiff(gg,2) = min((source_pts(gg,1)-gp_tmp1(:,1)).^2 + (source_pts(gg,2)-gp_tmp1(:,2)).^2);
                end
            end
        
            mindiff=sortrows(mindiff,2);%sort rows into increasing distance
        
            mindiff(find(mindiff(:,1)==0),:)=[];%remove zeros
            %take closest points to rst_sp
            if size(mindiff,1)>g_need
                gp_take=mindiff(1:g_need,1);
                gp_tmp1=vertcat(gp_tmp1,source_pts(gp_take,1:4)); clear g_need
                gp=gp_tmp1(:,1:3);
                gp_check = size(gp_tmp1,1)*(1/cos(dip*pi/180))/r_area;
            %if g_need larger than size of avaialble grid points, can assume full fault rupture    
            else
                gp=source_pts;
            end
        end
            
    else %if gp_tmp1 not indexed yet, find closes points to rst_sp in source id
        
        g_need = round(r_area/(1/cos(dip*pi/180))); %number of grid points as function or r_area
        mindiff=zeros(length(source_pts),2);
        for gg=1:length(source_pts)
               mindiff(gg,1)=gg;
               mindiff(gg,2) = (source_pts(gg,1)-rst_sp(1)).^2 + (source_pts(gg,2)-rst_sp(2)).^2;
        end   
              mindiff=sortrows(mindiff,2);%sort rows into increasing distance
              gp_take=mindiff(1:g_need,1);
              gp_tmp1=source_pts(gp_take,1:4);
              gp=gp_tmp1(:,1:3);  
              gp_check= size(gp_tmp1,1)*(1/cos(dip*pi/180))/r_area;
    end
        
    end %end try statement
       
    %gp(:,3) = gp(:,3) + td; %Not vertically projecting, but can change
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTION FOR SOURCE TO SITE USING SYNCAT_MSSD CATALOG %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomly float individual section ruptures in fault plane
% For syncat_MSSD catalog only
% Use source2site_v3 for syncatYC85_MSSD
% Rupture position constrain by fault width and section width

function [gp]= source2site_v1(syncat_ijk,faultgrid_wgs84_1,secgrid_wgs84_1,num_sec,num_fault,grid_indx1);

%if multifault rupture

if syncat_ijk(4)>600
    % Don't float multifault ruptures
     mflt_index = find(syncat_ijk(4)==grid_indx1);
     %Only index the first value, as other ones are empty cells
     gp(:,:) = faultgrid_wgs84_1{mflt_index(1)}(:,1:3);
  
    
elseif syncat_ijk(4)>300 %if fault rupture
        
        flt_indx= find(syncat_ijk(4)==num_fault(:,1));
        
       
        % if fault is not part of a multifault system
        if num_fault(flt_indx,23)==0 
            gp(:,:) = faultgrid_wgs84_1{flt_indx}(:,1:3);
            
       %if fault part of multi fault system, float rupture     
        else

            mflt_id = num_fault(flt_indx,23); %Find multifault id fault is part of
            flt_grid_indx=find(mflt_id==grid_indx1);%Use multi_fault ID to index fault grid
            
            %find grid points in multi-fault system associated with
            %particular fault
            fault_gp(:,:) = faultgrid_wgs84_1{flt_grid_indx}(find(faultgrid_wgs84_1{flt_grid_indx}(:,4)==syncat_ijk(4)),1:3);
            
            ttd = num_fault(flt_indx,11); %toptop depth
            %bottom top depth. Difference between multi fault and fault width
            %projected through fault dip
            btd = max(fault_gp(:,3))-(num_fault(flt_indx,12)*sin(num_fault(flt_indx,8)*pi/180));
            %float top depth randomly between ttd and btd
            ftd = (btd-ttd).*rand(1,1);
            %New bottom depth of rupture given random assignment of ftd
            std_r = ftd+sin((num_fault(flt_indx,8)*pi/180))*num_fault(flt_indx,12);
            indx=find(fault_gp(:,3) > ftd & fault_gp(:,3) <std_r);
            
            %Some faults cannot be extrapolated down dip due to geometrical tapering effects.
            %Otherwise keep with original geometry
            
            if length(indx)>0
                gp(:,:) = fault_gp(indx,1:3);
            else
                gp(:,:) = fault_gp(:,1:3);
            end
            
        end
        
else 
        %find section ruptures
        
        s_id = find(syncat_ijk(4)==num_sec(:,1));%find section num_id so can numerically index in num_sec
        f_id = find(num_sec(s_id,15)==num_fault(:,1));%find fault num_id so can numerically index in num_sec
        
        if num_sec(s_id,18)>0 %if sec is part of multi-fault, index that too
              mflt_id = num_sec(s_id,18); %Find multifault id sec is part of
              flt_grid_indx=find(mflt_id==grid_indx1);
        else
              flt_grid_indx=f_id;

        end
        
        
        %check indexing between num_sec, EQCAT, and secgrid is consistent
        if secgrid_wgs84_1{s_id}(1,5)==syncat_ijk(4);
            ttd = num_sec(s_id,13); %toptop depth
            %bottom top depth. Difference between fault and section width
            %projected through fault dip
            btd = max(faultgrid_wgs84_1{f_id}(:,3))-(num_sec(s_id,14)*sin(num_sec(s_id,10)*pi/180));
            %float top depth randomly between ttd and btd
            ftd = (btd-ttd).*rand(1,1);
            %New bottom depth of rupture given random assignment of ftd
            std_r = ftd+(sin(num_sec(s_id,10)*pi/180)*num_sec(s_id,14));
            %Index cells in sec grid between ftd and std
            indx=find(secgrid_wgs84_1{s_id}(:,3) > ftd & secgrid_wgs84_1{s_id}(:,3) <std_r);
            
            %Some sections cannot be extrapolated down dip due to geometrical tapering effects.
            %Otherwise keep with original geometry
            
            if length(indx)>0
                gp(:,:) = secgrid_wgs84_1{s_id}(indx,1:3);
            else
                gp(:,:) = secgrid_wgs84_1{s_id}(:,1:3); 
            end
                       
        else
            %error in indexing
            print('error in indexing')
        end 
end

%gp(:,3) = gp(:,3) + td; %Not vertically projecting, but can change
%td = syncat_ijk(7);

end


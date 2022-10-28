%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find intersection of two faults and remove points %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flt1 is relatively long fault whose geometry is not chanved
% flt2 is cut off by intersection with fault 1 
% Returns grid points of new fault in flt2a

% option to return vert as varargout for use in faultintersect_figure

function [area2,flt2a,newpoints_f] = faultintersect(flt1,flt2,geometry_top,geometry_bottom,fault_geometry_points,faultgrid_wgs84,fault_geom_,num_fault,fig_option);

clear dt1_ s1_ s1_faces s1_vertices dt2_ dt2_utm s2_ s2_faces s2_vertices ss1 ss2 tmp1 tmp3 tmp4 geometry_top_utm geometry_bottom_utm area2
clear vert in on remove_points remove_points_mat cutoffpoints cutoffpoints_utm cutoff_alphaShape flt2a flt2area flt2surfacearea id_sec vert_tmp

%Find geometry of Flt1 in terms of triangles
for i=1:height(geometry_top{flt1})-1   
    %Find geometry of each distinct section in terms of delunay triangles
         dt1_{i}=delaunayTriangulation(vertcat(geometry_top{flt1}(i:(i+1),1),geometry_bottom{flt1}(i:(i+1),1)),vertcat(geometry_top{flt1}(i:(i+1),2),geometry_bottom{flt1}(i:(i+1),2)),vertcat(geometry_top{flt1}(i:(i+1),3),geometry_bottom{flt1}(i:(i+1),3)));
         if isempty(dt1_{i}.ConnectivityList)==0 %Remove geometry from short link sections
         s1_{i}=tetramesh(dt1_{i},'Visible','off');
         s1_vertices{i}=s1_{i}.Vertices;
         s1_faces{i}=s1_{i}.Faces;
           
       end
end


%Repeat for Flt2
for i=1:height(geometry_top{flt2})-1
    dt2_{i}=delaunayTriangulation(vertcat(geometry_top{flt2}(i:(i+1),1),geometry_bottom{flt2}(i:(i+1),1)),...
        vertcat(geometry_top{flt2}(i:(i+1),2),geometry_bottom{flt2}(i:(i+1),2)),vertcat(geometry_top{flt2}(i:(i+1),3),geometry_bottom{flt2}(i:(i+1),3)));
    if isempty(dt2_{i}.ConnectivityList)==0
       s2_{i}=tetramesh(dt2_{i},'Visible','off');
       s2_vertices{i}=s2_{i}.Vertices;
       s2_faces{i}=s2_{i}.Faces;
    end
end

%find area of Flt2 using alpha shape
% Need to convert to utm to find area in m^2
[geometry_top_utm(:,2),geometry_top_utm(:,1)]=deg2utm(geometry_top{flt2}(:,1),geometry_top{flt2}(:,2));
[geometry_bottom_utm(:,2),geometry_bottom_utm(:,1)]=deg2utm(geometry_bottom{flt2}(:,1),geometry_bottom{flt2}(:,2));
for i=1:height(geometry_top{flt2})-1
        flt2shape{i}=alphaShape(vertcat(geometry_top_utm(i:(i+1),2),geometry_bottom_utm(i:(i+1),2)),vertcat(geometry_top_utm(i:(i+1),1),geometry_bottom_utm(i:(i+1),1)),...
        vertcat(10^3*(geometry_top{flt2}(i:(i+1),3)),10^3*(geometry_bottom{flt2}(i:(i+1),3))));
        flt2surfacearea{i}=surfaceArea(flt2shape{i})/(2*10^6);%divide by two as only considering plane area
        %divide by 10^6 to get to km^2
end 
%Sum area from individual sections to get total fault area
flt2area=sum(cell2mat(flt2surfacearea));


% Find intersection of two fault surfaces

vert_tmp={}; 
vert={}; tmp_count1=0;
for i=1:height(geometry_top{flt2})-1

    tmp_count=0;
    
    if isempty(s2_faces{i}) ==0 
    ss2.faces=s2_faces{i};
    ss2.vertices=s2_vertices{i};
 for j=1:height(geometry_top{flt1})-1
      if isempty(s1_faces{j}) ==0 
        ss1.faces=s1_faces{j};
        ss1.vertices=s1_vertices{j};
        
        [a,SurfAB]= SurfaceIntersection(ss1,ss2);
            
        if isempty(SurfAB.vertices)==0
            svertices_=SurfAB.vertices;
            vert_tmp{i}{j}=svertices_;
        else
            % Put in dummy values for surfaces that don't intersect and 
            % can be subsequently removed
            vert_tmp{i}{j}=[0, 0, 0];
            tmp_count=tmp_count+1;
        end
            
        clear ss1 a SurfAB  svertices_
      end
      
  
  end %end j loop
      if tmp_count==length(j)
         tmp_count1=1;%record sections with no intersections
      end
    clear ss2
 end %end if statement for s2_faces{i}
 
 %clean intersection vector
    vert{i}=cell2mat(vert_tmp{i}');
    indx_clr=find(vert{i}(:,1)==0);
    vert{i}(indx_clr,:)=[];
    vert{i} = unique(vert{i},'rows');
end



%Plot fault planes with line of intersection (vert)  

if fig_option==1

    figure(100);
    patch('XData',fault_geometry_points{flt1}(:,2),'YData',fault_geometry_points{flt1}(:,1),'ZData',fault_geometry_points{flt1}(:,3),'FaceColor','g');hold on
    patch('XData',fault_geometry_points{flt2}(:,2),'YData',fault_geometry_points{flt2}(:,1),'ZData',fault_geometry_points{flt2}(:,3),'FaceColor','c');hold on
    if isempty(vert)==0
        for i=1:length(vert)
            if isempty(vert{i})==0
            plot3(vert{i}(:,2),vert{i}(:,1),vert{i}(:,3),'r-','LineWidth',4) 
            end
        end
    end; hold on
    ax=gca; ax.ZDir = 'reverse'; daspect([0.01 0.01 1]); view(3); grid on; hold off %
    xlabel('Longitude'); ylabel('Latitude'); zlabel('Depth (km)');
end

% If no intersection, function ends here

if isempty(vert)==1
 flt2a = "These faults do not intersect";
 area2 = "no need for this measurement";
 disp(flt2a)
else    
    
% Define polygon cutoff that is part of fault plane from fl2 to be removed 
% Define revised polygon for flt2

cutoffpoints={}; geometry_bottom_tmp={};
ellipsoid = almanac('earth','wgs84','meters');

if fault_geom_{flt2}(1,10) >90 && fault_geom_{flt2}(1,10) <270
               
    for mm=1:length(vert)
        
        if isempty(vert{mm})==0
            
                    [~,mmindx]=max(vert{mm}(:,1));
                    tmp_distance=geometry_bottom{flt2}(mm+1,1)-vert{mm}(mmindx,1);
            
             if mm>1 & tmp_count1==1 & isempty(vert{mm-1})==1%if adjoining section has no intersection
                geometry_bottom_tmp{mm}=vertcat(vert{mm},geometry_bottom{flt2}(mm-1,:));
                
              %if only 1 section considered and its northern most point is cutoff, add original coordiante back in        
             elseif tmp_count1==0 & length(vert)==1 & max(vert{mm}(:,1))<geometry_bottom{flt2}(mm+1,1) & tmp_distance>0.1
                 geometry_bottom_tmp{mm}=vertcat(vert{mm},geometry_bottom{flt2}(mm+1,:)); 
            else
                geometry_bottom_tmp{mm}=vert{mm};
             end
           cutoffpoints{mm}=[vert{mm}; flipud(geometry_bottom{flt2}(mm:(mm+1),:)); vert{mm}(1,:)];
        else
            if mm==length(vert)%if at northernmost section tip 
                geometry_bottom_tmp{mm}=geometry_bottom{flt2}(mm:(mm+1),:);
            else
                geometry_bottom_tmp{mm}=geometry_bottom{flt2}(mm,:);
            end
        cutoffpoints{mm}=[];%remove points for flt2 sections that don't intersect with flt1
        end
    end 
    
      
    
    geometry_bottom_mat=cell2mat(geometry_bottom_tmp');
    geometry_bottom_mat=sortrows(geometry_bottom_mat,1,'descend');
    
else
   
    for mm=1:length(vert) 
        if isempty(vert{mm})==0
            
              [~,mmindx]=max(vert{mm}(:,1));
              tmp_distance=geometry_bottom{flt2}(mm+1,1)-vert{mm}(mmindx,1);
            
            if mm ~=length(vert) & tmp_count1==1 & isempty(vert{mm+1})==1%if adjoining section has no intersection
                geometry_bottom_tmp{mm}=vertcat(vert{mm},geometry_bottom{flt2}(mm+1,:));
                
            %if only 1 section considered and its southern most point is cutoff, add original coordiante back in
            elseif tmp_count1==0 & length(vert)==1 & min(vert{mm}(:,1))>min(geometry_bottom{flt2}(:,1)) & tmp_distance>0.1
                geometry_bottom_tmp{mm}=vertcat(vert{mm},geometry_bottom{flt2}(mm+1,:));  
             
                
            else
                geometry_bottom_tmp{mm}=vert{mm};
            
            end
            cutoffpoints{mm}=[vert{mm}; geometry_bottom{flt2}(mm:(mm+1),:); vert{mm}(1,:)];
            
        else
            if mm==1 || mm==length(vert)%if at southernmost section tip
                geometry_bottom_tmp{mm}=geometry_bottom{flt2}(mm:(mm+1),:);
               
            else
                geometry_bottom_tmp{mm}=geometry_bottom{flt2}(mm,:);  
            end
        cutoffpoints{mm}=[];%remove points for flt2 sections that don't intersect with flt1
        end
    end 
    
    
    
    geometry_bottom_mat=cell2mat(geometry_bottom_tmp');
    geometry_bottom_mat=sortrows(geometry_bottom_mat,1,'ascend');
end
       %Correction for Lipilli fault which gets cut by Usisya Twice
    if flt2==86 && flt1==70
        for mm=1:length(vert)
            if isempty(vert{mm})==0
                geometry_bottom_tmp{mm}=geometry_bottom_tmp{mm}(1:length(geometry_bottom_tmp{mm})-1,:);
            else
                tmp=find(fault_geometry_points{flt2}(:,1) > geometry_bottom{flt2}(mm,1) & fault_geometry_points{flt2}(:,1) < geometry_bottom{flt2}(mm+1,1) & fault_geometry_points{flt2}(:,3)>2 ); 
                geometry_bottom_tmp{mm}=fault_geometry_points{flt2}(tmp,1:3);
            end
        end
        geometry_bottom_mat=cell2mat(geometry_bottom_tmp');
        geometry_bottom_mat=sortrows(geometry_bottom_mat,1,'descend');
    end



%create new geometry bottom variable with updated points

new_points=vertcat(geometry_top{flt2},geometry_bottom_mat,geometry_top{flt2}(1,:));

cutoffpoints(cellfun('isempty',cutoffpoints)) = [];
newpoints_f=horzcat(new_points,ones(length(new_points),1)*num_fault(flt2,1));%Returns fault geometry points

% Find points from wgs84 grid to be removed from cutoffpoints using inpolygon
in={};
on={};
remove_points={};

for ii=1:length(cutoffpoints)
    [in{ii},on{ii}]=inpolygon(faultgrid_wgs84{flt2}(:,2),faultgrid_wgs84{flt2}(:,1),cutoffpoints{ii}(:,2),cutoffpoints{ii}(:,1));
    remove_points{ii}=[faultgrid_wgs84{flt2}(in{ii},2) faultgrid_wgs84{flt2}(in{ii},1) faultgrid_wgs84{flt2}(in{ii},3)];
end

remove_points_mat=cell2mat(remove_points');
tmp=ismember(faultgrid_wgs84{flt2},remove_points_mat);
tmp1=find(tmp(:,1)==0);%Find points to be removed from fault grid
flt2a=faultgrid_wgs84{flt2}(tmp1,:);%Extract only points to be kept  

%new area of fault (as a fraction of grid points removed)
area2=flt2area*(length(flt2a)/length(faultgrid_wgs84{flt2}));
    
    if fig_option==1
    % Plot fault planes with points to be kept
    figure(508);
        patch('XData',fault_geometry_points{flt1}(:,2),'YData',fault_geometry_points{flt1}(:,1),'ZData',fault_geometry_points{flt1}(:,3),'FaceColor','g');hold on
        patch('XData',fault_geometry_points{flt2}(:,2),'YData',fault_geometry_points{flt2}(:,1),'ZData',fault_geometry_points{flt2}(:,3),'FaceColor','c');hold on
        for ii=1:length(remove_points)
        patch('XData', cutoffpoints{ii}(:,2),'YData', cutoffpoints{ii}(:,1),'ZData', cutoffpoints{ii}(:,3),'FaceColor','r');hold on    
        plot3(flt2a(:,2),flt2a(:,1),flt2a(:,3),'k*'); hold on;
        end; hold on;
    ax=gca; ax.ZDir = 'reverse'; daspect([0.01 0.01 1]); view(3); grid on %
    xlabel('Longitude'); ylabel('Latitude'); zlabel('Depth (km)');
    end
    
end

end



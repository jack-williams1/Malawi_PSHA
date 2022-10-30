%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Probabilistic Seismic Hazard Analysis for Malawi-1   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plot PSHA figures for site and ground motions
% estimated in PSHA_MSSM_site

clear
close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%add path to GM from Blue Crystal Runs Updated as necessary
addpath([mydir(1:idcs(end-1)-1),'/psha_stored_gm/20221021_GM']);

addpath([mydir(1:idcs(end)-1) '/misc_functions']);

load('syncat_PSHA_MSSM_input','prob_level','GMPEindex','num_GMPEindex','vs30_site_ref',...
    'obs_duration','lower_return_period','num_need','max_bg','par_opts',...
    't_limit','NumSimu','T','num_T','syncat_name','site','site_name');

num_syncat=5; %number of MSSM catalogs

num_site=length(site); num_need_fs=num_need*5; obs_duration_fs=obs_duration*5;

%if code run before, load sorted and clean ground motions and skip straight
%to 'Summary results: Sort data'
load GM_20221021 

%% load and combine results from each parallization

clear  tmp_GM_fs tmp_GM_fs_ref

tmp_GM_fs=cell(1,num_syncat);  tmp_GM_fs_ref=cell(1,num_syncat);
tmp_GM_fs_mean=cell(1,num_site); tmp_GM_fs_mean_ref=cell(1,num_site);

for pp=1:length(par_opts)
    
   disp(['on parellization ',num2str(pp) ' out of ', num2str(length(par_opts))])
    
   load(strcat('AMAX_GM_bg_',string(par_opts(pp)))); load(strcat('AMAX_GM_',string(par_opts(pp))))
   
   for  gg=1:num_GMPEindex
   
        for ii=1:num_site
       
           if pp==1 %create dummy variables
            tmp_GM_bg{gg}{ii}={}; tmp_GM_bg_ref{gg}{ii}={}; 
            tmp_GM_bg{gg}{ii}=single(zeros(1,num_T)); tmp_GM_bg_ref{gg}{ii}=single(zeros(1,num_T));
            
            tmp_GM_fs{gg}{ii}={}; tmp_GM_fs_ref{ii}={}; 
            tmp_GM_fs{gg}{ii}=single(zeros(1,num_T)); tmp_GM_fs_ref{gg}{ii}=single(zeros(1,num_T));
           end
           
           %Note converting values from double to single to stop variables becoming too large
           tmp_GM_bg{gg}{ii}=vertcat(tmp_GM_bg{gg}{ii}, eval(strcat('AMAX_GM_bg_',string(par_opts(pp)),'{',string(num2str(gg)),'}','{',string(num2str(ii)),'}')));
           tmp_GM_bg_ref{gg}{ii}=vertcat(tmp_GM_bg_ref{gg}{ii}, eval(strcat('AMAX_GM_bg_ref_',string(par_opts(pp)),'{',string(num2str(gg)),'}','{',string(num2str(ii)),'}')));
           
           if gg==num_GMPEindex
               
                for kk=1:num_T
                    
                    if pp==1    
                        tmp_GM_mean_bg{ii}{kk}={}; tmp_GM_mean_bg{ii}{kk}=single(zeros(1,4));
                        tmp_GM_mean_bg_ref{ii}{kk}={}; tmp_GM_mean_bg_ref{ii}{kk}=single(zeros(1,4));
                    end
              
                    tmp_GM_mean_bg{ii}{kk}=vertcat(tmp_GM_mean_bg{ii}{kk}, eval(strcat('AMAX_GM_mean_bg_',string(par_opts(pp)),'{',string(num2str(ii)),'}','{',string(num2str(kk)),'}')));
                    tmp_GM_mean_bg_ref{ii}{kk}=vertcat(tmp_GM_mean_bg_ref{ii}{kk}, eval(strcat('AMAX_GM_mean_bg_ref_',string(par_opts(pp)),'{',string(num2str(ii)),'}','{',string(num2str(kk)),'}')));
            end
               
           end
           
    if pp==length(par_opts)
        tmp_GM_bg{gg}{ii}(1,:)=[]; tmp_GM_bg_ref{gg}{ii}(1,:)=[];%Remove dummy row of zeros
    end
    
        end %end ii loop for num_site
   end %end gg loop for num_GMPEindex
   
   for ss=1:num_syncat %for each MSSM syncat
         
       if pp==1
            tmp_GM_fs{ss}=repmat({cell(1,num_site)},num_GMPEindex,1);  tmp_GM_fs_ref{ss}=repmat({cell(1,num_site)},num_GMPEindex,1);
       end
       
       for gg=1:num_GMPEindex
       
        for ii=1:num_site
           
            if pp==1 %create dummy variable
                tmp_GM_fs{ss}{gg}{ii}=single(zeros(1,num_T)); tmp_GM_fs_ref{ss}{gg}{ii}=single(zeros(1,num_T));
            end
       
            tmp_GM_fs{ss}{gg}{ii}=vertcat(tmp_GM_fs{ss}{gg}{ii}, eval(strcat('AMAX_GM_',string(par_opts(pp)),'{',string(num2str(ss)),'}','{',string(num2str(gg)),'}','{',string(num2str(ii)),'}')));
            tmp_GM_fs_ref{ss}{gg}{ii}=vertcat(tmp_GM_fs_ref{ss}{gg}{ii}, eval(strcat('AMAX_GM_ref_',string(par_opts(pp)),'{',string(num2str(ss)),'}','{',string(num2str(gg)),'}','{',string(num2str(ii)),'}')));
        
            if gg==num_GMPEindex && ss==num_syncat
                
                for kk=1:num_T
                
                    if pp==1
                        tmp_GM_fs_mean{ii}{kk}=single(zeros(1,4));
                        tmp_GM_fs_mean_ref{ii}{kk}=single(zeros(1,4));
                    end
                
                tmp_GM_fs_mean{ii}{kk}=vertcat(tmp_GM_fs_mean{ii}{kk}, squeeze(eval(strcat('AMAX_GM_mean_',string(par_opts(pp)),'{',string(num2str(ii)),'}','{',string(num2str(kk)),'}'))));
                tmp_GM_fs_mean_ref{ii}{kk}=vertcat(tmp_GM_fs_mean_ref{ii}{kk}, squeeze(eval(strcat('AMAX_GM_mean_ref_',string(par_opts(pp)),'{',string(num2str(ii)),'}','{',string(num2str(kk)),'}'))));
                end
               
            end %end loop for GM mean
            
        end %end ii loop for site
        
       end %end gg loop for GMPE     
         
  end %end ss loop for syncat
   
   clear(strcat('AMAX_GM_bg_',string(par_opts(pp)))); clear(strcat('AMAX_GM_bg_ref_',string(par_opts(pp))));
   clear(strcat('AMAX_GM_',string(par_opts(pp)))); clear(strcat('AMAX_GM_ref_',string(par_opts(pp))));
   clear(strcat('AMAX_GM_mean_bg_',string(par_opts(pp)))); clear(strcat('AMAX_GM_mean_bg_ref_',string(par_opts(pp))));
   clear(strcat('AMAX_GM_mean_',string(par_opts(pp)))); clear(strcat('AMAX_GM_mean_ref_',string(par_opts(pp))));
end

%% Extract only num_need from each syncat

AMAX_GM_bg=cell(num_GMPEindex,1); AMAX_GM_bg_ref=cell(num_GMPEindex,1); 
AMAX_GM_fs=cell(num_syncat,1); AMAX_GM_fs_ref=cell(num_syncat,1); 
AMAX_GM_mean_fs=cell(num_site,1); AMAX_GM_mean_fs_ref=cell(num_site,1); 

for gg=1:num_GMPEindex
    
     for ii=1:num_site
         
       for kk=1:num_T  
        tmp=sortrows(tmp_GM_bg{gg}{ii}(:,kk),1,'descend');
        tmp_ref=sortrows(tmp_GM_bg_ref{gg}{ii}(:,kk),1,'descend');
   
        AMAX_GM_bg{gg}(:,kk,ii)=flip(tmp(1:num_need,:));
        AMAX_GM_bg_ref{gg}(:,kk,ii)=flip(tmp_ref(1:num_need,:));
        
        if gg==num_GMPEindex
             tmp=sortrows(tmp_GM_mean_bg{ii}{kk},1,'descend'); tmp_ref=sortrows(tmp_GM_mean_bg_ref{ii}{kk},1,'descend');
             AMAX_GM_mean_bg{ii}{kk}=tmp(1:num_need,:);   AMAX_GM_mean_bg{ii}{kk}=flip(AMAX_GM_mean_bg{ii}{kk});
             AMAX_GM_mean_bg_ref{ii}{kk}=tmp_ref(1:num_need,:); AMAX_GM_mean_bg_ref{ii}{kk}=flip(AMAX_GM_mean_bg_ref{ii}{kk});
        end
        
       end
       
     end
  
end

%Sort fs GM
for ss=1:num_syncat
    
    AMAX_GM_fs{ss}=cell(num_GMPEindex,1);  AMAX_GM_fs_ref{ss}=cell(num_GMPEindex,1);
 
    for gg=1:num_GMPEindex
        
    for ii=1:num_site
        
        for kk=1:num_T
            
        tmp = sortrows(tmp_GM_fs{ss}{gg}{ii}(:,kk),'descend'); tmp_ref = sortrows(tmp_GM_fs_ref{ss}{gg}{ii}(:,kk),'descend');
        %select only num_need from combination of all parallelizations, make sure values are ascending
        AMAX_GM_fs{ss}{gg}(:,kk,ii) = flip(tmp(1:num_need));AMAX_GM_fs_ref{ss}{gg}(:,kk,ii) = flip(tmp_ref(1:num_need)); 
        end
     
    end
    end
           
end

%Extract num_need mean GMM
for ii=1:num_site
    for kk=1:num_T
        tmp=sortrows(tmp_GM_fs_mean{ii}{kk},1,'descend');
        tmp_ref=sortrows(tmp_GM_fs_mean_ref{ii}{kk},1,'descend');
        AMAX_GM_mean_fs{ii}{kk}=flip(tmp(1:num_need_fs,:));
        AMAX_GM_mean_fs_ref{ii}{kk}=flip(tmp_ref(1:num_need_fs,:));
    end
end

%% Combine bg and fs ground motions

AMAX_GM_mean_cb = zeros(num_need_fs,num_T,num_site);
AMAX_GM_mean_cb_ref = zeros(num_need_fs,num_T,num_site);

for ss=1:num_syncat
   
    AMAX_GM_cb{ss}=cell(num_GMPEindex,1);  AMAX_GM_cb_ref{ss}=cell(num_GMPEindex,1);
    
    for gg = 1:num_GMPEindex
    
        for ii=1:num_site
      
            for kk=1:num_T
            %combine and sort AMAX for each syncat and GMPE combinations
            tmp = sort([AMAX_GM_fs{ss}{gg}(:,kk,ii); AMAX_GM_bg{gg}(:,kk,ii)]);%Combine and sort bg and fs GM values
            AMAX_GM_cb{ss}{gg}(:,kk,ii) = tmp(num_need+1:2*num_need,:);%Only select num_need highest values from combined record

            tmp = sort([AMAX_GM_fs_ref{ss}{gg}(:,kk,ii); AMAX_GM_bg_ref{gg}(:,kk,ii)]);%Combine and sort bg and fs GM values
            AMAX_GM_cb_ref{ss}{gg}(:,kk,ii) = tmp(num_need+1:2*num_need,:);%Only select num_need highest values from combined record
            
                %combine mean GMM from all realisations
                if gg==num_GMPEindex && ss==num_syncat
                    tmp_bg_dup=repmat(AMAX_GM_mean_bg{ii}{kk}(:,1),obs_duration_fs/obs_duration,1);
                    tmp=sort([AMAX_GM_mean_fs{ii}{kk}(:,1); tmp_bg_dup],'descend');
                    AMAX_GM_mean_cb(:,kk,ii)= flip(tmp(1:num_need_fs,:));
                    
                    tmp_bg_dup_ref=repmat(AMAX_GM_mean_bg_ref{ii}{kk}(:,1),obs_duration_fs/obs_duration,1);
                    tmp=sort([AMAX_GM_mean_fs_ref{ii}{kk}(:,1); tmp_bg_dup_ref],'descend');
                    AMAX_GM_mean_cb_ref(:,kk,ii)= flip(tmp(1:num_need_fs,:));
                end
            
            end

        end
    
    end
end


save('GM_20221021','AMAX_GM_bg','AMAX_GM_bg_ref','AMAX_GM_fs','AMAX_GM_fs_ref',...
   'AMAX_GM_mean_bg','AMAX_GM_mean_bg_ref','AMAX_GM_mean_fs','AMAX_GM_mean_fs_ref',...
    'AMAX_GM_mean_cb','AMAX_GM_mean_cb_ref','AMAX_GM_cb','AMAX_GM_cb_ref');



%% Summary results: Sort data

F(:,1) = ((1:num_need)'+obs_duration-num_need)/(obs_duration+1); 
F_fs(:,1) = ((1:num_need_fs)'+obs_duration_fs-num_need_fs)/(obs_duration_fs+1); 

PSA_fractile_CDF_bg=zeros(num_GMPEindex,length(prob_level),num_T,num_site);
PSA_fractile_CDF_bg_ref=zeros(num_GMPEindex,length(prob_level),num_T,num_site);
PSA_fractile_CDF_fs=zeros(num_syncat,num_GMPEindex,length(prob_level),num_T,num_site);
PSA_fractile_CDF_fs_ref=zeros(num_syncat,num_GMPEindex,length(prob_level),num_T,num_site);
PSA_fractile_CDF_cb=zeros(num_syncat,num_GMPEindex,length(prob_level),num_T,num_site);
PSA_fractile_CDF_cb_ref=zeros(num_syncat,num_GMPEindex,length(prob_level),num_T,num_site);


for ss=1:num_syncat
    
    for gg=1:num_GMPEindex
        
         for pp=1:length(prob_level)
             for ii=1:num_site
             
                if ss==num_syncat
                    PSA_fractile_CDF_bg(gg,pp,:,ii) = AMAX_GM_bg{gg}(num_need-round(obs_duration/prob_level(pp))+1,:,ii);
                    PSA_fractile_CDF_bg_ref(gg,pp,:,ii) = AMAX_GM_bg_ref{gg}(num_need-round(obs_duration/prob_level(pp))+1,:,ii);         
                end
                
            PSA_fractile_CDF_fs(ss,gg,pp,:,ii) = AMAX_GM_fs{ss}{gg}(num_need-round(obs_duration/prob_level(pp))+1,:,ii);
            PSA_fractile_CDF_fs_ref(ss,gg,pp,:,ii)  = AMAX_GM_fs_ref{ss}{gg}(num_need-round(obs_duration/prob_level(pp))+1,:,ii);         
            PSA_fractile_CDF_cb(ss,gg,pp,:,ii) = AMAX_GM_cb{ss}{gg}(num_need-round(obs_duration/prob_level(pp))+1,:,ii);
            PSA_fractile_CDF_cb_ref(ss,gg,pp,:,ii) = AMAX_GM_cb_ref{ss}{gg}(num_need-round(obs_duration/prob_level(pp))+1,:,ii);
                
             
             end
         end
    end
end

%% Figure options

label_opt=vertcat(["(a) ","(b) ","(c) "],["(d) ","(e) ","(f) "],["(g) ","(h) ","(i) "]);
prob_level_label=["10% PoE in 50 years","2% PoE in 50 years"];%update if necessary
prob_level_label2=vertcat(["10% PoE " " in 50 years"],["2% PoE " " in 50 years"]);

%Select vs30 option: 1 =300 m/s or USGS value, 2 = ref value
vs30_option=2;

if vs30_option ==1
     AMAX_GM_bg_tmp =  AMAX_GM_bg; AMAX_GM_fs_tmp =  AMAX_GM_fs;
     AMAX_GM_cb_tmp =  AMAX_GM_cb;
     PSA_fractile_bg_tmp=PSA_fractile_CDF_bg; PSA_fractile_fs_tmp=PSA_fractile_CDF_fs;
     PSA_fractile_cb_tmp=PSA_fractile_CDF_cb;
     AMAX_GM_mean_fs_tmp = AMAX_GM_mean_fs; AMAX_GM_mean_bg_tmp = AMAX_GM_mean_bg;
     AMAX_GM_mean_cb_tmp = AMAX_GM_mean_cb;
    
elseif vs30_option==2
     AMAX_GM_bg_tmp =  AMAX_GM_bg_ref; AMAX_GM_fs_tmp =  AMAX_GM_fs_ref;
     AMAX_GM_cb_tmp =  AMAX_GM_cb_ref;
     PSA_fractile_bg_tmp=PSA_fractile_CDF_bg_ref; PSA_fractile_fs_tmp=PSA_fractile_CDF_fs_ref;
     PSA_fractile_cb_tmp=PSA_fractile_CDF_cb_ref;
     AMAX_GM_mean_fs_tmp = AMAX_GM_mean_fs_ref; AMAX_GM_mean_bg_tmp = AMAX_GM_mean_bg_ref;
     AMAX_GM_mean_cb_tmp = AMAX_GM_mean_cb_ref;
 
end


%% Seismic hazard curve for each spectal acceleration
% Also obtain ground acceleration for prob_level(pp) for each SA

figure (100);
tiledlayout(3,num_site,'tilespacing','compact')

gmpe_marker=(["o","x","v","^"]);
m_shift=1:2000:10000;%REPLACE AND COMMENT OUT WHEN USING BC
%m_shift=400:2500;
syncat_col=vertcat([0.8500, 0.3250, 0.0980],[0, 0.5, 0],[0.75 0.75 0],...
    [0.4940, 0.1840, 0.5560],[0, 0, 1]);
%generate 5 equally space intervals between markers

poe_text=["10% PoE" "2% PoE"]; tmp_shift=[0.001 0.0002];

%Select Spectal Acceleration for plot figure
select_T=8; %1 = PGA, 8 = 3 s

if select_T==1
     x_lim=[0 0.6]; tmp_shift1=0.2;
else select_T==8
     x_lim=[0 0.1];tmp_shift1=0.05;
end

for ii=1:num_site
   
    nexttile
    tmp = zeros(num_need,num_syncat*num_GMPEindex);
    count=1; 
    
    for gg=1:num_GMPEindex
        
        for ss=1:num_syncat
             
                marker_space = round(length(F)+1-logspace(1,log10(length(F)-m_shift(gg)),10));
              
                hazard_curve(ss,gg)=loglog(AMAX_GM_cb_tmp{ss}{gg}(:,select_T,ii),1-F,'Color',syncat_col(ss,:),'Marker',gmpe_marker(gg),...
                    'MarkerIndices',marker_space,'MarkerEdgeColor',syncat_col(ss,:),'LineWidth',1.2,'MarkerSize',4); hold on;
                tmp(:,count)=AMAX_GM_cb_tmp{ss}{gg}(:,select_T,ii);
                count=count+1;
               
         end
     end
    
    %plot mean curve
    m2=loglog(mean(tmp'),1-F,'k-','LineWidth',2); hold on;
    
    %plot horizontal line for prob_level
    for pp=1:length(prob_level)
        loglog([0.01 10],[1/prob_level(pp) 1/prob_level(pp)],'k--','LineWidth',1.5); hold on;
        text(0.75,(1/prob_level(pp)+tmp_shift(pp)),[poe_text(pp) 'in 50 years'],'FontSize',11);
    end
    
    grid on;
    xlabel('Acceleration (g)'); ylabel(['Annual probability' newline ' of exceedance']); title(strcat(label_opt(1,ii),site_name(ii)),'Fontweight','normal'); axis([0.01 10 0.0001 0.01]); axis square;
    set(gca,'fontsize',11); 
    %clear tmp
end


[legh,objh] =legend([line([1,1],[1,1],'Color',syncat_col(1,:),'LineWidth',1.2),...
    hazard_curve(2:num_syncat,1)',hazard_curve(1,1:num_GMPEindex),m2],...
    {'Direct MSSM','length-limited, \it{W}, Char','length-limited, \it{W}, G-R','layer-limited \it{W}, char','layer-limited, \it{W} G-R',...
    'Boore2014','Akkar2014','Chiou2014','Atkinson2013','mean hazard curve'},...
    'FontSize',9,'Location','eastoutside');
%Shrink markers so not visible on legend
lineh = findobj(objh(11:21),'type','line'); set(lineh,'MarkerSize',0.1);
%Change some line colors to grey
linei = findobj(objh(21:28),'type','line'); set(linei,'Color',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);

for pp=1:length(prob_level)

for ii=1:num_site
    
    PSA_fractile=zeros(num_syncat,num_GMPEindex);
    
    nexttile
    
    for ss=1:num_syncat
        for gg=1:num_GMPEindex
            PSA_fractile(ss,gg)=PSA_fractile_cb_tmp(ss,gg,pp,select_T,ii);
            PSA_fractile_col=reshape(PSA_fractile,num_syncat*num_GMPEindex,1);
        end
    end
    
    
    x_values=[0.001:0.001:max(PSA_fractile_col)+tmp_shift1];
    pd=fitdist(PSA_fractile_col,'beta');
    y_pdf=pdf(pd,x_values);
    
    plot(x_values,y_pdf,'k-','LineWidth',1.3); 
    xlabel('Acceleration (g)'); 
    xlim([0 max(PSA_fractile_col)+tmp_shift1])
    ylabel('PDF')
    
    
    %Plot discrete values
    yyaxis right
  
    for gg=1:num_GMPEindex
        for ss=1:num_syncat
        line([PSA_fractile(ss,gg),PSA_fractile(ss,gg)],[0,1/(num_syncat*num_GMPEindex)],'Color',syncat_col(ss,:),'LineStyle','-','LineWidth',1.2);hold on
        plot(PSA_fractile(ss,gg),ones(height(PSA_fractile(ss,gg)),1)/(num_syncat*num_GMPEindex),'MarkerEdgeColor',syncat_col(ss,:),...
            'MarkerFaceColor','w' ,'Marker',gmpe_marker(gg),'LineWidth',1.2);hold on  
        end
    end
    
    title(strcat(label_opt(pp+1,ii),site_name(ii)),'Fontweight','normal'); xlim(x_lim); ylim([0 0.2]); ylabel('Weighting');
    set(gca,'fontsize',11); axis square; set(gca,'YColor',[0.35 0.35 0.35]);
    
    clear x_values y_pdf 
    
end
end

set(gcf,'Position',[341 83 972 714]);

%% Plot ecdf of PoE values at 2% and 10% PoE in 50 years

NumSimu1=1000;

figure(101);
tiledlayout(2,num_site,'tilespacing','compact')

for pp=1:length(prob_level)
    
    for ii=1:num_site  

        nexttile
        
        count=0;

        PSA_fractile=zeros(num_syncat,num_GMPEindex);

        for ss=1:num_syncat
            for gg=1:num_GMPEindex
                PSA_fractile(ss,gg)=PSA_fractile_cb_tmp(ss,gg,pp,select_T,ii);
                PSA_fractile_col=reshape(PSA_fractile,num_syncat*num_GMPEindex,1);
            end
        end

        x_values=[0.001:0.001:max(PSA_fractile_col)+tmp_shift1];
        pd=fitdist(PSA_fractile_col,'beta');


        h1=cdfplot(PSA_fractile_col); hold on
        h1.LineWidth=1;
        h2_val=betacdf(sort(PSA_fractile_col),pd.a,pd.b);
        h2=plot(sort(PSA_fractile_col),h2_val,'LineWidth',1);

        for jj=1:NumSimu1

        [test_val,kstest_stat]=kstest2(betarnd(pd.a,pd.b,[20,1]),PSA_fractile_col);

            if test_val==0
                count=count+1;
            end
        end
       
        title(strcat(label_opt(pp,ii),site_name(ii),'; count: ',num2str(count)),'Fontweight','normal'); 
        set(gca,'fontsize',11); axis square; set(gca,'YColor',[0.35 0.35 0.35]); xlabel('Acceleration (g)');
        legend('ECDF','Beta CDF','Location','southeast');
    end
end

set(gcf,'Position',[440 156 866 641]);

%% Plot uniform hazard spectra
T_taken = [1 4 6 8]; 
T(T_taken);

T2 = T; T2(find(T==0)) = 0.05;

% Uniform hazard spectra
figure(1)
marker1=vertcat(["ko-" "ro-" "bo-"],["kx-" "rx-" "bx-"]);
marker2=(["o-","x-"]);
color2=vertcat([0.8500, 0.3250, 0.0980],[0 0.5 0.0],[0.4940 0.1840 0.5560]);


for pp = 1:length(prob_level)

  %Plot by source for different sites and return periods
    h1=subplot(2,4,1);
        
    for ii=1:num_site
        
        mean_fractile=squeeze(mean(PSA_fractile_fs_tmp(:,:,pp,:,ii),[1 2]));
        loglog(T2,mean_fractile,marker1(pp,ii),'Linewidth',1); hold on; grid on;
        if ii == 1; xlabel('Period (s)'); ylabel('Acceleration (g)'); axis([0.01 10 0.01 10]); axis square; title('(a) MSSM','Fontweight','normal');  end
        set(gca,'fontsize',12); hold on
    end
    
    subplot(2,4,2);
    for ii=1:num_site
        
        mean_fractile=squeeze(mean(PSA_fractile_bg_tmp(:,pp,:,ii),[1 2]));
        loglog(T2,mean_fractile,marker1(pp,ii),'Linewidth',1); hold on; grid on;
        if ii == 1; xlabel('Period (s)'); ylabel('Acceleration (g)'); axis([0.01 10 0.01 10]); axis square; title('(b) Areal Sources','Fontweight','normal'); end
        set(gca,'fontsize',12); hold on
    end
    
    subplot(2,4,3);
    for ii=1:num_site
        
        mean_fractile=squeeze(mean(PSA_fractile_cb_tmp(:,:,pp,:,ii),[1 2]));
        loglog(T2,mean_fractile,marker1(pp,ii),'Linewidth',1); hold on; grid on;
        if ii == 1; xlabel('Period (s)'); ylabel('Acceleration (g)'); axis([0.01 10 0.01 10]); axis square; title('(c) Combined','Fontweight','normal'); end
        set(gca,'fontsize',12);hold on
    end
    
    if pp==length(prob_level) && ii==num_site
        
     subplot(2,4,4);   
     ax1=gca; set(ax1,'visible','off');
     L1=legend(h1,{[site_name(1) + newline + prob_level_label(1)],[site_name(2) + newline + prob_level_label(1)],[site_name(3) + newline + prob_level_label(1)],...
        prob_level_label(2),prob_level_label(2),prob_level_label(2)},'Position',[0.7536 0.596 0.1848 0.15],'NumColumns',2,'FontSize',11);   
 
    end
    
    %Plot by sites for different sources and return periods
    
     for ii=1:num_site
         if pp==1
            h3=subplot(2,4,ii+4);
         else
            subplot(2,4,ii+4)
         end        
         
          mean_fractile_fs=squeeze(mean(PSA_fractile_fs_tmp(:,:,pp,:,ii),[1 2]));
          mean_fractile_bg=squeeze(mean(PSA_fractile_bg_tmp(:,pp,:,ii),[1 2]));
          mean_fractile_cb=squeeze(mean(PSA_fractile_cb_tmp(:,:,pp,:,ii),[1 2]));
         
        loglog(T2,mean_fractile_fs,marker2(pp),'Linewidth',1,'Color',color2(1,:)); hold on
        loglog(T2,mean_fractile_bg,marker2(pp),'Linewidth',1,'Color',color2(2,:)); hold on
        loglog(T2,mean_fractile_cb,marker2(pp),'Linewidth',1,'Color',color2(3,:)); hold on; grid on;
        if pp == 1; xlabel('Period (s)'); ylabel('Acceleration (g)'); axis([0.01 10 0.01 10]); axis square; title(strcat(label_opt(2,ii),site_name(ii)),'Fontweight','normal'); end
        set(gca,'fontsize',12);hold on
     end
     
     if pp==length(prob_level) && ii==num_site
     subplot(2,4,8);   
     ax1=gca; set(ax1,'visible','off');
     L2=legend(h3,{["MSSM Sources: " + newline + prob_level_label(1)],["Areal Sources: " + newline + prob_level_label(1)],["Combined: " + newline + prob_level_label(1)],...
        prob_level_label(2),prob_level_label(2),prob_level_label(2)},'Position',[0.7536 0.1214 0.1848 0.15],'NumColumns',2,'FontSize',11); hold on
     end

end

set(gcf,'Position',[166 243 1128 554]);

%% Seismic disaggregation (exceeding samples)

select_T=1;%SA for testing

fignum_disagg = 1000; tmp_idx=0;
figure(fignum_disagg);


for pp = 1:length(prob_level)
    
    for ii = 1:num_site
    
    tmp_idx=tmp_idx+1;
    
    subplot(length(prob_level),num_site,tmp_idx)
    
        title_opt=[label_opt(pp,ii),site_name(ii),prob_level_label2(pp,:),num2str(mean(mean(PSA_fractile_cb_tmp(:,:,pp,select_T,ii))),2)];
       
        % Combine sources
        AMAX_GM_bg_dup = sortrows(repmat(AMAX_GM_mean_bg_tmp{ii}{select_T},num_syncat,1),1,'descend');%repeat so are equivalen length
        AMAX_GMtmp_cb = [AMAX_GM_mean_fs_tmp{ii}{select_T}; AMAX_GM_bg_dup];
        AMAX_GMtmp_cb = sortrows(AMAX_GMtmp_cb,1,'ascend');%values need to be ascending so right prob_level indexed
        AMAX_GMtmp_cb = AMAX_GMtmp_cb(num_need_fs+1:2*num_need_fs,:);
        disagg_MSSM_JW(num_need_fs,F_fs,200,AMAX_GMtmp_cb,prob_level(pp),T(select_T),fignum_disagg,0.1,title_opt);
        
        hold on
    end
    hold on
end

set(gcf,'Position',[77 117 1323 680]);

%% Mean PGA values for different MSSM event catalogs

tmp=zeros(num_GMPEindex,num_site); source_pga=zeros(num_syncat,length(prob_level),num_site);
tmp_ref=zeros(num_GMPEindex,num_site); source_pga_ref=zeros(num_syncat,length(prob_level),num_site);

for ss=1:num_syncat

    for pp=1:length(prob_level)   
    
        for gg=1:num_GMPEindex
       
            tmp(gg,:) = squeeze(AMAX_GM_cb{ss}{gg}(num_need-round(obs_duration/prob_level(pp))+1,select_T,:));
            tmp_ref(gg,:) = squeeze(AMAX_GM_cb_ref{ss}{gg}(num_need-round(obs_duration/prob_level(pp))+1,select_T,:));
        end
        
        for ii=1:num_site
            
            source_pga(ss,pp,ii)=mean(tmp(:,ii));
            source_pga_ref(ss,pp,ii)=mean(tmp_ref(:,ii));
        end
        
    end 
end

function [points,points_matrix] = points_aggregator_capture_strategy(radiodata,sv_hist_1,sv_hist_2,sv_hist_3,az_del,el_del);

points_matrix =[]; % [Source ID, observation number, h, hmid, hmax, P, delta, points]

source_list=radiodata(:,1);
az_sources=radiodata(:,2);
el_sources=radiodata(:,3);

for i=1:(length(sv_hist_1)-1)
    
    % check for normal 1 alone
    az_1 = az_del(i);
    el_1 = el_del(i);
    
    az_2 = az_del(i+1);
    el_2 = el_del(i+1);
    
    % Reduce source search space
    min_az=min([az_1 az_2]);
    max_az=max([az_1 az_2]);
    id_az = ((az_sources>min_az-0.2)&(az_sources<max_az));
    
    min_el=min([el_1 el_2]);
    max_el=max([el_1 el_2]);
    id_el = ((el_sources>min_el-0.2)&(el_sources<max_el));
    
    id_common =id_el & id_az;
    
    az_temp_sources =az_sources(id_common);
    el_temp_sources =el_sources(id_common);
    
    % Solve equations to get the line eqn coeffs m,c of y= mx + c
    A=[az_1 1;az_2 1];
    B=[el_1;el_2];
    X=A\B;
    
    % In ax+by+c=0 form
    a=X(1);
    b=-1;
    c=X(2);
    
    % Find distance of all sources with respect to the lines
    if length(az_temp_sources)>0 && length(el_temp_sources)>0
        distance_matrix = abs(a*az_temp_sources+b*el_temp_sources+c)./norm([a b]);
    else
        distance_matrix = 0;
    end
    
    if sum(distance_matrix<0.1)~=0 && (abs(az_1-az_2)<90) && (abs(el_1-el_2)<45) && length(az_temp_sources)>0 && length(el_temp_sources)>0
        [~,idx_query]=min(distance_matrix);
        capture_az=az_temp_sources(idx_query);
        capture_el=el_temp_sources(idx_query);
                
        idx = (capture_az==az_sources) & (capture_el==el_sources);
        idx = source_list(idx);
        
        r1=sv_hist_1(i,1:3);
        r2=sv_hist_2(i,1:3);
        r3=sv_hist_3(i,1:3);
        
        h1= norm(cross(r2-r1,r3-r1))/norm(r3-r1);
        h2= norm(cross(r3-r2,r1-r2))/norm(r3-r2);
        h3= norm(cross(r1-r3,r2-r3))/norm(r1-r3);
        
        h_sort = sort( [h1 h2 h3] );
        
        h= h_sort(1);
        hmid= h_sort(2);
        hmax= h_sort(3);
        
        source_id = radiodata(idx,1);
        
        dec = radiodata(idx,3);
        
        % check if the least altitude is greater than 10,000 km
        if h > 1e4
            
            if length(points_matrix)==0
                ob_count = 0;
            else
                ob_count = sum(points_matrix(:,1)==source_id);
            end
            
            % Find if this the first observation
            if ob_count ==0
                
                P=1;
                
            elseif ob_count==1
                
                if  hmax/h >= 3
                    P=3;
                else
                    P=1;
                end
                
            elseif ob_count==2
                
                temp = find(points_matrix(:,1)==source_id);
                
                P_second_observation = points_matrix(temp(2),6);
                
                if hmax/hmid >=3 && hmid/h >=3
                    
                    P=6;
                    
                elseif hmax/h >=3 && P_second_observation == 1
                    
                    P=3;
                    
                else
                    
                    P=1;
                    
                end
                
            elseif ob_count==3
                
                P=0;
                
            end
            
        end
        
        if length(points_matrix)~=0
            points_prev= points_matrix(end,end);
        else
            points_prev=0;
        end
        
        
        points_curr= points_prev + P * h * (0.2 + cosd(dec)^2);
        points_matrix = [points_matrix; source_id  ob_count h hmid hmax P dec points_curr];
        
    end
end


if length(points_matrix) >0
    points = points_matrix(end,end);
else
    points = 0;
end

end


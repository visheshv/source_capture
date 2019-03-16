function [points,points_matrix] = points_aggregator(radiodata,sv_hist_1,sv_hist_2,sv_hist_3,az_del,el_del);

points_matrix =[]; % [Source ID, observation number, h, hmid, hmax, P, delta, points]

for i=1:length(sv_hist_1)
    
    for j=1:2
    % check for normal 1/2
    
    if j==1
        az_check = az_del(i);
        el_check = el_del(i);
    elseif j==2
        az_check = az_del(i+length(sv_hist_1));
        el_check = el_del(i+length(sv_hist_1));
    end
    
    angle = [az_check el_check];
    
    idx = knnsearch(radiodata(:,2:3),angle);
    
    error = norm(radiodata(idx,2:3) - angle);
   
    if error<0.1
        
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
    
    
end
        
        if length(points_matrix) >0
            points = points_matrix(end,end);
        else
            points = 0;
        end

end
    
   
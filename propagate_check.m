function [state_store,capture_states] = propagate_check(y1,y2,y3,tdur)

global mu_earth R_earth radiodata n_segments y_initial u_current

dtcsv = 60 * 24 * 3; % minutes in 3 days
ti = 0.0d0;
tf = 0.0d0;
tetol=1.0e-12;
tsim = 86400.0d0 * tdur;
neq=18;
state_store=zeros(floor(tsim/(dtcsv*60)),18);
counter=0;
yi = [y1;y2;y3];

u_current=zeros(9,1);

while(1)
    
    % step size guess (seconds)

    counter=counter+1;
    
    h = 10.0;
    
    ti = tf;
    
    tf = ti + 60.0d0 * dtcsv;
    
    % integrate from ti to tf
    
    yfinal = rkf78('sel_eqm_modified', neq, ti, tf, h, tetol, yi);
    
    % compute current state vector
    r1 = yfinal(1:3);
    r2 = yfinal(7:9);
    r3 = yfinal(13:15);
    
    [az_1,el_1]= calculate_celestial_track(r1',r2',r3');


    temp=(radiodata(:,2)-az_1).^2+(radiodata(:,3)-el_1).^2 - 5^2;        %Check to see if lying within 5 deg
    
    if sum(temp<0) ==0
        % No action
    else
        %% Capture the source
        % Find the source to be captured
        [~,idx]=min(temp);
        source_target=radiodata(idx,2:3); % Targeted source position
        
        y_initial = yi;
        
        % Call NLP
        xg = zeros(9*n_segments+1,1);
        xlb(1) = dtcsv * 60 * 0.1;                                  % Minimum transfer time 0.3 days
        xlb(2:(9*n_segments+1),1) = -0.1 * ones(9*n_segments,1);
        xub(1) = dtcsv * 60 * 1.1;                                  % Minimum transfer time 3.3 days
        xub(2:(9*n_segments+1),1) =  0.1 * ones(9*n_segments,1);
        
        % f = [transfer time (objective function); az_n_end; el_n_end;  magnitude of thrust S1;S2;S3]  
        
        flow(1) = 0;
        flow(2) = source_target(1)-0.07;
        flow(3) = source_target(2)-0.07;
        flow(4) = 0;
        flow(5) = 0;
        flow(6) = 0;
        
        fupp(1) = +Inf;
        fupp(2) = source_target(1)+0.07;
        fupp(3) = source_target(2)+0.07;
        fupp(4) = 0.1^2;
        fupp(5) = 0.1^2;
        fupp(6) = 0.1^2;     
        
        flow = flow';
        fupp = fupp';
        
        xmul = zeros(9 * n_segments + 1, 1); % No info about dual
        
        xstate = zeros(9 * n_segments + 1, 1); % No desire to provide special information
        
        fmul = zeros(6, 1); % No information about Lagrange multipliers
        
        fstate = zeros(6, 1); % No special information
        
        snspec('capture_specs.txt');

        snscreen('on');
        
        [x, f, inform, xmul, fmul] = snopt(xg, xlb, xub, xmul, xstate,flow, fupp, fmul, fstate, 'capture_track');
        
        % Assign values
        
        % compute duration of each time interval (non-dimensional)
        deltat = x(1) / n_segments;
        
        % specify number of differential equations
        neq = 18;
        
        % truncation error tolerance
        tetol = 1.0e-10;
        
        % initialize initial time
        ti = -deltat;
        
        % set final time to current thrust duration
        tof = x(1);
        
        % step size guess (seconds)
        h = 30.0;
        
        capture_states= zeros(n_segments , 18);
        
        % Propagate using the achieved control history
        for i = 1:1:n_segments
            
            % current thrust values for Tx,Ty,Tz for S1,S2,S3
            u_current = x( 9 * i - 7 : 9 * i + 1 );
            
            % increment initial and final times
            ti = ti + deltat;
            tf = ti + deltat;
            
            % integrate from current ti to tf
            yfinal = rkf78('sel_eqm_modified', neq, ti, tf, h, tetol, yi);
            
            yi = yfinal;
           
            capture_states(i,:) = yi';
            
            % check for end of simulation
            if (tf >= tof)
            
                break;
                
            end
            
        end
        
        r1 = yfinal(1:3);
        r2 = yfinal(7:9);
        r3 = yfinal(13:15);
        normal_1= cross(r2-r1,r3-r1);
        [az_1, el_1, ~]= cart2sph(normal_1(1),normal_1(2),normal_1(3));
        
        % Rad to deg conversion
        az_1(az_1<0) = 2*pi +az_1(az_1<0);
        az_1 = (180/pi) * az_1;  % Achieved source position
        el_1 = (180/pi) * el_1;
            
        
        state_store(counter,:)=yfinal';

        break % after the first source hit
       
    end
    
    state_store(counter,:)=yfinal';
    
    yi = yfinal;
    
    % check for end of simulation
    
    if (tf >= tsim)
        
        break;
        
    end
    
end
state_store( ~any(state_store,2), : ) = [];  %rows



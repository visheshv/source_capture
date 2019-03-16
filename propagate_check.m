function [state_store,capture_states] = propagate_check(y1,y2,y3,tdur)
% Function to propagate state space equations for satellite
% position-velocity vectors and plan path course corrections if a source is
% within a 5 deg band of the final point of propagation. The path plan is
% done by using the 0.1 N of thrust along all axes of the spacecraft
% S1/S2/S3 and achieves the source position within 0.1 deg of accuracy.

% Inputs: y1,y2,y3: Initial state [r,v] of S1/S2/S3, i.e. y1,y2,y3, (km,km/s)
%         tdur: Duration of propagation (days)
% Output: state_store: State history, 
%         capture_states: finer course change history if a new path is
%         planned  [km, km/s]

%%  Initialize states, time step, other propagation params
global mu_earth R_earth radiodata n_segments y_initial u_current

dtcsv = 60 * 24 * 3; % duration of propagation (minutes in 3 days)
ti = 0.0d0;          % Initial time 
tf = 0.0d0;          % Final time
tetol=1.0e-12;       % Tolerance
tsim = 86400.0d0 * tdur;    
neq=18;              % Number of equations, dr/dt [3x1] and dv/dt [3x1] for all spacecraft
state_store=zeros(floor(tsim/(dtcsv*60)),18); % Store states through propagation
counter=0;
yi = [y1;y2;y3];

u_current=zeros(9,1);

while(1)
    
    % step size guess (seconds)

    counter=counter+1;
    
    h = 10.0; % [s]
    
    ti = tf;
    
    tf = ti + 60.0d0 * dtcsv;
    
    % integrate from ti to tf
    
    yfinal = rkf78('sel_eqm_modified', neq, ti, tf, h, tetol, yi); % Propagate all states 
    
    % compute final state vector
    r1 = yfinal(1:3);
    r2 = yfinal(7:9);
    r3 = yfinal(13:15);
    
    [az_1,el_1]= calculate_celestial_track(r1',r2',r3');

    % Check if a radiosource is close to the final position of the triangle
    % normal
    temp=(radiodata(:,2)-az_1).^2+(radiodata(:,3)-el_1).^2 - 5^2;        %Check to see if lying within 5 deg circle of the final position of triangle normal
    
    if sum(temp<0) ==0
        % No action if there no source nearby
    else
        %% Capture the source if a source is found
        % Find the closest source to be captured
        [~,idx]=min(temp);
        source_target=radiodata(idx,2:3); % Targeted source position
        
        y_initial = yi;
        
        % Call NLP to plan a path to the source
        xg = zeros(9*n_segments+1,1);                               % Guess control history, [Transfer time (s), {Tx,Ty,Tz |S1/S2/S3 for 'n' segments} (N)]
        xlb(1) = dtcsv * 60 * 0.1;                                  % Minimum transfer time 0.3 days
        xlb(2:(9*n_segments+1),1) = -0.1 * ones(9*n_segments,1);    % Thrust bounds minimum
        xub(1) = dtcsv * 60 * 1.1;                                  % Minimum transfer time 3.3 days
        xub(2:(9*n_segments+1),1) =  0.1 * ones(9*n_segments,1);    % Thrust bounds maximum
        
        % f vector consists of objective function limits in the first entry
        % and linear/non linear constraints in the remaining entries
        % f = [transfer time (objective function); az_n_end; el_n_end;  magnitude of thrust S1;S2;S3]  
        
        flow(1) = 0;
        flow(2) = source_target(1)-0.07;    % Requirement of the capture within 0.1 deg, kept as 0.07, since 0.07 sized square has diagonal< 0.1 
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
        
        xmul = zeros(9 * n_segments + 1, 1);   % No info about dual variables available 
        
        xstate = zeros(9 * n_segments + 1, 1); % No desire to provide special information
        
        fmul = zeros(6, 1);                    % No information about Lagrange multipliers
        
        fstate = zeros(6, 1);                  % No special information
        
        snspec('capture_specs.txt');           % Empty file, no specs required

        snscreen('on');
        
        % Call SNOPT solver to find the optimal control history and the
        % value of the objective function
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
        
        %% Recreate the path 
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



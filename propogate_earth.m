function state_store = propogate_earth(yi,tdur)

global mu_earth R_earth

dtcsv = 20; % minutes
ti = 0.0d0;
tf = 0.0d0;
tetol=1.0e-12;
tsim = 86400.0d0 * tdur;
neq=6;
state_store=zeros(floor(tsim/(dtcsv*60)),6);
counter=0;


while(1)
    
    % step size guess (seconds)

    counter=counter+1;
    
    h = 10.0;
    
    ti = tf;
    
    tf = ti + 60.0d0 * dtcsv;
    
    % integrate from ti to tf
    
    yfinal = rkf78('sel_eqm', neq, ti, tf, h, tetol, yi);
    
    % compute current state vector
    
    for i = 1:1:3
        
        rf(i) = yfinal(i);
        
        vf(i) = yfinal(i + 3);
        
    end
    
    % compute current orbital elements
    
    oevf = eci2orb1(mu_earth, rf, vf);
    
    % write mission elapsed time and orbital elements to data file
    
    rp = oevf(1) * (1.0d0 - oevf(2));
    
    ra = oevf(1) * (1.0d0 + oevf(2));
    
    hp = rp - R_earth;
    
    ha = ra - R_earth;
    
    %state_store=[state_store; tf / 86400.0, oevf(1), oevf(2), rtd * oevf(3), rtd * oevf(4), rtd * oevf(5), rtd * oevf(6), hp, ha ];
    state_store(counter,:)=yi';
    
    yi = yfinal;
    
    % check for end of simulation
    
    if (tf >= tsim)
        
        break;
        
    end
    
end
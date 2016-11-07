function [window,wnrange,gasvec,invgas,ninvgas,sol_shift_wn,solar_line_file,mindep] = window_details(gas,window_number,varargin)
%  [window,wnrange,gasvec,invgas,sol_shift_wn,solar_line_file] = window_details(gas,window_number)

mindep = 0.125;
solar_line_file = 'input_data/solar_lines_1500-1700.dat';

switch lower(gas)
    
  case 'o2'

    window = [7765 8005];
    wnrange = [7870 7889.4];
    gasvec = {'o2','h2o','co2'};
    sol_shift_wn = 7884.82; 
    solar_line_file = 'input_data/solar_lines_1200-1300.dat';

    invgas = {'o2','h2o'};
    ninvgas = length(invgas);
    mindep = 0.05;

  case 'co2'
    
    switch window_number
        
      case 1
        
        window = [6180 6260];
        wnrange = window;
        gasvec = {'co2','h2o','hdo','ch4'};
        sol_shift_wn = 6208.55;

    end
    
  case 'ch4'
    
    switch window_number
        
      case 1
        
        window = [5996 6008];        
        wnrange = window;
        gasvec = {'ch4','h2o','co2'};
        sol_shift_wn = 6005.84; %6000.262;
        
        if (nargin==3)
            
            switch varargin{1};
                
              case 1
                
              case 2
                
            end
            
        end
        
      case 2        
        
        window = [5883 5990];
        wnrange = window;
        gasvec = {'ch4','h2o','co2','n2o'}; % no hdo lines
        sol_shift_wn = 5968.3;
        
        % individual lines in this window
        if (nargin==3)
            
            switch varargin{1};
                
              case 1 % ch4, h2o
                
                gasvec = {'ch4','h2o'};
                wnrange = [5971.2 5972.9];
                sol_shift_wn = 5972.33;
                
              case 2 % all gases
                
                wnrange = [5959.4 5962.5];
                sol_shift_wn = 5939.03;
                
              case 3 % 
                
                gasvec = {'ch4','h2o'};
                wnrange = [5937.6 5938.6];
                % sol_shift_wn = 5937.36;
                
              case 4 % ch4
                
                gasvec = {'ch4'};
                wnrange = [5949 5950];
                mindep =  0.15;

              case 5

                gasvec = {'ch4','co2'};
                wnrange = [5982.5 5984];
                sol_shift_wn = 5983.535;
                mindep =  0.18;
                
              case 6
                
                gasvec = {'ch4'};
                wnrange = [5926 5927];
                mindep =  0.15;
                
            end

            invgas = gasvec;
            ninvgas = length(invgas);
            
        end
        
      case 3 % suitable ch4 window for profile retrieval
        
        window = [5996 6008];
        wnrange = [6003 6005.5];
        gasvec = {'ch4','h2o'};
        sol_shift_wn = 6005.85;
      
        % retrieve all gases
        invgas = gasvec;
        ninvgas = length(invgas);        
  
    end

end


function [window,wnrange,gasvec,invgas,ninvgas,sol_shift_wn,solar_line_file] = window_details(gas,window_number,varargin)
%  [window,wnrange,gasvec,invgas,sol_shift_wn,solar_line_file] = window_details(gas,window_number)

solar_line_file = 'input_data/solar_lines_1500-1700.dat';

switch lower(gas)
    
  case 'co'
    
    switch window_number
        
      case 1
        
        window = [4210 4240];
        wnrange = [4210 4240];
        gasvec = {'co'};
        sol_shift_wn = 4232.83;
        solar_line_file = 'input_data/solar_lines_2100-2500.dat';
        
        invgas = {'co'};
        ninvgas = length(invgas);
        
      case 2
        
        window = [4280 4300];
        wnrange = [4280 4300];
        gasvec = {'co','ch4','h2o'};
        sol_shift_wn = 4291.292;
        solar_line_file = 'input_data/solar_lines_2100-2500.dat';
        
        invgas = {'co','ch4','h2o'};
        ninvgas = length(invgas);
        
      case 3
        
        window = [2055 2060];
        wnrange = [2055 2060];
        gasvec = {'co'};
        sol_shift_wn = 4291.292;
        solar_line_file = 'input_data/solar_lines_2100-2500.dat';
        
        invgas = {'co'};
        ninvgas = length(invgas);
    end
    
    
  case 'o2'
    
    window = [7870 7890];
    %wnrange = [7870 7876];
    %wnrange = [7872 7875.5];
    wnrange = [7870 7878];
    gasvec = {'o2','h2o'};
    %sol_shift_wn = 7884.81; 
    sol_shift_wn = 7882.745;
    solar_line_file = 'input_data/solar_lines_1200-1300.dat';

    invgas = {'o2','h2o'};
    ninvgas = length(invgas);

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

              case 5

                gasvec = {'ch4','co2'};
                wnrange = [5982.5 5984];
                sol_shift_wn = 5983.535;
                
              case 6
                
                gasvec = {'ch4'};
                wnrange = [5926 5927];
                
            end

            invgas = gasvec;
            ninvgas = length(invgas);
            
        end
        
      case 3 % suitable ch4 window for profile retrieval
        
        window = [5996 6008];
        wnrange = [6003 6005.15];
        gasvec = {'ch4','h2o'};
        sol_shift_wn = 6005.85;
      
        % retrieve all gases
        invgas = gasvec;
        ninvgas = length(invgas);        
  
    end

end


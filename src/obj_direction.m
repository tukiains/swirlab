function[dir]=obj_direction(sza,saa)

dir(1) = cos(saa)*sin(sza); % x
dir(2) = sin(saa)*sin(sza); % y
dir(3) = cos(sza);          % z
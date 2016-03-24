function [a] = ggg_ils(apo,ns,freq,grid_step,L,fovd)
% [a] = ggg_ils(apo,ns,freq,grid_step,L,fovd)

% FTIR instrument function (windowed) from GGG
%
% apo = level of apozidataion (1-4)
% ns = number of points to be calculated
% freq = center frequency of the window [cm-1]
% grid = desired grid spacing [cm-1]
%
% GGG:
% In normal use, PROFZL will be called twice: Once for the synthetic
% spectrum with the actual values of RESNOG and RECTOG, and once for the measured
% spectrum with RECTOG=0. Convolving the measured spectrum with its own
% (infinite) SINC function would be a do-nothing operation. However,
% convolving the measured spectrum with a finite and weakly apodized version of its
% own SINC function will improve the agreement with the synthetic spectrum.

    resnog = 0.5/L/grid_step;                                                                                  

    rectog = freq*fovd^2/8/grid_step;

    c(:,1) = [1.0, 0.5480, 0.2600, 0.0900];
    c(:,2) = [0.0, -.0833, -.154838, 0.00];
    c(:,3) = [0.0, 0.5353, .894838, .5875];
    c(:,4) = [0.0, 0.0000, 0.0000, 0.3225];
    
    hwid = 0.5*(ns-1);
    
    np = 2 + round(4*rectog/resnog);
    
    del = rectog/sqrt(np*np-1.0);
    
    if (rectog<0.0001) 
        np=1;
    end
    
    can = pi/resnog;
    
    a = zeros(1,ns);
    
    for k=1:ns       
        xx = k-1-hwid;
        for jp = -np+1:2:np-1
            t  = can*(xx+jp*del/2);
            t2 = t*t;
            t4 = t2*t2;
            if (t2 >= 1.2) 
                q0 = sin(t)/t;
                p  = cos(t);
                q1 = 3.0*(q0-p)/t2;
                tr = 2.0*(1-p)/t2;
                q2 = -(15.0*((1.0-3.0/t2)*q0+3.0*p/t2)/t2);
                q4 = 945.0*((1.0-45.0/t2+105.0/t4)*q0+5.0*(2.0-21.0/t2)*p/t2)/t4;
            else
                t6 = t2*t4;
                t8 = t2*t6;
                q0 = 1.0 - t2/6.0  + t4/120.0  - t6/5040.0   + t8/362880.0;
                q1 = 1.0 - t2/10.0 + t4/280.0  - t6/15120.0  + t8/1330560.0;
                tr = 1.0 - t2/12.0 + t4/360.0  - t6/20160.0  + t8/1814400.0;
                q2 = 1.0 - t2/14.0 + t4/504.0  - t6/33264.0  + t8/3459456.0;
                q4 = 1.0 - t2/22.0 + t4/1144.0 - t6/102960.0 + t8/14002560.0;
            end
            
            if(apo>=4) 
                a(k) = a(k) + tr;
            elseif (apo==0) 
                a(k) = a(k) + Q0;
            else
                a(k) = a(k) + c(apo,1)*q0 + c(apo,2)*q1 + c(apo,3)*q2 + c(apo,4)*q4;
            end
            
        end 
      
        a(k) = a(k) * ((1.0-(xx/(hwid+0.0))^2)^2);
        
    end 
    
vngrid = freq-(ns-1)/2*grid_step:grid_step:freq+(ns-1)/2*grid_step;

%ai = interp1(wn2wl(vngrid),a,wl_grid);

%ai(isnan(ai)==1) = 0;

%a = a';
%ai = ai';
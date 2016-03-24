function [p1,p2] = crossing_point(v,vec,re,alt)
% [p1,p2] = crossing_point(v,vec,re,alt)
%
% crossing of vec and the layer (re+alt)

% angle between v and vec
theta = acos(dot(-v,vec)/(norm(-v)*norm(vec)));

% crossing of vec and the layer (re+alt)
a = norm(v);
c = re+alt;
t = roots([1,-2*a*cos(theta),(a^2-c^2)]);
p1 = v + min(t)*vec; 
p2 = v + max(t)*vec;


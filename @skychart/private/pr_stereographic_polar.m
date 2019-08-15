function [X,Y]=pr_stereographic_polar(Az,ZenithDist);
%------------------------------------------------------------------------------
% pr_stereographic_polar function                                     AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              Stereographic polar projection.
%              This projection preserves angles.
% Input  : - Vector of Azimuth, in deg.
%          - Vector of Zenith-distance, in deg.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                  November 2004  
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------
  d2r   = pi/180;
  Az    = Az*d2r; ZenithDist=ZenithDist*d2r;
  X     = cos(Az).*tan(0.5.*ZenithDist);
  Y     = sin(Az).*tan(0.5.*ZenithDist);
end % pr_stereographic_polar

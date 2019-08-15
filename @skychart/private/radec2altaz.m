function [Az Alt] = radec2altaz(RA,DEC, julianday, place)
  % radec2altaz: convert RA/DEC coordinates to horizontal frame
  %
  % [az,alt] = radec2altaz(ra, dec, julianday, place): compute the Alt Az 
  %   (horizontal coordinates) from the RA DEC ones, all in [deg].
  %
  % DEC is the      'height' respective to North pole
  % Altitude is the 'height' respective to Zenith
  %
  % uses: horiz_coo, refraction (MAAT by E. Ofek)
  %       URL : http://wise-obs.tau.ac.il/~eran/matlab.html
  
  d2r = pi/180;
  r2d = 1/d2r;

  % taken from MAAT: horiz_coo and refraction (take [rad])
  % Eran O. Ofek http://weizmann.ac.il/home/eofek/matlab/
  HorizCoo          = horiz_coo([ RA(:) DEC(:) ]*d2r, ...
    julianday, place*d2r, 'h');
  HorizCoo(:,2)     = HorizCoo(:,2) + refraction(HorizCoo(:,2));
  HorizCoo          = HorizCoo*r2d;
  Az                = HorizCoo(:,1);
  Alt               = 90-HorizCoo(:,2);

end % radec2altaz

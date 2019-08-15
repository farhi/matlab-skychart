function LST=lst(JD,EastLong,dummy)
  %--------------------------------------------------------------------------
  % lst function                                                       ephem
  % Description: Local Sidereal Time, (mean or apparent), for vector of
  %              JDs and a given East Longitude.
  % Input  : - Vector of JD [days], in UT1 time scale.
  %          - East Longitude in radians.
  % Output : - vector of LST in fraction of day.
  % Tested : Matlab 5.3
  %     By : Eran O. Ofek                    Aug 1999
  %    URL : http://weizmann.ac.il/home/eofek/matlab/
  % Example: LST=lst(2451545+[0:1:5]',0);  % LST at Greenwhich 0 UT1
  % Reliable: 1
  %--------------------------------------------------------------------------

    % convert JD to integer day + fraction of day
    TJD     = floor(JD - 0.5) + 0.5;
    DayFrac = JD - TJD;
    T       = (TJD - 2451545.0)./36525.0;
    GMST0UT = 24110.54841 + 8640184.812866.*T + 0.093104.*T.*T - 6.2e-6.*T.*T.*T;

    % convert to fraction of day in range [0 1)
    GMST0UT = GMST0UT./86400.0;
    GMST0UT = GMST0UT - floor(GMST0UT);
    LST     = GMST0UT + 1.0027379093.*DayFrac + EastLong./(2.*pi);
    LST     = LST - floor(LST);
end % lst

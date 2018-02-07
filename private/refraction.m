function [R]=refraction(Alt,MetoData,Formula)
%------------------------------------------------------------------------------
% refraction function                                                    ephem
% Description: Estimate atmospheric refraction, in visible light.
% Input  : - Vector of altitude [radians].
% Output : - Refraction correction in radians. (Add to Alt).
% Reference : AA 2001
% Tested : Matlab 5.3
%     By : Eran O. Ofek                   Nov 2000
%    URL : http://weizmann.ac.il/home/eran/matlab/
% See also: refraction_wave.m
% Example: [R]=refraction([20;30]./RAD);
% Reliable: 1
%------------------------------------------------------------------------------
  RAD       = 180./pi;
  MetoData  = [20, 1000, 0.8];
  AltD      = Alt.*RAD;

  T   = MetoData(1) + 273;  % convert to K
  P   = MetoData(2);
  Ia  = find(AltD>15);
  Ib  = find(AltD<=15 & AltD>-0.6);
  In  = find(AltD<=-0.6);
  R     = zeros(size(Alt));
  R(Ia) = 0.00452.*P./(T.*tan(Alt(Ia)).*RAD);
  R(Ib) = P.*(0.1594 + 0.0196.*AltD(Ib) + 0.00002.*AltD(Ib).*AltD(Ib)) ...
          ./(T.*(1 + 0.505.*AltD(Ib) + 0.0845.*AltD(Ib).*AltD(Ib)).*RAD);
  R(In) = 0;
end % refraction

function [OutData]=convertdms(InData,InType,OutType)
%--------------------------------------------------------------------------
% convertdms function                                                ephem
% Description: Convert between various representations of coordinates and
%              time as sexagesimal coordinates, degrees and radians.
% Input  : - Input data in one of the following formats:
%            degrees, radians, [H M S.S], [Sign D M S.S],
%            [H M.M], [Sign D M.M], 
%            'HH:MM:SS.S', '+DD:MM:SS.S'.
%          - Input type:
%            'D'  - [Sign D M S.S] format.
%            'H'  - [H M S.S] format.
%            'f'  - fraction in range [0,1].
%          - Output type:
%            'r'  - radians, in range [0,2*pi] (default).
%            'R'  - radians, in range [-pi,pi].
% Output : - Requested output.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jun 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: convertdms([1;2;3],'H','d');
% Reliable: 2
%--------------------------------------------------------------------------
  RAD     = 180./pi;
  Epsilon = 1e-8;
  IsCell  = 0;
  SignVec = ['-';'+'];

  if (nargin==2),
     OutType = 'r';
  elseif (nargin==3),
     % no default
  else
     error('Illegal number of input arguments');
  end

  % convert InData to radians
  switch InType
  case {'D'}
    RadData = (InData(:,2)+InData(:,3)./60+InData(:,4)./3600).*InData(:,1)./RAD;
  case {'H'}
    RadData = (InData(:,1)+InData(:,2)./60+InData(:,3)./3600).*15./RAD;
  case {'f'}
    RadData = InData.*2.*pi;
  otherwise
    error([ 'Illegal InType ' InType ]);
  end

  switch OutType
  case {'r'}
    % do nothing, allready in range [0,2*pi]
    OutData = RadData;
  case {'R'}
    RadData = 2.*pi.*(RadData./(2.*pi) - floor(RadData./(2.*pi)));
    I = find(RadData>pi);
    OutData = RadData;
    OutData(I) = RadData(I) - 2.*pi;
  case {'f'}
    OutData = RadData./(2.*pi);
  otherwise
    error([ 'Illegal OutType ' OutType ]);
  end
end % convertdms

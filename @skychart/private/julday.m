function JD=julday(Date,Output)
%--------------------------------------------------------------------------
% julday function                                                    ephem 
% Description: Convert Julian/Gregorian date to Julian Day.
% Input  : - Gregorian of Julian date in one of the following formats
%            [Day, Month, Year, Day_Fraction]
%            or [Day, Month, Year, Hour, Min, Sec]
%            or [Day, Month, Year] - in this case set Day_Fraction to 0.
% Output : - Row vector of Julian days.
% Tested : Matlab 3.5
%     By : Eran O. Ofek                    Jan 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: julday1.m, jd2date.m, date_axis.m
% Example: julday([1 1 2000 10 30 0]);
%          julday([1 1 2000; 2 2 2000]);
%          julday;  % JD of now
% Reliable: 1
%--------------------------------------------------------------------------

  if (nargin==0),
     Date = clock;
     Date = Date(:,[3 2 1 4 5 6]);
  end
  if (isempty(Date)),
     Date = clock;
     Date = Date(:,[3 2 1 4 5 6]);
  end

  Y = Date(:,3);
  M = Date(:,2);
  D = Date(:,1);

  [Lines, Rows] = size(Date);
  switch Rows
   case 3
      F = zeros(Lines,1);
   case 4
      F = Date(:,4);
   case 6
      F = convertdms(Date(:,4:6),'H','f');
   otherwise
      error('Illegal number of column in Date matrix');
  end

  B       = zeros(Lines,1);
  Im3     = find(M<3);
  Y(Im3)  = Y(Im3) - 1;
  M(Im3)  = M(Im3) + 12;

  Iy    = find(Y>1582 | (Y==1582 & M>10) | (Y==1582 & M==10 & D>=15));
  A     = floor(Y(Iy).*0.01);
  B(Iy) = 2 - A + floor(0.25.*A);
  JD    = floor(365.25.*(Y + 4716)) + floor(30.6001.*(M + 1)) + D + B - 1524.5 + F;

end % julday

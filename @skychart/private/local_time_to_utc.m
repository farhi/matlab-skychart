function dateStringUtc = local_time_to_utc( dateVectorLocal, format )
%
% local_time_to_utc -- convert local time/date information to UTC.
%
% dateStringUtc = local_time_to_utc( dateVectorLocal ) converts a vector of local
%    date/time information to UTC dates and times.  The format of dateVectorLocal can be
%    either Matlab datenumbers (N x 1 or 1 x N) or else a vector of Matlab date vectors
%    (ie, N x 6, for example [2008 10 22 13 41 00 ; 2008 10 22 13 41 12]).  The time/date
%    values will be converted to UTC and returned in 'dd-mmm-yyyy HH:MM:SS Z' format.
%
% dateStringUtc = local_time_to_utc( dateVectorLocal, format ) returns the UTC time / date
%    values in a format specified by numeric argument format.  The format values match
%    those used by the Matlab datestr function.  Supported formats are as follows:
%
%    Number           String                     Example
%    ===========================================================================
%       0             'dd-mmm-yyyy HH:MM:SS Z'   01-Mar-2000 15:45:17 Z
%      13             'HH:MM:SS Z'               15:45:17 Z    
%      15             'HH:MM Z'                  15:45 Z       
%      21             'mmm.dd,yyyy HH:MM:SS Z'   Mar.01,2000 15:45:17 Z
%      30 (ISO 8601)  'yyyymmddTHHMMSSZ'         20000301T154517Z 
%      31             'yyyy-mm-dd HH:MM:SS Z'    2000-03-01 15:45:17 Z
%
% Version date:  2008-November-10.
%

%Copyright (c) 2009, Peter Tenenbaum
%All rights reserved.

%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:

%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%    * Neither the name of the National Aeronautics and Space Administration nor the names
%      of its contributors may be used to endorse or promote products derived
%      from this software without specific prior written permission.

% https://fr.mathworks.com/matlabcentral/fileexchange/22295-local-time-to-utc
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.

% Modification History:
%
%    2008-November-12, PT:
%        switch to use of Java methods and eliminate check for unix vs PC (java methods
%        work on all platforms).
%    2008-November-10, PT:
%        added a bunch of disp statements which should help diagnose mysterious smoke test
%        failures in PMD and PA.
%
%=========================================================================================

% arguments -- if the argument is a vector of datenums, convert to a vector of time/date
% vectors

  if ( size(dateVectorLocal,2) ~= 6 || max(dateVectorLocal(:,2)) > 12 )
      dateVectorLocal = datevec(dateVectorLocal) ;
  end
  
% check that the format number is supported

  if ( nargin == 1)
      format = 0 ;
  end
  supportedFormats = [0 13 15 21 30 31] ;
  if ~ismember(format,supportedFormats)
      error('common:localTimeToUtc:unsupportedFormat' , ...
          'local_time_to_utc:  unsupported format requested' ) ;
  end
  
% call the converter vector

  dateVectorUtc = local_date_vector_to_utc( dateVectorLocal ) ;
  
% convert the date vectors to the desired format

  dateStringUtc = datestr(dateVectorUtc,format) ;
  stringLength = size(dateStringUtc,2) ;
  nStrings = size(dateStringUtc,1) ;
  
% append the 'Z' 

  dateStringUtc = append_z_to_date_string( dateStringUtc, format ) ;
  
return  
  
% and that's it!

%
%
%

%=========================================================================================

% subfunction which converts a Matlab date vector in local time to one in UTC

function dateVectorUtc = local_date_vector_to_utc( dateVectorLocal )
      
  nDates = size(dateVectorLocal,1) ;
  dateVectorUtc = zeros(nDates,6) ;
  
% import the Java classes needed for this process

  import java.text.SimpleDateFormat ;
  import java.util.Date ;
  import java.util.TimeZone ;
  
% instantiate a SimpleDateFormat object with a fixed time/date format and UTC time zone

  utcFormatObject = SimpleDateFormat('yyyy-MM-dd HH:mm:ss') ;
  utcFormatObject.setTimeZone(TimeZone.getTimeZone('UTC')) ;
      
% loop over date strings

  for iDate = 1:nDates

      dateVec = dateVectorLocal(iDate,:) ;
      
%     instantiate a Java Date class object with the local time.  Note that Java year is
%     year since 1900, and Java month is zero-based
      
      localDateObject = Date(dateVec(1)-1900, dateVec(2)-1, dateVec(3), ...
                             dateVec(4), dateVec(5), dateVec(6)) ;                        
                         
%     convert the date object to a string in the correct format and in UTC

      dateStringUtc = char(utcFormatObject.format(localDateObject)) ;
         
%     pick through the resulting string and extract the data we want, converting to
%     numbers as we go

      dateVectorUtc(iDate,1) = str2num(dateStringUtc(1:4)) ;
      dateVectorUtc(iDate,2) = str2num(dateStringUtc(6:7)) ;
      dateVectorUtc(iDate,3) = str2num(dateStringUtc(9:10)) ;
      dateVectorUtc(iDate,4) = str2num(dateStringUtc(12:13)) ;
      dateVectorUtc(iDate,5) = str2num(dateStringUtc(15:16)) ;
      dateVectorUtc(iDate,6) = str2num(dateStringUtc(18:19)) ;
          
  end % loop over dates
      
return

%=========================================================================================

% append a 'Z', signifying UTC, to a Matlab date string.

function datestrZ = append_z_to_date_string( datestr, format )

  nStrings = size(datestr,1) ;
  stringLength = size(datestr,2) ;
  if (format ~= 30)
      addString = repmat(' Z',nStrings,1) ;
      stringSizeIncrease = 2 ;
  else
      addString = repmat('Z',nStrings,1) ;
      stringSizeIncrease = 1 ;
  end
  datestrZ = [datestr(:) ; addString(:)] ;
  datestrZ = reshape(datestrZ,nStrings,stringLength+stringSizeIncrease) ;

return

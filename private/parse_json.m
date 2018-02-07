function data = parse_json(string)
% DATA = PARSE_JSON(string)
% This function parses a JSON string and returns a cell array with the
% parsed data. JSON objects are converted to structures and JSON arrays are
% converted to cell arrays.

% F. Glineur, 2009
% (inspired by the JSON parser by Joel Feenstra on MATLAB File Exchange
% (http://www.mathworks.com/matlabcentral/fileexchange/20565) but with 
% faster handling of strings)

% https://fr.mathworks.com/matlabcentral/fileexchange/23393--another--json-parser
% Copyright (c) 2011, François Glineur
%  All rights reserved.

%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:

%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the distribution

%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

pos = 1;
len = length(string);
% String delimiters and escape characters are identified beforehand to improve speed
esc = regexp(string, '["\\]'); index_esc = 1; len_esc = length(esc);

if pos <= len
    switch(next_char)
        case '{'
            data = parse_object;
        case '['
            data = parse_array;
        otherwise
            error_pos('Outer level structure must be an object or an array');
    end
end

    function object = parse_object
        parse_char('{');
        object = [];
        if next_char ~= '}'
            while 1
                str = parse_string;
                if isempty(str)
                    error_pos('Name of value at position %d cannot be empty');
                end
                parse_char(':');
                val = parse_value;
                object.(valid_field(str)) = val;
                if next_char == '}'
                    break;
                end
                parse_char(',');
            end
        end
        parse_char('}');
    end

    function object = parse_array
        parse_char('[');
        object = cell(0, 1);
        if next_char ~= ']'
            while 1
                val = parse_value;
                object{end+1} = val;
                if next_char == ']'
                    break;
                end
                parse_char(',');
            end
        end
        parse_char(']');
    end

    function parse_char(c)
        skip_whitespace;
        if pos > len || string(pos) ~= c
            error_pos(sprintf('Expected %c at position %%d', c));
        else
            pos = pos + 1;
            skip_whitespace;
        end
    end

    function c = next_char
        skip_whitespace;
        if pos > len
            c = [];
        else
            c = string(pos);
        end        
    end
    
    function skip_whitespace
        while pos <= len && isspace(string(pos))
            pos = pos + 1;
        end
    end

     function str = parse_string
        if string(pos) ~= '"'
            error_pos('String starting with " expected at position %d');
        else
            pos = pos + 1;
        end
        str = '';
        while pos <= len
            while index_esc <= len_esc && esc(index_esc) < pos 
                index_esc = index_esc + 1;
            end
            if index_esc > len_esc
                str = [str string(pos:end)];
                pos = len + 1;
                break;
            else
                str = [str string(pos:esc(index_esc)-1)];
                pos = esc(index_esc);
            end
            switch string(pos)
                case '"' 
                    pos = pos + 1;
                    return;
                case '\'
                    if pos+1 > len
                        error_pos('End of file reached right after escape character');
                    end
                    pos = pos + 1;
                    switch string(pos)
                        case {'"' '\' '/'}
                            str(end+1) = string(pos);
                            pos = pos + 1;
                        case {'b' 'f' 'n' 'r' 't'}
                            str(end+1) = sprintf(['\' string(pos)]);
                            pos = pos + 1;
                        case 'u'
                            if pos+4 > len
                                error_pos('End of file reached in escaped unicode character');
                            end
                            str(end+1:end+6) = string(pos-1:pos+4);
                            pos = pos + 5;
                    end
                otherwise % should never happen
                    str(end+1) = string(pos);
                    pos = pos + 1;
            end
        end
        error_pos('End of file while expecting end of string');
    end

    function num = parse_number
        [num, one, err, delta] = sscanf(string(pos:min(len,pos+20)), '%f', 1); % TODO : compare with json(pos:end)
        if ~isempty(err)
            error_pos('Error reading number at position %d');
        end
        pos = pos + delta-1;
    end

    function val = parse_value
        switch(string(pos))
            case '"'
                val = parse_string;
                return;
            case '['
                val = parse_array;
                return;
            case '{'
                val = parse_object;
                return;
            case {'-','0','1','2','3','4','5','6','7','8','9'}
                val = parse_number;
                return;
            case 't'
                if pos+3 <= len && strcmpi(string(pos:pos+3), 'true')
                    val = true;
                    pos = pos + 4;
                    return;
                end
            case 'f'
                if pos+4 <= len && strcmpi(string(pos:pos+4), 'false')
                    val = false;
                    pos = pos + 5;
                    return;
                end
            case 'n'
                if pos+3 <= len && strcmpi(string(pos:pos+3), 'null')
                    val = [];
                    pos = pos + 4;
                    return;
                end
        end
        error_pos('Value expected at position %d');
    end

    function error_pos(msg)
        poss = max(min([pos-15 pos-1 pos pos+20],len),1);
        if poss(3) == poss(2)
            poss(3:4) = poss(2)+[0 -1];         % display nothing after
        end
        msg = [sprintf(msg, pos) ' : ... ' string(poss(1):poss(2)) '<error>' string(poss(3):poss(4)) ' ... '];
        ME = MException('JSONparser:invalidFormat', msg);
        throw(ME);
    end

    function str = valid_field(str)   
    % From MATLAB doc: field names must begin with a letter, which may be
    % followed by any combination of letters, digits, and underscores.
    % Invalid characters will be converted to underscores, and the prefix
    % "alpha_" will be added if first character is not a letter.
        if ~isletter(str(1))
            str = ['alpha_' str];
        end
        str(~isletter(str) & ~('0' <= str & str <= '9')) = '_';   
    end

end

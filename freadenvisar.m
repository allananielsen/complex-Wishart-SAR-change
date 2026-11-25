function [image,p,t,hdr_size] = freadenvisar(fname, varargin)

% freadenvisar - read fully polarimetric SAR image with ENVI like header
%       (V. Guissard, Apr 29, 2004)
%
% 				Reads an image of ENVI standard type 
%				to a [(lines x samples) x bands] MATLAB array
%               adapted to polarimetric SAR images
%               (by Allan A. Nielsen, 23 Dec 2011)
%
% Syntax
%
% image = freadenvisar(fname);
% [image,p] = freadenvisar(fname);
% [image,p,t] = freadenvisar(fname,'b',0);
%
% INPUT :
%
% fname	string	giving the full pathname of the ENVI image to read
% 'b' or 'l' for big or little endian (defaults to 'n')
% read_hdr_only (defaults to false)
% hdr_size
%
% OUTPUT :
%
% image	lines by samples by bands array containing the ENVI image values organised in
%				(lines x samples) x bands.
%             output order is alphabetical in indexes with
%             real parts before imaginary parts, fx. for full polarimetry:
%               hhhh (real), hhhv (complex), hhvv (complex),
%               hvhv (real), hvvv (complex), vvvv (real)
% p		1 by 3	vector that contains (1) the number of samples, (2) the
%       number of liness and (3) the number of bands of the opened image.
%
% t		string	describing the image data type in MATLAB convention.
%
% NOTE : freadenvi needs the corresponding image header file generated
%		 automatically by ENVI. The ENVI header file must have the same name
%		 as the ENVI image file + the '.hdr' extension.
% Example hdr-file (interleave *must* be 'bsq':
% samples = 1024
% lines   = 1024
% bands   = 6
% data type = 4
% interleave = bsq
%
% Some of the first files I worked with were written row-major, big endian
% Matlab reads/writes column-major.
% There may be issues with transpose/permute, #rows/#cols.
%
% After reading the hdr-file '*.hdr' we first try to open a file called fx '*hhhh'.
% If that file does not exist we try to open a file called '*hhhh.emi'.
% This file contains in its first line a header which we skip.

% Modified for polarimetric SAR data by
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk
% 23 Dec 2011 and 27 May 2016

endian = 'n'; % 'native' to your system
if nargin>1
%     endian = 'b'; % big endian: 'b' (typically UNIX)
%     endian = 'l'; % little endian: 'b' (PC)
    endian = varargin{1}; % endian big: 'b' (UNIX)
end
% hdr_size = 4*1024; % 4*2048
read_hdr_only = false;
if nargin>2
%     endian = 'b'; % big endian: 'b' (typically UNIX)
%     endian = 'l'; % little endian: 'b' (PC)
    read_hdr_only = varargin{2};
end

% Parameters initialization
elements = {'samples ' 'lines   ' 'bands   ' 'data type ' 'interleave '};
d = {'uint8' 'int16' 'int32' 'float32' 'float64' 'uint16' 'uint32' 'int64' 'uint64'};
%d = {'*uint8' 'int16' 'int32' 'float32' 'float64' 'uint16' 'uint32' 'int64' 'uint64'};
% Check user input
if ~ischar(fname)
    error('fname should be a char string');
end

% Open ENVI header file to retreive s, l, b & d variables
rfid = fopen(strcat(fname,'.hdr'),'r');

% Check if the header file is correctely opened
if rfid == -1
    error('Input header file does not exist');
end;

% Read ENVI image header file and get p(1) : nb samples,
% p(2) : nb lines, p(3) : nb bands and t : data type
while 1
    tline = fgetl(rfid);
    if ~ischar(tline), break, end
    [first,second] = strtok(tline,'=');
    
    switch first
        case elements(1)
            [f,s] = strtok(second);
            p(1) = str2num(s);
        case elements(2)
            [f,s] = strtok(second);
            p(2) = str2num(s);
        case elements(3)
            [f,s] = strtok(second);
            p(3) = str2num(s);
        case elements(4)
            [f,s] = strtok(second);
            t = str2num(s);
            switch t
                case 1
                    t = d(1);
                case 2
                    t = d(2);
                case 3
                    t = d(3);
                case 4
                    t = d(4);
                case 5
                    t = d(5);
                case 12
                    t = d(6);
                case 13
                    t = d(7);
                case 14
                    t = d(8);
                case 15
                    t = d(9);
                otherwise
                    error('Unknown image data type');
            end
        case elements(5) % input must be BSQ here
            [f,interleave] = strtok(second);
            if ~strcmp(interleave,' bsq')
                error('''interleave'' in .hdr file  must be '' bsq''')
            end
    end
end
fclose(rfid); 

hdr_size = 4*p(1); %4*1024; % 4*2048
if nargin>3
    hdr_size = varargin{3};
end

t = t{1,1};

if read_hdr_only
    image = NaN;
    return
end

% Open the ENVI image and store it in the 'image' MATLAB array
disp(['Input is (',(num2str(p(2))),' lines) x (', ...
        (num2str(p(1))),' samples) x (',(num2str(p(3))),' bands)' ...
        ' type ', (t), (' image ...')]);

if p(3)==6 % full polarimetric matrix
    
    % read hhhh (real)
    % fid = fopen(strcat(fname,'hhhh'),'r');
    fid = fopen(strcat(fname,'hhhh'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hhhh.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hhhh.emi',' not found'));
        end
        fseek(fid,hdr_size,-1);
    end
    [hhhh,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1))
        error('data in fname do not match nrows, ncols, nvars');
    end
    % hhhh = reshape(hhhh,[p(1),p(2),p(3)]); % original freadenvi reshape
    hhhh = reshape(hhhh,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);

    % read hhhv (complex)
    % fid = fopen(strcat(fname,'hhhv'),'r');
    fid = fopen(strcat(fname,'hhhv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hhhv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hhhv.emi',' not found'));
        end
        fseek(fid,2*hdr_size,-1);
    end
    [hhhv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1)*2)
        error('data in fname do not match nrows, ncols, nvars');
    end
    hhhvRE = hhhv(1:2:end-1);
    hhhvIM = hhhv(2:2:end);
    clear hhhv
    hhhvRE = reshape(hhhvRE,[p(1),p(2)]); % original freadenvi reshape
    hhhvIM = reshape(hhhvIM,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);

    % read hhvv (complex)
    % fid = fopen(strcat(fname,'hhvv'),'r');
    fid = fopen(strcat(fname,'hhvv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hhvv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hhvv.emi',' not found'));
        end
        fseek(fid,2*hdr_size,-1);
    end
    [hhvv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1)*2)
        error('data in fname do not match nrows, ncols, nvars');
    end
    hhvvRE = hhvv(1:2:end-1);
    hhvvIM = hhvv(2:2:end);
    clear hhvv
    hhvvRE = reshape(hhvvRE,[p(1),p(2)]); % original freadenvi reshape
    hhvvIM = reshape(hhvvIM,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);

    % read hvhv (real)
    % fid = fopen(strcat(fname,'hvhv'),'r');
    fid = fopen(strcat(fname,'hvhv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hvhv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hvhv.emi',' not found'));
        end
        fseek(fid,hdr_size,-1);
    end
    [hvhv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1))
        error('data in fname do not match nrows, ncols, nvars');
    end
    % hvhv = reshape(hvhv,[p(1),p(2),p(3)]); % original freadenvi reshape
    hvhv = reshape(hvhv,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);

    % read hvvv (complex)
    % fid = fopen(strcat(fname,'hvvv'),'r');
    fid = fopen(strcat(fname,'hvvv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hvvv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hvvv.emi',' not found'));
        end
        fseek(fid,2*hdr_size,-1);
    end
    [hvvv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1)*2)
        error('data in fname do not match nrows, ncols, nvars');
    end
    hvvvRE = hvvv(1:2:end-1);
    hvvvIM = hvvv(2:2:end);
    clear hvvv
    hvvvRE = reshape(hvvvRE,[p(1),p(2)]); % original freadenvi reshape
    hvvvIM = reshape(hvvvIM,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);

    % read vvvv (real)
    % fid = fopen(strcat(fname,'vvvv'),'r');
    fid = fopen(strcat(fname,'vvvv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'vvvv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'vvvv.emi',' not found'));
        end
        fseek(fid,hdr_size,-1);
    end
    [vvvv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1))
        error('data in fname do not match nrows, ncols, nvars');
    end
    % vvvv = reshape(vvvv,[p(1),p(2),p(3)]); % original freadenvi reshape
    vvvv = reshape(vvvv,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);

    image = cat(3,hhhh,hhhvRE,hhhvIM,hhvvRE,hhvvIM,hvhv,hvvvRE,hvvvIM,vvvv);
    % imshowsarrgb(image,[6 1 9])
    
elseif p(3)==4 % azimuthal symmetry
    
    % read hhhh (real)
    % fid = fopen(strcat(fname,'hhhh'),'r');
    fid = fopen(strcat(fname,'hhhh'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hhhh.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hhhh.emi',' not found'));
        end
        fseek(fid,hdr_size,-1);
    end
    [hhhh,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1))
        error('data in fname do not match nrows, ncols, nvars');
    end
    % hhhh = reshape(hhhh,[p(1),p(2),p(3)]); % original freadenvi reshape
    hhhh = reshape(hhhh,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);
    
    % read hhvv (complex)
    % fid = fopen(strcat(fname,'hhvv'),'r');
    fid = fopen(strcat(fname,'hhvv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hhvv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hhvv.emi',' not found'));
        end
        fseek(fid,2*hdr_size,-1);
    end
    [hhvv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1)*2)
        error('data in fname do not match nrows, ncols, nvars');
    end
    hhvvRE = hhvv(1:2:end-1);
    hhvvIM = hhvv(2:2:end);
    clear hhvv
    hhvvRE = reshape(hhvvRE,[p(1),p(2)]); % original freadenvi reshape
    hhvvIM = reshape(hhvvIM,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);
    
    % read hvhv (real)
    % fid = fopen(strcat(fname,'hvhv'),'r');
    fid = fopen(strcat(fname,'hvhv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hvhv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hvhv.emi',' not found'));
        end
        fseek(fid,hdr_size,-1);
    end
    [hvhv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1))
        error('data in fname do not match nrows, ncols, nvars');
    end
    % hvhv = reshape(hvhv,[p(1),p(2),p(3)]); % original freadenvi reshape
    hvhv = reshape(hvhv,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);
    
    % read vvvv (real)
    % fid = fopen(strcat(fname,'vvvv'),'r');
    fid = fopen(strcat(fname,'vvvv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'vvvv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'vvvv.emi',' not found'));
        end
        fseek(fid,hdr_size,-1);
    end
    [vvvv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1))
        error('data in fname do not match nrows, ncols, nvars');
    end
    % vvvv = reshape(vvvv,[p(1),p(2),p(3)]); % original freadenvi reshape
    vvvv = reshape(vvvv,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);
    
    image = cat(3,hhhh,hhvvRE,hhvvIM,hvhv,vvvv);
    % imshowsarrgb(image,[4 1 5])
    
elseif p(3)==3 % dual polarimetry, here HH/HV
    
    % read hhhh (real)
    % fid = fopen(strcat(fname,'hhhh'),'r');
    fid = fopen(strcat(fname,'hhhh'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hhhh.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hhhh.emi',' not found'));
        end
        fseek(fid,hdr_size,-1);
    end
    [hhhh,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1))
        error('data in fname do not match nrows, ncols, nvars');
    end
    % hhhh = reshape(hhhh,[p(1),p(2),p(3)]); % original freadenvi reshape
    hhhh = reshape(hhhh,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);
    
    % read hhvv (complex)
    % fid = fopen(strcat(fname,'hhvv'),'r');
    fid = fopen(strcat(fname,'hhvv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'hhvv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'hhvv.emi',' not found'));
        end
        fseek(fid,2*hdr_size,-1);
    end
    [hhvv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1)*2)
        error('data in fname do not match nrows, ncols, nvars');
    end
    hhvvRE = hhvv(1:2:end-1);
    hhvvIM = hhvv(2:2:end);
    clear hhvv
    hhvvRE = reshape(hhvvRE,[p(1),p(2)]); % original freadenvi reshape
    hhvvIM = reshape(hhvvIM,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);
    
    % read vvvv (real)
    % fid = fopen(strcat(fname,'vvvv'),'r');
    fid = fopen(strcat(fname,'vvvv'),'r',endian);
    if fid==-1
        fid = fopen(strcat(fname,'vvvv.emi'),'r',endian);
        if fid==-1
            error(strcat(fname,'vvvv.emi',' not found'));
        end
        fseek(fid,hdr_size,-1);
    end
    [vvvv,count] = fread(fid,t);
    % if count~=(p(2)*p(1)*p(3))
    if count~=(p(2)*p(1))
        error('data in fname do not match nrows, ncols, nvars');
    end
    % vvvv = reshape(vvvv,[p(1),p(2),p(3)]); % original freadenvi reshape
    vvvv = reshape(vvvv,[p(1),p(2)]); % original freadenvi reshape
    fclose(fid);
    
    image = cat(3,hhhh,hhvvRE,hhvvIM,vvvv);
    % imshowsarrgb(image,[4 1 5])
    
else
    
    image = NaN;
    error('freadenvisar: input must contain 6, 4 or 3 real and complex frames')

end

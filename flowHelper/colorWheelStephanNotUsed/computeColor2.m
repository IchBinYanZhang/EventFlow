function img = computeColor2(u,v, opt)
% img = computeColor2(u,v, opt) - colory images for flow field
%
% INPUT
%   u,v             flow field u,v components
%   opt.ordered     structured or random arrangement
%   opt.overlap     allow overlapping of arrows,lines
%   opt.arrowtype   0=none, 1=arrows, 2=lines
%   opt.speedfak    speed gain
%   opt.usecolor    color code 0=off 1=on
%   opt.maxspeed=1; 
% OUTPUT
%   img[h,w,3]
%
% EXAMPLE
% opt=struct();
% ccode2 = computeColor2(u,v,opt);
%
% (2010.08)

%
%   According to the c++ source code of Daniel Scharstein
%   Contact: schar@middlebury.edu
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-10-31 21:20:30 (Wed, 31 Oct 2006) $

% Copyright 2007, Deqing Sun.
%
%                         All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author and Brown University not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
%
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
% PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR BROWN UNIVERSITY BE LIABLE FOR
% ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

if nargin()<3, opt=struct(); end

[h,w]=size(u);

speed = sqrt( u.^2 + v.^2 );
maxspeed = max(speed(:));
if (maxspeed==0) maxspeed=1; end;

if isfield(opt,'arrowtype'), arrowtype=opt.arrowtype; else arrowtype=0; end
if isfield(opt,'usecolor'), usecolor=opt.usecolor; else usecolor=1; end

if isfield(opt,'maxspeed'), maxspeed=opt.maxspeed; else maxspeed=1; end

if isfield(opt,'overlap'), overlap=opt.overlap; else overlap=1; end
if isfield(opt,'ordered'), ordered=opt.ordered; else ordered=1; end
%if isfield(opt,'arrows'), arrows=opt.arrows; else arrows=1; end

% for arrows 
if isfield(opt,'speedfak'), speedfak=opt.speedfak; else speedfak=10; end
if isfield(opt,'stepsize'), stepsize=opt.stepsize; else stepsize=speedfak*2; end

if isfield(opt,'noise'), noise=opt.noise; else noise=0; end
if isfield(opt,'active'), active=opt.active; else active=1; end

if active<1,
    disp('WARNING: computeColor2() is OFFLINE. ');
    img=randn(h,w,3);
    return;
end

nanIdx = isnan(u) | isnan(v);
u(nanIdx) = 0;
v(nanIdx) = 0;

colorwheel = makeColorwheel();
ncols = size(colorwheel, 1);

%% speed normalisation

rad = sqrt(u.^2+v.^2)  ; % speed 
rad_log=real(log2(1+rad));

% rad=rad./max(rad(:));
% rad = rad ./ rad_log;
rad=rad./maxspeed;

% rad(rad<0.5850)=0;  % og2(1+0.125) % less than 1/8 pixel
rad(isnan(rad))=0;

a = atan2(-v, -u)/pi;

fk = (a+1) /2 * (ncols-1) + 1;  % -1~1 maped to 1~ncols

k0 = floor(fk);                 % 1, 2, ..., ncols

k1 = k0+1;
k1(k1==ncols+1) = 1;

f = fk - k0;




if (usecolor>0)

	for i = 1:size(colorwheel,2)
		tmp = colorwheel(:,i);
		col0 = tmp(k0)/255;
		col1 = tmp(k1)/255;
		col = (1-f).*col0 + f.*col1;
		
		
		
		idx = rad <= 1;
		col(idx) = 1-rad(idx).*(1-col(idx));    % increase saturation with radius
		
		col(~idx) = col(~idx)*1;             % out of range
		
		
		img(:,:, i) = uint8(floor(255*col.*(1-nanIdx)));
	end
else
	img = ones(h,w,3).*255;
end

usedmask = zeros(h,w);
arrowmask = zeros(h,w,3);
arrowmask1 = zeros(h,w);

itan = atan2(v,u);

noisefield = randn(h,w).*noise;

itan = itan + noisefield;

speedi = sqrt(u.^2 + v.^2);

for z1=1:stepsize:h
	for s1=1:stepsize:w
		
		if ordered==1
			z=round(z1);
			s=round(s1);
		else
			z = round(rand()*(h-1))+1;
			s = round(rand()*(w-1))+1;
		end
		
		rad = speedi(z,s)/maxspeed*speedfak;
		
		% build an arrow here
		if (usedmask(z,s)==0) && (rad>1) && (arrowtype>0)
			
			switch arrowtype
				
				case 1
					
					
					% 			dir = itan(z,s)+pi/2;
					dir = itan(z,s);
					
					py1 = round(   z + sin(dir)*rad        );
					px1 = round(   s + cos(dir)*rad        );
					
					
					
					
					py3 = round(   z + sin(dir+0.75*pi)*rad*0.5         );
					px3 = round(   s + cos(dir+0.75*pi)*rad*0.5         );
					
					py4 = round(   z + sin(dir-0.75*pi)*rad*0.5        );
					px4 = round(   s + cos(dir-0.75*pi)*rad*0.5       );
					
					
					% 		arrowmask1 = arrowmask1 + bresline([h  w], [px1 py1; px2 py2] );
					arrowmask1 = arrowmask1 + bresline([h  w], [px3 py3; px1 py1] );
					arrowmask1 = arrowmask1 + bresline([h  w], [px4 py4; px1 py1] );
					arrowmask1 = arrowmask1 + bresline([h  w], [px4 py4; px3 py3] );
					
					
					
				case 2
					
					% tangential
					
					speed=speedi(z,s);
					dir = itan(z,s);
					
					py1 = round(   z + sin(dir)*rad        );
					px1 = round(   s + cos(dir)*rad        );
					
					py2 = round(   z + sin(dir+3.14)*rad      );
					px2 = round(   s + cos(dir+3.14)*rad        );
					
					arrowmask1 = arrowmask1 + bresline([h  w], [px1 py1; px2 py2] );
					
				otherwise
					% none
					
			end
			
			% 			min([py1 py2 py3 py4]):max([py1 py2 py3 py4]), min([px1 px2 px3 px4]):max([px1 px2 px3 px4]))=1;
			
			range = max(5, round(rad*1.2) );
			
			% do not paint arrows over others
			if overlap<1
				usedmask(max(1, z-range):min(z+range,h),  max(1, s-range):min(s+range,w)) = 1;
			end
		end
		
		% occupy used-mask
		
		% 		arrowmask(5,5:10,:)=1;
	end
end

arrowmask = repmat(arrowmask1,[1 1 3]);

img(arrowmask>0)=0;


end



% %%
% function colorwheel = makeColorwheel()
% 
% %   color encoding scheme
% 
% %   adapted from the color circle idea described at
% %   http://members.shaw.ca/quadibloc/other/colint.htm
% 
% 
% RY = 15;
% YG = 6;
% GC = 4;
% CB = 11;
% BM = 13;
% MR = 6;
% 
% ncols = RY + YG + GC + CB + BM + MR;
% 
% colorwheel = zeros(ncols, 3); % r g b
% 
% col = 0;
% %RY
% colorwheel(1:RY, 1) = 255;
% colorwheel(1:RY, 2) = floor(255*(0:RY-1)/RY)';
% col = col+RY;
% 
% %YG
% colorwheel(col+(1:YG), 1) = 255 - floor(255*(0:YG-1)/YG)';
% colorwheel(col+(1:YG), 2) = 255;
% col = col+YG;
% 
% %GC
% colorwheel(col+(1:GC), 2) = 255;
% colorwheel(col+(1:GC), 3) = floor(255*(0:GC-1)/GC)';
% col = col+GC;
% 
% %CB
% colorwheel(col+(1:CB), 2) = 255 - floor(255*(0:CB-1)/CB)';
% colorwheel(col+(1:CB), 3) = 255;
% col = col+CB;
% 
% %BM
% colorwheel(col+(1:BM), 3) = 255;
% colorwheel(col+(1:BM), 1) = floor(255*(0:BM-1)/BM)';
% col = col+BM;
% 
% %MR
% colorwheel(col+(1:MR), 3) = 255 - floor(255*(0:MR-1)/MR)';
% colorwheel(col+(1:MR), 1) = 255;
% 
% end
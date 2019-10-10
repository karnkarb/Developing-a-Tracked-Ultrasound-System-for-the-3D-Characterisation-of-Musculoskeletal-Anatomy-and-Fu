function [ Range_Of_Motion ] = RangeOfMotion ( Angular_Variance )
 
%-----------------------------------------------------------------------
%
%     Written by:
%     R. Ehrig
%     Zuse Institute Berlin (ZIB)
%     Takustrasse 7
%     D-14195 Berlin-Dahlem
%     phone : +49-30-84185-0
%     e-mail: ehrig@zib.de
%     URL   : http://www.zib.de
%
%-----------------------------------------------------------------------
%
% Set initial data
%
%-----------------------------------------------------------------------

Range_Of_Motion = zeros ( 3, 1 );

ROMfun = @(x,r) sin ( 0.5 * x ) - 0.5 * r * x;

for i = 1 : 3;
 
  if ( Angular_Variance(i) < 0.000001 )

     Range_Of_Motion(i) = 0.0;

  else

     r = 1.0 - Angular_Variance(i);

     Range_Of_Motion(i)  = ...
         fzero ( @(x) ROMfun ( x, r ), [ 0.000001, 2.0 * pi ] ) ...
         * 180.0 / pi;

   end;
   
end;

%-----------------------------------------------------------------------

return;
end

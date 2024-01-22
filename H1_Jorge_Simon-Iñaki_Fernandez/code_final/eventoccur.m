function [eventpos, ister, direct] = eventoccur(t,y)
 % Compute Mach number
 Mach = y(1) / sqrt((5 * y(3)) / 3);
 
 % Check when Mach number is 1
 eventpos = Mach - 1;
 %  Function cases
 ister = 1;
 direct = 0;
end
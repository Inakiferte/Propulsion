function [eventpos, ister, direct] = eventoccur(t,y)
 % Mars Case============================
 SOI_M = 5.84e8;                        % in m
 eventpos = y(1) - SOI_M;
 
 % Earth Case===========================
 %SOI_E = 9.27e8;                         % in m
 %eventpos = y(1) - SOI_E;

 % Both Cases
 ister = 1;
 direct = 0;
end
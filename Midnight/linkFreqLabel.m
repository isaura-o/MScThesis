function [C1, C1sor] = linkFreqLabel(tCC,tPC,tM3)
%
% Gets the frequencies of the 5 strongest links on from the table files.
%
tCC = textCCD.t1(:,2:3,:); tCC2 = string(tCC); numCC = str2double(regexp(tCC2,'\d*','Match'));
tmp = permute(numCC, [1 3 2]); concCC = reshape(tmp, [],2);
         
tM3 = textM3D.t1(:,2:3,:); tM32 = string(tM3); numM3 = str2double(regexp(tM32,'\d*','Match')); 
tmp = permute(numM3, [1 3 2]); concM3 = reshape(tmp, [],2);
 
tPC = textPCD.t1(:,2:3,:); tPC2 = string(tPC); numPC = str2double(regexp(tPC2,'\d*','Match'));
tmp = permute(numPC, [1 3 2]); concPC = reshape(tmp, [],2);


addc = [concCC; concM3; concPC];
[c, ia, ic] = unique(addc, 'rows');
freq=accumarray(ic,1)
C1 = [c, freq]; C1sor = sortrows(C1,3);


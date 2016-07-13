function []=micomsteps(baclin,batrop)
% []=MICOMSTEPS(BACLIN,BATROP) displays valid time steps for the model MICOM.
%    BACLIN and BATROP are two element vectors containing the lower and upper
%    baroclinic and barotropic time steps, respectively.

sd=86400;
bc=baclin(1):baclin(2);
bc=bc(find(rem(sd,bc)==0));

for n=1:length(bc)
  bt=batrop(1):batrop(2);
  bt=bt(find(rem(bc(n)./bt,2)==0));
  if ~isempty(bt)
    disp(' ')
    disp('baclin  batrop  fac')
    disp(sprintf('%6d  %6d   %d',bc(n),bt(1),bc(n)/bt(1)))
    for j=2:length(bt)
      disp(sprintf('        %6d   %d',bt(j),bc(n)/bt(j)))
    end
  end
end

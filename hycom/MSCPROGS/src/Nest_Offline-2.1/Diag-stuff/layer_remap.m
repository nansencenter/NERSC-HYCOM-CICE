function newint=layer_remap(oldint,olddens,newdens,newdp0);
% Fit one staircase profile to another -- Old staircase profile 
% is oldint/olddens. New staircase densities given by newdens,
% new staircase interfaces given by newint. Should handle 
% new staircase which has densities not in set (min(olddens), max(olddens))
%
% 

tol=1e-3;
oldkdm=prod(size(olddens));
newkdm=prod(size(newdens));
klist=zeros(newkdm,2);
newint=zeros(newkdm,1);
for k=1:newkdm

   % -- Go through old layers, find densities which are lower than 
   % -- target density (1:klist)
   for kold=1:oldkdm
      if (olddens(kold)<=newdens(k))
         klist(k,1)=kold;
      end 
   end 

   % -- Go through old layers, find densities which are higher than 
   % -- target density.
   for kold=oldkdm:-1:1
      if (olddens(kold)>=newdens(k))
         klist(k,2)=kold;
      end 
   end 
end 

figure(1) ; clf;
plot(klist(:,1)); hold on
plot(klist(:,2),'Color','r')

figure(2) ; clf;
I=find(klist(:,1)>0);
plot(I,olddens(klist(I,1))); hold on
I=find(klist(:,2)>0);
plot(I,olddens(klist(I,2)),'Color','r');
plot(newdens(:),'Color','g');

% Cycle new layers
kolast=1;
for k=1:newkdm

   % Last new layer upper interface
   if (k==1) 
      newup=0.;
   else
      newup=newint(k-1);
   end 
   disp([ 'newup:' num2str(newup) ])


   if (klist(k,1)==0) 
      % klist is zero - could not find layers lighter 
      % than the target density -> layer outcrops at surface, dp
      % should be zero in an isopycnic model, here we must consider
      % layer thickness restrictions
      newint(k)=newup+newdp0(k);

   elseif (klist(k,1)<kolast) 
      % There are water masses lighter than this one, but they are
      % above the current interface. We can not reach the
      % target densty. The next iteration may succeed, however, so
      newint(k)=newup+newdp0(k);


   elseif (klist(k,2)==0) 
      % klist2 is zero - no water masses are heavier than the target 
      % density. There is no hope of reaching our target density.... This
      % layer fills  the water column
      newint(k)=oldint(oldkdm);

   % There are elements heavier and lighter than the target density 
   % in the currently "available" old water column - further processing needed...
   else


      tmpdens=0.;
      initdens=0;
      sumdp=0.;
      match=0;
      for kold=1:oldkdm

         % Last old layer upper interface
         if (kold==1) 
            oldup=0.;
         else
            oldup=oldint(kold-1);
         end 
         disp([ '   oldup/dnw:' num2str(oldup) ' ' num2str(oldint(kold)) ])

         
         % fraction of this old layer available for mixing
         oldpfrac=max(0,oldint(kold) - max(newup,oldup));
         disp([ '   olddpfrac:' num2str(oldpfrac)])

         if (sumdp==0.) 
            sumdp=oldpfrac;
            tmpdens=olddens(kold);
            initdens=tmpdens;
         else
            % mixing a "newdp" layer of density olddens(kold) will result
            % in target density

            % Various checks
            if (abs(newdens(k)-tmpdens)<tol) 
               disp('   match');
               newdp=0.;
               %if (match==0)
               %   kolast=kold;
               %end
               match=1;
            elseif (newdens(k)<tmpdens) 
               newdp=oldpfrac;
            elseif (newdens(k)>olddens(kold)) 
               newdp=oldpfrac;
            else
               newdp=sumdp*(newdens(k)-tmpdens)/(olddens(kold)-newdens(k));
               newdp=max(0.,min(newdp,oldint(kold)-oldup));
            end 

            tmpdens=(olddens(kold)*newdp + tmpdens*sumdp)/(sumdp+newdp);
            sumdp=sumdp+newdp;
         end


      end

      newint(k)=newup+sumdp;
      if (match==0 & sumdp>1e-4) % No match -- use dp0
         newint(k)=newint(k)-sumdp;
         newint(k)=newint(k)+newdp0(k);
      end



   end

   if (k==newkdm) 
      newint(k)=oldint(oldkdm);
   end 
end 
figure(3) ; clf; 
stairs(olddens,-oldint) ; hold on
stairs(newdens,-newint,'Color','r') ; hold on

function RelDif = ModelChange(M1,M0)
% Check the changes in the transition probabilities.  Looks at the
% change in N.C(0,i,j) and the change in M.C(~0,i,j)
% relative to the maximum nonzero element.  The maximum of all these is
% returned.
epsi=.001;
[nu ns ns1]=size(M1.C);
A0=squeeze(M0.C(1,:,:));
A1=squeeze(M1.C(1,:,:));
Ad=max(max((abs(2*(A1-A0)/(A1+A0+epsi)))));
% Now we zero out the no-step entry and compute the rest.
M0.C(1,:,:)=0;
M1.C(1,:,:)=0;
N0=squeeze(sum(M0.C));
N1=squeeze(sum(M1.C));
Nd=max(max( squeeze(max(2*abs(M1.C-M0.C)))./(N1+N0+epsi) ));
% Ad
% Nd
RelDif=max(Ad,Nd);

%Shanmei Liu
% April 18, 2016

% Not specify the derivatives analytically, relying on built in function df_deriv.


% Additional Folder for Path
% The utility functions are stored in a separate folder, which has to be
% loaded in addition to |ddebiftool|. We also include the folders for
% periodic orbits and for normal form functions.
clear;
close all
%addpath('d:/My Documents/MATLAB/dde_biftool_v3.1/dde_biftool_v3.1/ddebiftool/',...
%    'd:/My Documents/MATLAB/dde_biftool_v3.1/dde_biftool_v3.1/ddebiftool_extra_psol/',...
%    'd:/My Documents/MATLAB/dde_biftool_v3.1/dde_biftool_v3.1/ddebiftool_utilities/',...
%    'd:/My Documents/MATLAB/dde_biftool_v3.1/dde_biftool_v3.1/ddebiftool_extra_nmfm/');
addpath('H:/My Documents/Desktop/dde_biftool_v3.1/dde_biftool_v3.1/ddebiftool/',...
    'H:/My Documents/Desktop/dde_biftool_v3.1/dde_biftool_v3.1/ddebiftool_extra_psol/',...
    'H:/My Documents/Desktop/dde_biftool_v3.1/dde_biftool_v3.1/ddebiftool_utilities/',...
    'H:/My Documents/Desktop/dde_biftool_v3.1/dde_biftool_v3.1/ddebiftool_extra_nmfm/');


% Definition of right-hand side
% Right hand side , xx is 1 (variables) x 2 (1 delay + 1 )
% Parameter vector has the order $[\delta, p, \alpha,\gamma, \tau]
Nicholson_sys_rhs=@(xx,par)[...
    exp(-par(1).*par(5))*par(2).*xx(1,2,:).*exp(-par(3).*xx(1,2,:))-par(4)*xx(1,1,:)];

funcs=set_funcs(...
    'sys_rhs',Nicholson_sys_rhs,...
    'sys_tau',@()[5],...
    'x_vectorized',true);
ind_tau=5;  % the index of delay tau --used later for continuation


% general continuation parameters, kept in a cell list
%original
parbd={'min_bound',[ind_tau, 0],'max_bound',[ind_tau,80],...
   'max_step',[ind_tau,0.2]};

% Steady state branches
%The convenience function |SetupStst| can be used to define the initial
% piece of a steady-state branch. Its first arguments is |funcs|, the
% others are name-value pairs. Important parameters:
% 
% * |'parameter'|: row vector of initial parameters
% * |'x'|: column vector of initial equilibrium
% * |'contpar'|: index of continuation parameter (or vector of indices)
% * |'step'|: initial step along branch (default 0.01)
%
% All other name-value pairs are appended as fields to the structures in
% branch1 if their names match. Of course, the |branch| structure can also
% be manipulated manually afterwards. The subsequent continuation extends
% the branch in both directions up to the boundaries |min_bound| and
% |max_bound|.

%
[branch1,suc]=SetupStst(funcs,...
  'parameter',[0.01, 15, 1,1, 65.428],'x',[1.5],...
     'contpar',ind_tau,'step',0.05,parbd{:});
 
branch1.method.continuation.plot=0; % don't plot prgress
[branch1,s1,f1,r1]=br_contn(funcs,branch1,100);
branch1=br_rvers(branch1);
[branch1,s1,f1,r1]=br_contn(funcs,branch1,100)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% continue in the second branch:
[branch2,suc]=SetupStst(funcs,...
  'parameter',[0.01, 15, 1,1, 6.2],'x',[1.5],...
     'contpar',ind_tau,'step',0.05,parbd{:});
 
branch2.method.continuation.plot=0; % don't plot prgress
[branch2,s2,f2,r2]=br_contn(funcs,branch2,100);
branch2=br_rvers(branch2);
[branch2,s2,f2,r2]=br_contn(funcs,branch2,100)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% continue in the third branch:
[branch3,suc]=SetupStst(funcs,...
  'parameter',[0.01, 15, 1,1, 0.214],'x',[1.5],...
     'contpar',ind_tau,'step',0.005,parbd{:});
 
branch3.method.continuation.plot=0; % don't plot prgress
[branch3,s3,f3,r3]=br_contn(funcs,branch3,100);
branch3=br_rvers(branch3);
[branch3,s3,f3,r3]=br_contn(funcs,branch3,100)



% Stability along branch
branch1.method.stability.minimal_real_part=-2;
%first output nunst, the number of unstable eigenvalues for bifurcation detection.
%nunst (vector of integers): number of unstable eigenvalues
nunst1=GetStability(branch1,'funcs',funcs); 
%point at which the number of unstable eigenvalues changes by 2
indhopf1=find(abs(diff(nunst1))==2)
if length(indhopf1)>1
    indhopf11=indhopf1(1)
else
    indhopf11=indhopf1
end

branch2.method.stability.minimal_real_part=-2;
%first output nunst, the number of unstable eigenvalues for bifurcation detection.
nunst2=GetStability(branch2,'funcs',funcs); 
%point at which the number of unstable eigenvalues changes by 2
indhopf2=find(abs(diff(nunst2))==2)
if length(indhopf2)>1
    indhopf21=indhopf2(1)
else
    indhopf21=indhopf2
end

branch3.method.stability.minimal_real_part=-2;
%first output nunst, the number of unstable eigenvalues for bifurcation detection.
nunst3=GetStability(branch3,'funcs',funcs); 
%point at which the number of unstable eigenvalues changes by 2
indhopf3=find(abs(diff(nunst3))==2)
if length(indhopf3)>1
    indhopf31=indhopf3(1)
else
    indhopf31=indhopf3
end


%Periodic orbit continuation
[branch4,suc]=SetupPsol(funcs,branch1,indhopf11,'degree',4,'intervals',40);
branch4.parameter.max_step(1,:)=[ind_tau,0.5];
[branch5,suc]=SetupPsol(funcs,branch2,indhopf21,'degree',4,'intervals',40);
branch5.parameter.max_step(1,:)=[ind_tau,0.5];
[branch6,suc]=SetupPsol(funcs,branch3,indhopf31,'degree',4,'intervals',40);
branch6.parameter.max_step(1,:)=[ind_tau,0.5];
figure(2); clf;
[branch4,s,f,r]=br_contn(funcs,branch4,250);
hold on;
[branch5,s,f,r]=br_contn(funcs,branch5,250);
[branch6,s,f,r]=br_contn(funcs,branch6,250);
xlabel('\tau');
ylabel('Amplitude');
hold off;

% Floquet multipliers

% Compute stability information
st_branch4=br_stabl(funcs,branch4,0,1);
[bifstab,crit]=GetStability(st_branch4,...
        'exclude_trivial',true,'critical',true);


floq_tau=arrayfun(@(x)x.parameter(ind_tau),st_branch4.point);
floq_re=real(crit);
floq=exp(floq_re);
figure(3); clf;
plot(floq_tau,floq)


% look at the period along the branch:
tau_per4=arrayfun(@(x)x.parameter(ind_tau),branch4.point);
periods4=[branch4.point.period];
figure(4); clf;
plot(tau_per4,periods4,'b.-');
hold on;
tau_per5=arrayfun(@(x)x.parameter(ind_tau),branch5.point);
periods5=[branch5.point.period];
plot(tau_per5,periods5,'b.-');
tau_per6=arrayfun(@(x)x.parameter(ind_tau),branch6.point);
periods6=[branch6.point.period];
plot(tau_per6,periods6,'b.-');
hold off;
xlabel('\tau');
ylabel('period');

% Periodic solution
figure(5);clf;
%psol=branch4.point(end)
psol=branch5.point(10)
p_pplot(psol);
xlabel('time/period');ylabel('x');


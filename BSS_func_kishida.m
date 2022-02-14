function BSS_func_kishida(filename,double_b,smp)

char_sign='_kk';
x=double_b;
[n,t] = size(x);
x=x-(ones(t,1)*mean(x'))';
double_b=x;

%mica=1;
%mbss=2;

%%%a_make_M_hs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
remove_b=double_b;
[bbb, bline]=size(remove_b);
remove_b=remove_b-(ones(bline,1)*mean(remove_b'))';
M = cov(remove_b');
S = inv(sqrtm(M));
xs=S*remove_b;
clear remove_b
average_xs=(ones(bline,1)*mean(xs'))';
xs=xs-average_xs;
clear average_xs
 
[n,t] = size(xs);


rf=input('  rf = ');
k=input('  k = ');
for tau=1:k;
    wait_str=strcat('Calculating time delayed correlation-matrices.@@',num2str((tau/k)*100),'%');
    h=waitbar(0,wait_str);
    waitbar(tau/k);
    close(h)
    Mk(:,:,tau)=corrm(xs',(smp/rf)/tau);
end

M = reshape(Mk(:,:,:),n,n*k);

% joint diagonalization
        C=joint_diag_real(M,0.001); % joint_diag gives back C^-1

clear xs

% %%%a_make_B_z0_hs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=C'*S; %joint_diag_real
z0=B*double_b; %B=inv(A) x=As => s=Bx
A=inv(B);


 filename6=strcat(filename,'_B_k',int2str(k),char_sign)
 eval(['save ',filename6,' A B smp rf k filename'])
 filename7=strcat(filename,'_z0_k',int2str(k),char_sign)
 eval(['save ',filename7,' z0 smp rf k filename'])



function BSS_func_kk(filename,double_b,smp)

char_sign='_kk';
x=double_b;
mmean0=menu('Do you take the zero mean average','yes','no')
if(mmean0==1)
[n,t] = size(x);
x=x-(ones(t,1)*mean(x'))';
end
double_b=x;

mica=menu('ICA menu','realî≈','complexî≈')
mbss=menu('BSS menu','kTå^','T/kå^','constå^','Tangó¨','T/kå^Å{kTå^','specialå^','double T/k-type of eplipsy','double T/k-type of SEF','double kT-type of SEF')

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
        wait_str=strcat('Calculating time delayed correlation-matrices.Å@Å@',num2str((tau/k)*100),'%');
        h=waitbar(0,wait_str);
        waitbar(tau/k);
        close(h)
        if (mbss==1)
            Mk(:,:,tau)=corrm(xs',tau*(smp/rf));
        elseif(mbss==2)
            Mk(:,:,tau)=corrm(xs',(smp/rf)/tau);
        elseif(mbss==3)
            Mk(:,:,tau)=corrm(xs',(smp/rf)+tau-1);
        elseif(mbss==5)
            Mk(:,:,tau)=corrm(xs',(smp/rf)/tau);
        elseif(mbss==6)
            Mk(:,:,tau)=corrm(xs',tau);
        elseif(mbss==7)
            Mk(:,:,tau)=corrm(xs',(smp/rf)/tau);
            if(rf==10)
                Mk(:,:,tau+k)=corrm(xs',(smp/(rf/10))/tau);
            elseif(rf==20)
                Mk(:,:,tau+k)=corrm(xs',(smp/(rf/10))/tau);
            end
        elseif(mbss==8)
            Mk(:,:,tau)=corrm(xs',(smp/rf)/tau);
        elseif(mbss==9)
            Mk(:,:,tau)=corrm(xs',tau*(smp/rf));
        end
    end
end



% joint diagonalization
    bssmm=menu('value of jthresh','0.001','0.0000001','special','0.01')
    if(bssmm==1)
        C=joint_diag_real(M,0.001); % joint_diag gives back C^-1
    elseif(bssmm==2)
        C=joint_diag_real(M,0.0000001); % joint_diag gives back C^-1
    elseif(bssmm==3)
        sess=input('please input sess value : ');
        [C D count]=joint_diag_real(M,sess);
        disp('disp(count)=')
        disp(count)
     elseif(bssmm==4)
        C=joint_diag_real(M,0.01); % joint_diag gives back C^-1
    end
 %   B=C'*S;

clear xs

% %%%a_make_B_z0_hs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[n,p] = size(double_b);
%[brank bline]=size(double_b);
%double_b=double_b-(ones(bline,1)*mean(double_b'))';
% B=real(C'*S); %joint_diag
B=C'*S; %joint_diag_real
z0=B*double_b; %B=inv(A) x=As => s=Bx
A=inv(B);


à»â∫ÇÕÉfÅ[É^ï€ë∂ÇÃÇΩÇﬂ
    filename6=strcat(filename,'_B_k',int2str(k),char_sign)
    eval(['save ',filename6,' A B smp rf k filename'])
    filename7=strcat(filename,'_z0_k',int2str(k),char_sign)
    eval(['save ',filename7,' z0 smp rf k filename'])



function [V, D, count] =  joint_diag_real(A,jthresh)

%[Input]
%  A       : m*mn matrix [A1 A2 ... An]. each Ai is m*m matrix.
%  jthresh : threshold an optional small number (typically = 1.0e-8).
%
%[Output]
%  V       : an m*m unitary matrix, which diagnalize each matrix A as below.
%  D       : V'*A1*V , ... , V'*An*V, which was approximately diagnalized.
%  count   : loop count
%

[m,nm] = size(A); %% n mxm matrix

V = eye(m);  %%
count=1;
encore = 1; 
while encore, 
	encore=0;
	for p=1:m-1, %% 
		Ip = p:m:nm; %% 1xn
		for q=p+1:m, 
			Iq = q:m:nm; %% 1xn
			%% Computing the Givens angles
			g       = [ A(p,Ip)-A(q,Iq)  ; A(p,Iq)+A(q,Ip) ]; %% h(A)^T for m matrix (2xn)
			[vcp,D] = eig(g*g');
			[la, K] = sort(diag(D));
			angles  = vcp(:,K(2));
			if angles(1)<0 , angles= -angles ; end ;
			c       = sqrt(0.5+angles(1)/2);
			s       = 0.5*angles(2)/c; 
			if abs(s)>jthresh, %%% updates matrices A and V by a Givens rotation
				encore          = 1 ;
				pair            = [p;q] ;
				G               = [ c -s ; s c ] ;
				V(:,pair)       = V(:,pair)*G ;
				A(pair,:)       = G' * A(pair,:) ;
				A(:,[Ip Iq])    = [ c*A(:,Ip)+s*A(:,Iq) -s*A(:,Ip)+c*A(:,Iq) ]; 
			end
        end%% q loop
    end%% p loop
    count=count+1;
    if count>2000
        break
    end
end%% while
D = A ;
return


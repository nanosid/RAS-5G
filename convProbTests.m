%% Initialization
clear;
M = 16;
K = 8;
gamma = 10;
e = 10^(-220/10);

%% Power minimization using SOCP
H = randn(K,M,'like',1i);
cvx_begin
	variable tau nonnegative
	variable W(M,K) complex
	minimize(tau)
	subject to
		for i = 1:K
			norm([H(i,:)*W 1]) <= sqrt(1+(1/gamma))*real(H(i,:)*W(:,i))
			imag(H(i,:)*W(:,i)) == 0
		end
		norm(vec(W)) <= sqrt(tau)
cvx_end

%% Power minimization using SDP
H = randn(M,M,K,'like',1i);
cvx_begin
	variable X(M,M,K) hermitian
	variable s(K,1) % slack variable for conformity of SDP with CVX
	obj = 0;
	for i = 1:K
		obj = obj + trace(X(:,:,i));
	end
	minimize(obj)
	subject to
		for i = 1:K
			const = 0;
			for j = 1:K
				if j ~= i
					const = const + trace(H(:,:,i)*X(:,:,j));
				end
			end
			trace(H(:,:,i)*X(:,:,i)) - gamma*const - s(i,1) == gamma % in CVX, SDP constraints need to be equalities
			s(i,1) >= 0
			X(:,:,i) == hermitian_semidefinite(M)
		end
cvx_end

%% Multicast power minimization SDP
H = randn(M,M,K,'like',1i);
cvx_begin
	variable X(M,M) hermitian
	variable s(K,1)
	minimize(trace(X))
	subject to
		for i = 1:K
			trace(X*H(:,:,i)) - s(i,1) == 1
			s(i,1) >= 0
		end
		X == hermitian_semidefinite(M)
cvx_end

%% Multicast beamforming max-min SNR with RAS int const
h_k = randn(M,K,'like',1i);
H = zeros(M,M,K);
for i = 1:K
	H(:,:,i) = h_k(:,i)*h_k(:,i)';
end
hgR = randn(M,1,'like',1i);
HgR = hgR*hgR';
cvx_begin
	variable tau nonnegative
	variable X(M,M) hermitian
	variable s(K,1)
	%variable t(1,1)
	minimize(-tau)
	subject to
		norm(diag(X*HgR)) <= e
		for i = 1:K
			trace(X*H(:,:,i)) - s(i) == tau
			s(i) >= 0
		end
		trace(X) <= 1
		%t <= 0
		X == hermitian_semidefinite(M)
cvx_end

%% Gaussian Randomization for approx rank-1 solution
L = 10000;
x = zeros(M,L);
% obj = zeros(1,L);
for i = 1:L
	% x(:,i) = randn(M,1,'like',1i);
	% x(:,i) = sqrt(2)*(chol(X)*[real(x(:,i)) imag(x(:,i))])*[1; 1i];
	x(:,i) = randn(M,1,'like',1i);
	[U,D,~] = eig(X);
	x(:,i) = U*sqrt(D)*x(:,i);
	x(:,i) = x(:,i)./norm(x(:,i));
	scale = e./(norm(x(:,i)'*hgR).^2);
	x(:,i) = x(:,i).*(scale.^(0.5));
	% min_idx = 1;
	% for j = 2:K
	% 	if trace(x(:,i)*x(:,i)'*H(:,:,j)) < trace(x(:,i)*x(:,i)'*H(:,:,min_idx))
	% 		min_idx = j;
	% 	end
	% end
	% scale = tau/trace(x(:,i)*x(:,i)'*H(:,:,min_idx));
	% x(:,i) = x(:,i)*scale;
end
[~,idx] = max(min(abs(x'*h_k),[],2));
% temp = zeros(1,L);
% for i = 1:L
% 	temp(i) = 10*log10(norm(x(:,i)'*hgR).^2);
% end
% min(temp)
% [~,idx] = min(temp);
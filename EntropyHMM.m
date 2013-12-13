clear all
N=3; 
% threshold = 0;
% q = max(rand(N) - threshold, 0); 
% q = q./repmat(sum(q, 2), 1, N);
% q = [2/3 1/6 1/6; 2/3 1/6 1/6; 2/3 1/6 1/6];
% q = [2/3 1/3 0; 2/3 1/3 0; 2/3 1/3 0];
% q = [1 0 0; 1 0 0; 1 0 0];
% q = ones(N)/N;
% p = ones(N)/N;

% q = [2/3 1/6]; % 2 free parameters
% q = [2/3 1/3]; % 2 free parameters
% q = [1/2 1/2]; % 1 free parameter
% q = [1 0]; % 1 free parameter
% q = ones(1, N - 1) / N; % 1 free parameter
p = ones(1, N - 1) / N;

k_entropy = q.*log2(q)
k_rel_entropy = q.*log2(q./p)
k_entropy_total = k_entropy .* k_rel_entropy
k_entropy_total(find(isnan(k_entropy_total))) = 0
k_entropy_sum = exp(-sum(k_entropy_total(:)))*size(q, 1)

% k_cross_entropy = q.*log2(p) + q.*log2(q)
% k_cross_entropy(find(isnan(k_cross_entropy))) = 0
% k_cross_entropy = -sum(k_cross_entropy(:))

% k_nnz = nnz(q) - size(q, 2) + sum(q(:) == 1);

% k_test = (min(1 - max(q'), max(q') - 1/N) + 1/N)*N
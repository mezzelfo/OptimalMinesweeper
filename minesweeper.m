clear all
%clc

A = [[1 0 0 0 0 0 0 0 0 0 0 0]
 [1 0 0 0 0 0 0 0 0 0 0 0]
 [1 1 0 0 0 0 0 0 0 0 0 0]
 [1 1 1 1 1 0 0 0 0 0 0 0]
 [0 0 0 1 1 1 0 0 0 0 0 0]
 [0 0 0 0 0 0 1 1 0 0 0 0]
 [0 0 0 0 0 0 1 1 1 0 0 0]
 [0 0 0 0 0 0 0 1 1 1 0 0]
 [0 0 0 0 0 0 0 0 1 1 1 0]
 [0 0 0 0 1 1 0 0 0 1 1 1]] ;
b = [1 1 1 2 1 1 1 2 1 1] ;

[infocount, varscount] = size(A);
size(b)

X = optimvar('X',varscount,'LowerBound',0,'UpperBound',1);
prob = optimproblem;
%prob.Constraints.cons1 = sum(X) <= 10;
prob.Constraints.cons1 = A*X == b';
prob.Objective = sum(X.*(X-1),'all');
prob.ObjectiveSense = 'min';
sol = solve(prob);
sol.X

%lsqlin(A,b,eye(25),ones(1,25))


%%
H = 8;
W = 8;
N = 10;
true_bombs = zeros(H,W);
while sum(true_bombs,'all') < N
    r = randi([1,H]);
    c = randi([1,W]);
    if true_bombs(r,c) == 0
        true_bombs(r,c) = 1;
    end
end
ker = [1 1 1;1 1 1;1 1 1];

%all_infos = ifft2(fft2(true_bombs,H+2,W+2).*fft2(ker,H+2,W+2));
%all_infos = all_infos(2:1+H,2:1+W);
%norm(all_infos-conv2(true_bombs,ker,'same'))
A = toeplitz([1 1 0 0 0 0 0 0]);
P = kron(A,A);
all_infos = reshape(P*true_bombs(:),[H,W]);
isequal(all_infos,conv2(true_bombs,ker,'same'))

%[tuValue, indices] = tu(P)
disp("P is not totally unimodular:")
disp("    det(P([1,3,18],[2,9,11])) = "+det(P([1,3,18],[2,9,11])))


%figure(1)
%imagesc(reshape(pinv(P)*all_infos(:),[H,W]))
%figure(2)
%imagesc(true_bombs)
%%
F = diag(randi([0,100],1,H*W)>80);

bombs = optimvar('bombs',[H,W],'Type','integer','LowerBound',0,'UpperBound',1);
prob = optimproblem;
prob.Constraints.cons1 = sum(bombs,'all') == N;
prob.Constraints.cons1 = F*P*bombs(:) == F*all_infos(:);
tot = zeros(H,W);
for i=1:H*W
    prob.Objective = bombs(i);
    prob.ObjectiveSense = 'max';
    sol = solve(prob,'options',optimoptions('intlinprog','Display','none'));
    tot = tot + sol.bombs;
end
figure(1)
imagesc(tot > 0)
figure(2)
imagesc(true_bombs)
tot >= true_bombs

%%
clc

F = diag(randi([0,100],1,H*W)>80);

bombs = optimvar('bombs',[H,W],'LowerBound',0,'UpperBound',1);
prob = optimproblem;
prob.Constraints.cons1 = sum(bombs,'all') == N;
prob.Constraints.cons1 = F*P*bombs(:) == F*all_infos(:);
prob.Objective = sum(bombs.*(bombs-1),'all');
prob.ObjectiveSense = 'min';
sol = solve(prob);
figure(1)
imagesc(sol.bombs)
figure(2)
imagesc(true_bombs)


%%
% bombs = optimvar('bombs',[H,W],'Type','integer','LowerBound',0,'UpperBound',1);
% prob = optimproblem('Objective',sum(bombs,'all'),'ObjectiveSense','max');
% prob.Constraints.cons1 = F*P*bombs(:) >= F*all_infos(:);
% sol = solve(prob);
% sum(sol.bombs,'all')
% sol.bombs >= true_bombs

% figure(1)
% imagesc(sol.bombs)
% figure(2)
% imagesc(true_bombs)
%figure(3)
%imagesc(reshape(diag(F),[H,W]))

%yolo = ifft2(fft2(all_infos,H+2,W+2)./fft2(ker,H+2,W+2));
%yolo = yolo(1:H,1:W);
%norm(yolo-true_bombs)
%figure(1)
%imagesc(yolo)
%figure(2)
%imagesc(true_bombs)
% figure(3)
% imagesc(all_infos)
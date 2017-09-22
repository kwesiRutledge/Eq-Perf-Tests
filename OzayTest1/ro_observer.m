%function to design observer

function L = ro_observer(Auu, Aau, Eu, Ea, M, W, p)
% p is the norm type

n = size(Auu,1);
m = size(Aau,1);

Lvar = sdpvar(n,m);

vertW = W.V;
const = [];
for i=1:size(vertW,1)
    v = vertW(i,:)';
    const = [const, norm(Auu-Lvar*Aau,p)+(norm((Eu-Lvar*Ea)*v,p)/M) <=1];
end

% objective is a bit arbitrary -- better objectives can be picked to
% improve synthesis performance
objective = norm(Auu-Lvar*Aau,p);

diag = optimize(const, objective);
if diag.problem == 0;
    L = value(Lvar);
    L(isnan(L))=0; 
else
    display('infeasible or problem with optimization');
    L = nan; %infeasible
end

%sanity check
for i=1:size(vertW,1)
    v = vertW(i,:)';
    norm(Auu-L*Aau,p)+(norm((Eu-L*Ea)*v,p)/M)
end
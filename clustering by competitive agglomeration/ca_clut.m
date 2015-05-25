function [ru,rv] = ca_clut(ej,no,cj,xj)
%% Variable Definition 
[len rr]= size(xj);
cmax = cj;
cf = 1;
ch = 1;
e_val = ej;
n_val = no;
t_val = 10;
breakval = 0.8;
iterval = 5;
%%
[center,U,obj_fcn] = fcm(xj,cmax,[2 iterval NaN NaN]);

%% calculate Ni
Ni = ones(1,cmax);
for j = 1:cmax
    Ni(1,j) = sum(U(j,:));
end
%%
V = ccenter(U,xj);

k = 0;
Nj = ones(1,len); 
%%
c_max=cmax;
while( 1 )
    [cmax r] = size(U);
    dis = zeros(cmax,len);
    for i = 1 : cmax
        for j = 1 : len
            dis(i,j) = norm(xj(j) - V(i)).^2;
        end
    end
    
U = U.^2;
%% calculate alpha
    temp = n_val * exp(-k/t_val);
    sum1 = 0;
    for i = 1 : cmax
        for j = 1:len
        sum1 = sum1 + U(i,j)*dis(i,j);
        end
    end
    sum2 = sum(Ni.^2);
    alpha = temp * (sum1 / sum2); % alpha change depent on time
 % alpha =1; % alpah is static
%%
    UF = zeros(cmax,len);
    UB = zeros(cmax,len);
    for i = 1 : cmax
        for j = 1 : len
            UF(i,j) = (dis(i,j).^-1)/ sum((dis(:,j).^-1));
        end
    end
    for j = 1 : len
         Nj(j) = dot(dis(:,j).^-1,Ni)/sum((dis(:,j).^-1));
    end
    for i = 1 : cmax
        for j = 1 : len
            UB(i,j) = alpha*(Ni(i) - Nj(j))/dis(i,j);
        end
    end
    %%
    U = UF + UB;
    negpos = U <= 0;
    U(negpos) = 0;
    lar1pos = U >= 1;
    U(lar1pos) = 1;
    
     
%% Calculate Ni again 
    for i = 1:cmax
        Ni(i) = sum(U(i,:));
    end
    
    pos = Ni < e_val;
    kn = sum(pos);
    flag = zeros(1,kn);
    ii = 1;
    for i = 1 : cmax
        if (pos(i) == 1)
            flag(ii) = i;
            ii = ii + 1;
        end           
    end
    %% Discard cluster
    U(flag,:) = [];
    Ni(flag) = [];
    [cf rc] = size(U);
    %%
    chgV = ccenter(U,xj);
    if ( size(V) == size(chgV))
        ch = sum(abs(V - chgV));
    end
    V = chgV;  
    if (cmax == cf && ch <=breakval) % break condition
        break;
    end
    k = k + 1;
end

sortV = sort(V);
disp('centers :');
disp(sortV);
%% sorting U and V matrix
% xu = zeros(1,cf);
% temp = zeros(1,cf);
% for i = 1:cf
%    temp = find( V == sortV(i));
%    xu(i) = temp(1,1);
% end
% sortU = zeros(cf,len);
% for i = 1:cf
%     sortU(i,:) = U(xu(i),:);
% end
% %%
% 
  ru = U;
  rv = sortV;
 



%% CCenter function  
function out = ccenter(u,data)
[c r] = size(u);
temp = ones(1,c);
for i = 1 : c
    temp(i) = dot((u(i,:).^2),data)/dot(u(i,:),u(i,:));
end
out = temp;
end 

end 

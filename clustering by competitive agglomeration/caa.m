function caa(d)


Data =d;
[r c k] = size(Data);

parameters = [7 3.5 18];%;5 5 20;3.4 5 20;4 1 20;3 4 18;4 7 18;4 5.5 18;6.2 7 22;5 5 18;4 5 26];
[property_size data_n] = size(parameters);
dx=Data(:,1)
dy=Data(:,2)  

%for i = 1 : property_size
   % call ca_clut(7,3.5,18 ,Data(:,1))
    [U,V] = ca_clut(parameters(1,1),parameters(1,2),parameters(1,3),dy);
     
    center_length = length(V);
    [uc ur] = size(U);
 
  
%end


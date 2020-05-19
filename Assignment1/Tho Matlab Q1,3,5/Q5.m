clear all
close all
load -ascii traffic.mat
load -ascii capacities.mat
load -ascii traveltime.mat
load -ascii flow.mat
W = zeros(13,1);
s = [];
t = [];
 for i = 1:13
    for j = 1:28
       if traffic(i,j) == 1
           s(1,length(s)+1) = i;
        for k = 1:17
           if traffic(k,j) == -1
               t(1,length(t)+1) = k;
               W(i,k) = 1;
           end
        end
       end
    end
end


rowswitch1 = eye(28);
rowswitch2(1,:) = rowswitch1(1,:);
rowswitch2(3,:) = rowswitch1(2,:);
rowswitch2(5,:) = rowswitch1(3,:);
rowswitch2(8,:) = rowswitch1(4,:);
rowswitch2(2,:) = rowswitch1(5,:);
rowswitch2(11,:) = rowswitch1(6,:);
rowswitch2(13,:) = rowswitch1(7,:);
rowswitch2(15,:) = rowswitch1(8,:);
rowswitch2(18,:) = rowswitch1(9,:);
rowswitch2(4,:) = rowswitch1(10,:);

rowswitch2(6,:) = rowswitch1(11,:);
rowswitch2(7,:) = rowswitch1(12,:);
rowswitch2(9,:) = rowswitch1(13,:);
rowswitch2(10,:) = rowswitch1(14,:);
rowswitch2(12,:) = rowswitch1(15,:);
rowswitch2(19,:) = rowswitch1(16,:);
rowswitch2(20,:) = rowswitch1(17,:);
rowswitch2(14,:) = rowswitch1(18,:);
rowswitch2(16,:) = rowswitch1(19,:);
rowswitch2(17,:) = rowswitch1(20,:);

rowswitch2(21,:) = rowswitch1(21,:);
rowswitch2(23,:) = rowswitch1(22,:);
rowswitch2(24,:) = rowswitch1(23,:);
rowswitch2(22,:) = rowswitch1(24,:);
rowswitch2(25,:) = rowswitch1(25,:);
rowswitch2(26,:) = rowswitch1(26,:);
rowswitch2(27,:) = rowswitch1(27,:);
rowswitch2(28,:) = rowswitch1(28,:);

capacities2 = rowswitch2*capacities;
flow2 = rowswitch2*flow;
traveltime2 = rowswitch2*traveltime;

G = digraph(s,t,capacities2);

graphshortestpath(sparse(W),1,17,'Weights',traveltime2);
graphmaxflow(sparse(W),1,17,'CAPACITY',capacities2);


outflow = (traffic>0)*flow;
inflow = (traffic<0)*flow;
net_flow = inflow-outflow;


%% 5.d
clear f
e=28;

in_out_flow = zeros(17,1);
in_out_flow(1) = outflow(1);
in_out_flow(end) = -outflow(1);


cvx_begin
    variable f(e)
    minimize( sum(sum((traveltime.*capacities).*inv_pos(ones(e,1)-f./capacities))-traveltime.*capacities));
    subject to
        traffic*f == in_out_flow
        0 <= f <= capacities
cvx_end 

f_social = f;

%% 5.e
clear f
clc
e=28;
in_out_flow = zeros(17,1);
in_out_flow(1) = outflow(1);
in_out_flow(end) = -outflow(1);

cvx_begin
    variable f(e)   
    minimize( sum(-traveltime.*capacities.*log(1-f./capacities)));
    subject to
        traffic*f == in_out_flow
        0 <= f <= capacities
cvx_end 

f_wardrop = f;

%% 5.f
clear f
clc
omega = f_social.*(traveltime./capacities)./((ones(e,1)-f_social./capacities).^2);

e=28;
in_out_flow = zeros(17,1);
in_out_flow(1) = outflow(1);
in_out_flow(end) = -outflow(1);

cvx_begin
    variable f(e)   
    minimize( sum(omega-traveltime.*capacities.*log(1-f./capacities)));
    subject to
        traffic*f == in_out_flow
        0 <= f <= capacities
cvx_end 

f_tolls = f;

%% 5.d
clear f
e=28;
in_out_flow = zeros(17,1);
in_out_flow(1) = outflow(1);
in_out_flow(end) = -outflow(1);

cvx_begin
    variable f(e)
    minimize( sum(sum((traveltime.*capacities).*inv_pos(ones(e,1)-f./capacities))-traveltime.*capacities-f.*traveltime));
    subject to
        traffic*f == in_out_flow
        0 <= f <= capacities
cvx_end 

f_additional = f;


%% compare

delay_social = sum(traveltime.*inv_pos(ones(e,1)-f_social./capacities))
delay_additional = sum(traveltime.*inv_pos(ones(e,1)-f_additional./capacities))

delay_wardrop = sum(traveltime.*inv_pos(ones(e,1)-f_wardrop./capacities))
delay_tolls = sum(traveltime.*inv_pos(ones(e,1)-f_tolls./capacities))







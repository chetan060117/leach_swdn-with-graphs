%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
%                     LEACH-SWDN Implementation                        %
%                                                                      %                                  
%     "A clustering algorithm based on energy information              %
%                      and cluster heads                               %
%            expectation for wireless sensor networks "                %                                                             
%                                                                      %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBMITTED BY-                                                        %
%                SE20UCSE071- K.Chetan                     %
%                     (B.Tech CSE-1)                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Field Dimensions - x and y maximum (in meters)
xm = 100;
ym = 100;

% x and y Coordinates of the Sink
sink.x = 0.5 * xm;
sink.y = 1.75 * ym;

% Number of Nodes in the field
n = 100;

% Optimal Election Probability of a node to become cluster head
p = 0.05;

% Energy Model (all values in Joules)
% Initial Energy 
Eo = 0.5;
% Eelec=Etx=Erx
ETX = 50 * 0.000000001;
ERX = 50 * 0.000000001;
% Transmit Amplifier types
Efs = 10 * 0.000000000001;
Emp = 0.0013 * 0.000000000001;
% Data Aggregation Energy
EDA = 5 * 0.000000001;

% Values for Hetereogeneity
% Percentage of nodes that are advanced
m = 0.05;
% \alpha
a = 0.1;

% Maximum number of rounds
rmax = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

% Computation of do
do = sqrt(Efs / Emp);

% Creation of the random Sensor Network
figure(1);

for i = 1:n
    S(i).xd = rand(1, 1) * xm;
    XR(i) = S(i).xd;
    S(i).yd = rand(1, 1) * ym;
    YR(i) = S(i).yd;
    S(i).G = 0;
    % initially there are no cluster heads only nodes
    S(i).type = 'N';
   
    temp_rnd0 = i;
    % Random Election of Normal Nodes
    if (temp_rnd0 >= m * n + 1) 
        S(i).E = Eo;
        S(i).ENERGY = 0;
        plot(S(i).xd, S(i).yd, 'o');
        hold on;
    end
    % Random Election of Advanced Nodes
    if (temp_rnd0 < m * n + 1)  
        S(i).E = Eo * (1 + a);
        S(i).ENERGY = 1;
        plot(S(i).xd, S(i).yd, '+');
        hold on;
    end
end

S(n+1).xd = sink.x;
S(n+1).yd = sink.y;
plot(S(n+1).xd, S(n+1).yd, 'x');

% Preallocate arrays
enrgy_res = zeros(n, rmax) + Eo;
energy_moy = zeros(n, rmax) + Eo;
Energy_disp = zeros(1, rmax);
PACKETS_TO_CH = zeros(1, rmax);
PACKETS_TO_BS = zeros(1, rmax);
DEAD = zeros(1, rmax);
DEAD_N = zeros(1, rmax);
DEAD_A = zeros(1, rmax);
ALIVE_NODE = zeros(1, rmax);
CLUSTERHS = zeros(1, rmax);
avg_r = zeros(1, rmax);
energy_ = zeros(1, rmax);
PDR = zeros(1, rmax);

% First Iteration
figure(1);

% Counter for CHs per round
rcountCHs = 0;
% Counter for CHs
countCHs = 0;
cluster = 1;
countCHs;
rcountCHs = rcountCHs + countCHs;
flag_first_dead = 0;

% Last round
last = rmax;

for r = 1:rmax
    fprintf('Round %d\n', r);

    % Operation for epoch
    if mod(r-1, round(1/p)) == 0
        for i = 1:n
            S(i).G = 0;
            S(i).cl = 0;
        end
    end

    hold off;

    % Number of dead nodes
    dead = 0;
    % Number of dead Advanced Nodes
    dead_a = 0;
    % Number of dead Normal Nodes
    dead_n = 0;
    % Alive nodes
    alive = 0;

    % Counter for bits transmitted to Base Station and to Cluster Heads
    packets_TO_BS = 0;
    packets_TO_CH = 0;
    % Counter for bits transmitted to Base Station and to Cluster Heads per round
    PACKETS_TO_CH(r) = 0;
    PACKETS_TO_BS(r) = 0;

    figure(1);

    for i = 1:n
        % Checking if there is a dead node
        if S(i).E <= 0
            plot(S(i).xd, S(i).yd, 'black .');
            dead = dead + 1;
            if S(i).ENERGY == 1
                dead_a = dead_a + 1;
            end
            if S(i).ENERGY == 0
                dead_n = dead_n + 1;
            end
            hold on;
        end
        if S(i).E > 0
            S(i).type = 'N';
            if S(i).ENERGY == 0  
                plot(S(i).xd, S(i).yd, 'o');
            end
            if S(i).ENERGY == 1  
                plot(S(i).xd, S(i).yd, '+');
            end
            hold on;
        end
    end
    plot(S(n+1).xd, S(n+1).yd, 'red x');

    DEAD(r) = dead;
    DEAD_N(r) = dead_n;
    DEAD_A(r) = dead_a;
    alive = n - dead;
    ALIVE_NODE(r) = alive;

    % Early Termination (check for 10% alive nodes)
    if alive < (n * 0.1)
        last = r;
        break;
    end

    % When the first node dies
    if dead == 1 && flag_first_dead == 0
        first_dead = r;
        flag_first_dead = 1;
    end

    countCHs = 0;
    cluster = 1;

    for i = 1:n
        if S(i).E > 0
            % Average energy
            energy_moy(i, r) = (1/r) * sum(enrgy_res(i, 1:r));
            % Generate a random number within the interval [0, (Eaverage_energy(i) / Eo)]
            temp_rand = (energy_moy(i, r) / Eo) * rand(1, 1);

            if S(i).G <= 0
                % k = (alive nodes in the network * %c(number of cluster heads defined at the initial time))
                k = ALIVE_NODE(r) * p;
                % Election of Cluster Heads
                if temp_rand <= (k * (S(i).E / Eo) / (n - k * mod(r, round(n/k))))
                    countCHs = countCHs + 1;
                    packets_TO_BS = packets_TO_BS + 1;
                    PACKETS_TO_BS(r) = packets_TO_BS;
                    
                    S(i).type = 'C';
                    S(i).G = 1/k;
                    C(cluster).xd = S(i).xd;
                    C(cluster).yd = S(i).yd;
                    plot(S(i).xd, S(i).yd, 'k*');
                    
                    distance = sqrt((S(i).xd - S(n+1).xd)^2 + (S(i).yd - S(n+1).yd)^2);
                    C(cluster).distance = distance;
                    C(cluster).id = i;
                    X(cluster) = S(i).xd;
                    Y(cluster) = S(i).yd;
                    cluster = cluster + 1;

                    % Calculation of Energy dissipated
                    % Size of data package is taken as 4000 units in bits
                    if distance > do
                        S(i).E = S(i).E - ((ETX + EDA) * 4000 + Emp * 4000 * (distance^4));
                        enrgy_res(i, r) = S(i).E;
                        Energy_disp(r) = Energy_disp(r) + ((ETX + EDA) * 4000 + Emp * 4000 * (distance^4));
                    else
                        S(i).E = S(i).E - ((ETX + EDA) * 4000 + Efs * 4000 * (distance^2));
                        enrgy_res(i, r) = S(i).E;
                        Energy_disp(r) = Energy_disp(r) + ((ETX + EDA) * 4000 + Efs * 4000 * (distance^2));
                    end
                end  
            end
        end 
    end

    CLUSTERHS(r) = countCHs;

    % Election of Associated Cluster Head for Normal Nodes
    for i = 1:n
        if S(i).type == 'N' && S(i).E > 0
            if cluster-1 >= 1
                min_dis = inf;
                min_dis_cluster = 0;
                for c = 1:cluster-1
                    temp_dis = sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2);
                    if temp_dis < min_dis
                        min_dis = temp_dis;
                        min_dis_cluster = c;
                    end
                end

                % Energy dissipated by associated Cluster Head
                if min_dis > do
                    S(i).E = S(i).E - (ETX * 4000 + Emp * 4000 * (min_dis^4));
                    enrgy_res(i, r) = S(i).E;
                    Energy_disp(r) = Energy_disp(r) + (ETX * 4000 + Emp * 4000 * (min_dis^4));
                else
                    S(i).E = S(i).E - (ETX * 4000 + Efs * 4000 * (min_dis^2));
                    enrgy_res(i, r) = S(i).E;
                    Energy_disp(r) = Energy_disp(r) + (ETX * 4000 + Efs * 4000 * (min_dis^2));
                end

                % Energy dissipated by CH
                if min_dis > 0
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * 4000);
                    enrgy_res(C(min_dis_cluster).id, r) = S(C(min_dis_cluster).id).E;
                    Energy_disp(r) = Energy_disp(r) + ((ERX + EDA) * 4000);
                end

                packets_TO_CH = n - dead - countCHs;
                PACKETS_TO_CH(r) = packets_TO_CH;
            end
        end
    end

    % Total energy left in the network
    for i = 1:n
        if S(i).E > 0
            energy_(r) = energy_(r) + S(i).E;
        end
    end

    % Average energy of alive nodes
    avg_r(r) = energy_(r) / ALIVE_NODE(r);

    % Packet Delivery Ratio (PDR)
    if (PACKETS_TO_BS(r) + PACKETS_TO_CH(r)) > 0
        PDR(r) = (PACKETS_TO_BS(r) / (PACKETS_TO_BS(r) + PACKETS_TO_CH(r))) * 100;
    end

end

% Plotting the results
% Network Lifetime (Figure 1): Number of alive nodes per round
figure(2);
plot(1:last, ALIVE_NODE(1:last));
title('Network Lifetime');
xlabel('Rounds');
ylabel('Number of Alive Nodes');

% Average Node Energy (Figure 3): Average energy (normalized) of nodes per round
figure(3);
plot(1:last, avg_r(1:last) / Eo);
title('Average Node Energy');
xlabel('Rounds');
ylabel('Average Energy (Normalized)');

% Packet Delivery Ratio (PDR) (Figure 4): Percentage of packets successfully delivered to the base station per round
figure(4);
plot(1:last, PDR(1:last));
title('Packet Delivery Ratio (PDR)');
xlabel('Rounds');
ylabel('PDR (%)');

% Energy Consumption (Figure 5): Total energy consumed by the network (normalized) per round
figure(5);
plot(1:last, Energy_disp(1:last) / n);
title('Energy Consumption');
xlabel('Rounds');
ylabel('Total Energy Consumed (Normalized)');

% Cluster Head Distribution (Figure 6): Number of cluster heads per round
figure(6);
plot(1:last, CLUSTERHS(1:last));
title('Cluster Head Distribution');
xlabel('Rounds');
ylabel('Number of Cluster Heads');

% Throughput (Packets Delivered) (Figure 7): Total number of packets delivered to the base station per round
figure(7);
plot(1:last, cumsum(PACKETS_TO_BS(1:last)));
title('Throughput (Packets Delivered)');
xlabel('Rounds');
ylabel('Total Packets Delivered to BS');

% Dead Nodes vs Rounds (Figure 8): Number of dead nodes per round
figure(8);
plot(1:last, DEAD(1:last));
title('Dead Nodes vs Rounds');
xlabel('Rounds');
ylabel('Number of Dead Nodes');

% Packets Transmitted to BS vs Rounds (Figure 9): Number of packets transmitted to the base station per round
figure(9);
plot(1:last, PACKETS_TO_BS(1:last));
title('Packets Transmitted to BS vs Rounds');
xlabel('Rounds');
ylabel('Packets Transmitted to BS');

% Remaining Energy of Nodes (Figure 10): Remaining energy of each node at the end of simulation
figure(10);
bar([S.E]);
title('Remaining Energy of Nodes');
xlabel('Node Index');
ylabel('Remaining Energy (Joules)');

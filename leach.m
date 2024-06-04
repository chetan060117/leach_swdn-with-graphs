       %
%                                                                      %                                  
%     "A clustering algorithm based on energy information              %
%                      and cluster heads                               %
%            expectation for wireless sensor networks "                %                                                             
%                                                                      %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBMITTED BY-                                                        %
%                SE20UCSE071- K.Chetan                                 %
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
Eo = 0.7;  % Increase the initial energy from 0.5 to 0.7
% Eelec = ETX = ERX
ETX = 40 * 0.000000001; % Reduced from 50
ERX = 40 * 0.000000001; % Reduced from 50
% Transmit Amplifier types
Efs = 8 * 0.000000000001; % Reduced from 10
Emp = 0.0013 * 0.000000000001;
% Data Aggregation Energy
EDA = 4 * 0.000000001; % Reduced from 5

% Values for Heterogeneity
% Percentage of nodes that are advanced
m = 0.1; % Increased from 0.05
% \alpha
a = 0.1;

% Maximum number of rounds
rmax = 30000;

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
   
    candidate_ch_probability = i;
    % Random Election of Normal Nodes
    if (candidate_ch_probability >= m * n + 1) 
        S(i).E = Eo;
        S(i).ENERGY = 0;
        plot(S(i).xd, S(i).yd, 'o');
        hold on;
    end
    % Random Election of Advanced Nodes
    if (candidate_ch_probability < m * n + 1)  
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
Throughput = zeros(1, rmax);

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

% Initialize sliding window parameters
window_size = 15; % Adjust window size as needed
energy_window = zeros(1, window_size);
window_index = 1;

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

    % Print number of alive and dead nodes
    fprintf('Round %d: Alive Nodes: %d, Dead Nodes: %d\n', r, alive, dead);

    % When the first node dies
    if dead == 1 && flag_first_dead == 0
        first_dead = r;
        flag_first_dead = 1;
    end

    countCHs = 0;
    cluster = 1;

    % Cluster Head Selection
    for i = 1:n
        if S(i).E > 0
            % Average energy
            energy_moy(i, r) = (1/r) * sum(enrgy_res(i, 1:r));
            % Generate a random number within the interval [0, 1]
            temp_rand = rand;

            if S(i).G <= 0
                % Update sliding window with the average energy of non-CH nodes
                energy_window(window_index) = energy_moy(i, r);
                window_index = mod(window_index, window_size) + 1;
                avg_energy_window = mean(energy_window);

                % Election of Cluster Heads based on dynamic probability
                if temp_rand <= (p * (S(i).E / Eo) / (1 - p * mod(r, round(1/p))) * (S(i).E / avg_energy_window))
                    countCHs = countCHs + 1;
                    packets_TO_BS = packets_TO_BS + 1;
                    PACKETS_TO_BS(r) = packets_TO_BS;
                    
                    S(i).type = 'C';
                    S(i).G = 1/p;
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
                    distance;
                    if distance > do
                        energy_consumed = ETX * 4000 + Emp * 4000 * (distance^4);
                    else
                        energy_consumed = ETX * 4000 + Efs * 4000 * (distance^2);
                    end
                    S(i).E = S(i).E - energy_consumed;
                    
                    % Ensure energy doesn't drop below zero
                    if S(i).E < 0
                        S(i).E = 0;
                    end
                    
                    fprintf('Node %d, Round %d: CH with remaining energy %.4f J\n', i, r, S(i).E);
                end
            end
        end
    end

    CLUSTERHS(r) = countCHs;
    hold on;

    % If there are fewer than 1% cluster heads, select cluster heads randomly
    if countCHs < 0.01 * n
        for i = 1:n
            if S(i).E > 0
                S(i).G = 0;
            end
        end
    end

    % Reset clusters for re-election
    for i = 1:n
        if S(i).E > 0
            if S(i).type == 'N'
                % Resetting to cluster head probability
                S(i).G = 0;
            end
        end
    end

    % Cluster Formation and Transmission
    for i = 1:n
        if S(i).type == 'N' && S(i).E > 0
            % Find the closest cluster head
            min_dis = inf;
            min_dis_cluster = 0;

            for c = 1:countCHs
                temp = sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2);
                if temp < min_dis
                    min_dis = temp;
                    min_dis_cluster = c;
                end
            end

            if min_dis_cluster ~= 0
                % Calculate the energy dissipated for transmission to CH
                if min_dis > do
                    energy_consumed = ETX * 4000 + Emp * 4000 * (min_dis^4);
                else
                    energy_consumed = ETX * 4000 + Efs * 4000 * (min_dis^2);
                end
                S(i).E = S(i).E - energy_consumed;

                % Ensure energy doesn't drop below zero
                if S(i).E < 0
                    S(i).E = 0;
                end

                % Cluster Head receives the data
                energy_consumed = ERX * 4000;
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - energy_consumed;

                % Ensure energy doesn't drop below zero for CH
                if S(C(min_dis_cluster).id).E < 0
                    S(C(min_dis_cluster).id).E = 0;
                end

                % Cluster Head aggregates the data and sends to base station
                energy_consumed = EDA * 4000;
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - energy_consumed;

                % Ensure energy doesn't drop below zero for CH after aggregation
                if S(C(min_dis_cluster).id).E < 0
                    S(C(min_dis_cluster).id).E = 0;
                end

                PACKETS_TO_CH(r) = PACKETS_TO_CH(r) + 1;
                Throughput(r) = Throughput(r) + 1; % Increment throughput
            end
        end
    end

    % Calculate and display average remaining energy
    avg_remaining_energy = mean([S.E]);
    fprintf('Round %d: Average remaining energy: %.4f J\n', r, avg_remaining_energy);

    % Update energy_disp for the current round
    Energy_disp(r) = sum([S.E]);

    % Calculate Packet Delivery Ratio (PDR)
    if PACKETS_TO_CH(r) + PACKETS_TO_BS(r) > 0
        PDR(r) = (PACKETS_TO_CH(r) + PACKETS_TO_BS(r)) / (PACKETS_TO_CH(r) + PACKETS_TO_BS(r) + DEAD(r));
    else
        PDR(r) = 0;
    end
end

% Plot results
figure;
plot(1:last, DEAD(1:last), 'r');
xlabel('Rounds');
ylabel('Number of Dead Nodes');
title('Number of Dead Nodes over Time');

figure;
plot(1:last, Energy_disp(1:last), 'b');
xlabel('Rounds');
ylabel('Total Energy Dispersion');
title('Total Energy Dispersion over Time');

% Plot Network Lifetime
figure;
plot(1:last, ALIVE_NODE(1:last), 'g');
xlabel('Rounds');
ylabel('Number of Alive Nodes');
title('Network Lifetime');

% Plot Packet Delivery Ratio (PDR)
figure;
plot(1:last, PDR(1:last) * 100, 'm');
xlabel('Rounds');
ylabel('Packet Delivery Ratio (%)');
title('Packet Delivery Ratio over Time');

% Plot Total Energy Consumption
figure;
plot(1:last, Eo * n - Energy_disp(1:last), 'k');
xlabel('Rounds');
ylabel('Total Energy Consumption (J)');
title('Total Energy Consumption over Time');

% Plot Cluster Head Distribution
figure;
plot(1:last, CLUSTERHS(1:last), 'c');
xlabel('Rounds');
ylabel('Number of Cluster Heads');
title('Cluster Head Distribution over Time');

% Plot Throughput
figure;
plot(1:last, Throughput(1:last), 'b');
xlabel('Rounds');
ylabel('Total Data Packets Successfully Delivered');
title('Throughput over Time');

% Plot histogram of remaining energy
figure;
histogram([S.E]);
xlabel('Nodes');
ylabel('Remaining Energy (J)');
title('Remaining Energy of Nodes at the End of Simulation');

disp('Simulation Completed.');
%%
% CV2XMode3 is the script of the implementation of the analytical models
% of the communication performance of C-V2X Mode 3 described in the
% following paper:
%
%   Paper title: LTE-V2X Mode 3 scheduling based on adaptive spatial reuse of radio resources
%   Paper authors: Daniel Sempere-García, Miguel Sepulcre and Javier Gozalvez
%
% This model quantifies the PDR (Packet Delivery Ratio) that could be
% achieved with our proposed scheme as a function of the distance between
% the transmitter and the receiver. In order to model the PDR, the
% following four mutually exclusive errors present in C-V2X are quantified:
%   1)	Errors due to half-duplex transmissions (HD)
%   2)	Errors due to a received signal power below the sensing power threshold (SEN)
%   3)	Errors due to propagation effects (PRO)
%   4)	Errors due to packet collisions (COL)
%
% CV2XMode3.m is the main script you have to run to get the PDR curve as a
% function of the distance for a given set of parameters, and the
% probability of each of the four transmission errors.
%
% Input parameters (shall be set manually in the beginning of the code):
%    alpha: traffic density in veh/km. Values tested: 120.
%    lambda: pNumber of packets transmitted per second per vehicle. Values tested: 10, 20 and 50.
%    Psen: sensing threshold (dBm). Values tested: -90.5
%    Pt: transmission power in dBm. Values tested: 23.
%    S: number of sub-channels per sub-frame. Values tested: 2 and 4.
%    B: packet size in bytes. Values tested: 190.
% 
% Output metrics:
%    PDR: Packet Delivery Ratio for different Tx-Rx distances 
%    deltaHD: probability of packet loss due to half-duplex transmissions for different Tx-Rx distances
%    deltaSEN: probability of packet loss due to a received signal power below the sensing power threshold for different Tx-Rx distances
%    deltaPRO: probability of packet loss due to propagation effects for different Tx-Rx distances
%    deltaCOL: probability of packet loss due to packet collisions for different Tx-Rx distances

tic
global distance noise Psen step_dB Pt B MCS

% Set to 1 to get debug figures:
DEBUG = 0;
DEBUG_HD = 0;
DEBUG_reuse = 0;

% Input parameters:
alpha = 120;                % Traffic density in vehicles/km
lambda = 10;                % Number of packets transmitted per second per vehicle
Psen = -90.5;               % Sensing threshold (dBm)
Pt = 23;                    % Transmission power in dBm
S = 4;                      % Number of sub-channels per sub-frame
B  = 190;                   % Packet size in bytes

% Parameters to be calculated based on the input ones:
if S == 2
    MCS = 7;                % MCS of each packet
    RBs = 14;               % Number of RBs per packet
elseif S == 4
    MCS = 9;                % MCS of each packet
    RBs = 10;               % Number of RBs per packet
end
distance = 0:25:1000;       % Distances to the Tx at which the PDR will be calculated
Sp = 1;                     % Number of sub-channels per packet
step_dB = 0.1;              % Discrete steps to compute the PDF of the SNR and SINR (dB)
N = S * 1000/lambda;        % Total number of resources in 1000ms
ivs = 1000/alpha;           % Inter-vehicle spacing in meters
K = N/floor(S/2);           % Parameter used as an index of the resources of the pool

%% Probability of packet loss due to half duplex effect (deltaHD_pre):
% The proposed scheduling scheme will try to assign resources in the same
% subframe to two vehicles if they are located at a distance equal to dHD,
% but in reality the distance between them will sometimes be lower or
% higher. To embed the mobility of the vehicles in the analytical model, we
% propose to use triangle functions around the multiples of dHD instead of
% unit delta functions. As a result, the probability that vehicle vr cannot
% receive a packet transmitted by vehicle vt due to the HD effect is
% expressed as the sum of multiple triangular functions. The total
% probability of each triangle is equal to 1, which implies that each
% triangle represents a vehicle in HD with the tx node:
d = 0;
p = 0;
numberOfMultiples = 10;                     % Number of multiples of the HD optimum distance considered
d = (0:(1.5*numberOfMultiples*N-1))*ivs;    % Vector with the distances at which there is a vehicle
% Based on the number of multiples, define the centres of the triangles:
centers = zeros(1,numberOfMultiples);
for selectedMultiple = 1:numberOfMultiples
    centers(selectedMultiple) = selectedMultiple*K/2;
end

xcenters = d;                                           % Copy of the distances vector
d = [-fliplr(xcenters(2:length(xcenters))) xcenters];   % Replicate the distances vector to the left (negative values)
pdfs = zeros(numberOfMultiples,length(xcenters));       % Matrix to store the triangular distributions for each multiple

for i = 1:numberOfMultiples
    d1_hd = xcenters(centers(i)+1-(K)/2);               % Lower value of the triangle
    dpeak_hd = xcenters(centers(i)+1);                  % Central value of the triangle
    d2_hd = xcenters(centers(i)+1+(K)/2);               % Maximum value of the triangle
    
    pd = makedist('Triangular','a',d1_hd,'b',dpeak_hd,'c',d2_hd);   % Generate the triangular distribution associated with the selected multiple
    pdfs(i,:) = pdf(pd,xcenters);                       % Get the pdf of the triangular distribution
    pdfs(i,:) = pdfs(i,:)*ivs;                          % Multiply the triangular curve so that the weighted values sum a total probability equal to 1
    
    % Debug figures:
    if DEBUG_HD
        s = sum(pdfs(i,:));
        nonZeroValuesIndexes = find(pdfs(i,:));
        nonZeroValues = length(find(pdfs(i,:)));
        product = prod((1-pdfs(i,nonZeroValuesIndexes)));
        
        figure;
        bar(xcenters,pdfs(i,:));grid on
        xlabel('Distance (m)');ylabel('%')
        title(['Triangular curve (a = ' num2str(d1_hd) ' / b = ' num2str(dpeak_hd) ' / c = ' num2str(d2_hd) ' / sum = ' num2str(sum(pdfs(i,:))) ' / prod = ' num2str(prod((1-pdfs(i,nonZeroValuesIndexes)))) ')'])
        xlim([0 5000])
    end
end
% Get the sum of all the triangular distributions:
p = sum(pdfs);

% Debug figures:
if DEBUG_HD
    figure;
    bar(xcenters,p);grid on
    xlabel('Distance (m)');ylabel('%')
    title('Sum of triangular curves')
    xlim([0 5000])
end

% Replicate the probabilities vector to the left (negative values):
p = [fliplr(p(2:length(p))) p];

% Debug figures:
if DEBUG_HD
    figure;
    bar(d,p);grid on
    xlabel('Distance (m)');ylabel('%')
    title('Sum of triangular curves replicated')
    xlim([-5000 5000])
end

% Get the value of deltaHD for each value of distance:
for k=1:length(distance)
    d2hd = (abs(distance(k)-(d)));
    [d_min,i] = min(d2hd);
    deltaHD_pre(k) = mean(p((i-0):(i+0)));
end

% Debug figures:
if DEBUG
    figure;
    bar(d,p)
    hold on;
    plot(distance,deltaHD_pre)
    xlim([0 1000]);ylim([0 1])
    xlabel('Distance (m)'); ylabel('Probability of HD');
    legend('p','deltaHD\_pre')
end

%% Probability of packet loss due to signal level below threshold (deltaSEN_pre):
[PL, std_dev] = get_PL_SH(distance);                                        % Obtain pathloss and shadowing for different Tx-Rx distances
deltaSEN_pre = 0.5 * ( 1 - erf( (Pt - PL - Psen)./(std_dev*sqrt(2)) ) );    % Calculate deltaSEN

%% Probability of packet loss due to propagation (deltaPRO_pre):

noise = -95 - 10*log10(50/RBs); % Noise corresponding to the DATA field of each message. Assumes a noise figure of 9dB and 10MHz channel (background noise of -95dBm). The total number of RBs in 10MHz is 50

[SNR, PDF_SNR] = get_SINRdistribution(Pt - PL, -201, std_dev, 3, noise, Psen, step_dB); % Obtain the PDF of the SNR experienced by Rx (without interference)
deltaPRO_pre = get_BLER(SNR, PDF_SNR, B, MCS, step_dB);                              % Get deltaPRO_pre

%% Probability of packet loss due to collision (deltaCOL_pre):
% Our proposal assumes that vehicles that reuse a vehicle are located at a
% distance around the reuse optimum distance and its positive integer
% multiples. The probability of experiencing a collision is considered as
% several triangular distribution curves centered at the reuse optimum
% distance and its multiples. The total probability of each triangle is
% equal to 1, which implies that each triangle represents an interfering
% vehicle. In this method, p_SIM is individual for every multiple of the
% reuse optimum distance. To get an individual deltaCOL_interferer value
% for each multiple, each individual p_SIM and p_int are multiplied, and
% then the resulting values are summed so as to get the deltaCOL_interferer
% value.

% Calculate the probability pSIM:
d = 0;
p_SIM = 0;   

numberOfMultiples = 5;                      % Number of multiples of the reuse optimum distance considered
d = (0:(1.5*numberOfMultiples*N-1))*ivs;    % Vector with the distances at which there is a vehicle
% Based on the number of multiples, define the centres of the triangles:
centers = zeros(1,numberOfMultiples);
for selectedMultiple = 1:numberOfMultiples
    centers(selectedMultiple) = selectedMultiple*N;
end

xcenters = d;                                           % Copy of the distances vector
d = [-fliplr(xcenters(2:length(xcenters))) xcenters];   % Replicate the distances vector to the left (negative values)
pdfs = zeros(numberOfMultiples,length(xcenters));       % Matrix to store the triangular distributions for each multiple
p_SIM = zeros(numberOfMultiples,length(d));             % Matrix to store the triangular distributions replicated for each multiple

for i = 1:numberOfMultiples
    d1_r = xcenters(centers(i)+1-N/2);                  % Lower value of the triangle
    dpeak_r = xcenters(centers(i)+1);                   % Central value of the triangle
    d2_r = xcenters(centers(i)+1+N/2);                  % Maximum value of the triangle

    pd = makedist('Triangular','a',d1_r,'b',dpeak_r,'c',d2_r);  % Generate the triangular distribution associated with the selected multiple
    
    pdfs(i,:) = pdf(pd,xcenters);                       % Get the pdf of the triangular distribution
    pdfs(i,:) = pdfs(i,:)*ivs;                          % Multiply the triangular curve so that the weighted values sum a total probability equal to 1
    
    p_SIM(i,:) = [fliplr(pdfs(i,2:length(pdfs(i,:)))) pdfs(i,:)];   % Replicate the probabilities vector to the left (negative values)
    
    % Debug figures:
    if DEBUG_reuse
        s = sum(pdfs(i,:));
        nonZeroValuesIndexes = find(pdfs(i,:));
        nonZeroValues = length(find(pdfs(i,:)));
        product = prod((1-pdfs(i,nonZeroValuesIndexes)));
        
        figure;
        bar(xcenters,pdfs(i,:));grid on
        xlabel('Distance (m)');ylabel('%')
        title(['Triangular curve (a = ' num2str(d1_r) ' / b = ' num2str(dpeak_r) ' / c = ' num2str(d2_r) ' / sum = ' num2str(sum(pdfs(i,:))) ' / prod = ' num2str(prod((1-pdfs(i,nonZeroValuesIndexes)))) ')'])
        xlim([0 5000])
        
        figure(99);
        bar(d,p_SIM(i,:));grid on;hold on
        xlabel('Distance (m)');ylabel('%')
        title('Triangular curves replicated')
        xlim([-5000 5000])
    end
    
end
% Call 'get_COL' function to get deltaCOL_pre:
[deltaCOL_pre,p_aux,p_int] = get_COL(d, p_SIM, deltaPRO_pre, DEBUG_reuse);
        
% Debug figures:
if DEBUG_reuse
    figure;
    bar(d,sum(p_SIM))
    hold on;
    plot(distance,deltaCOL_pre)
    xlim([0 1000])
    xlabel('Distance (m)'); ylabel('Probability of collision');
    legend('p','deltaCOL_pre\')
end

%% Figure that represents the values of each of the calculated errors 
% before calculating the final ones:
figure;hold on;grid on
plot(distance,deltaHD_pre)
plot(distance,deltaSEN_pre)
plot(distance,deltaPRO_pre)
plot(distance,deltaCOL_pre)
legend('deltaHD_pre','deltaSEN_pre','deltaPRO_pre','deltaCOL_pre')

% Calculate the final value of each error taking into account that SEN
% errors exclude HD errors, PRO errors exclude HD and SEN errors, and that
% COL errors exclude HD, SEN and PRO errors:
deltaHD  = deltaHD_pre ;
deltaSEN = deltaSEN_pre .* (1 - deltaHD_pre);
deltaPRO = deltaPRO_pre .* (1 - deltaHD_pre) .* (1 - deltaSEN_pre);
deltaCOL = deltaCOL_pre .* (1 - deltaHD_pre) .* (1 - deltaSEN_pre) .* (1 - deltaPRO_pre);

% Figure that represents the values of each of the calculated errors
figure;hold on;grid on
plot(distance,deltaHD)
plot(distance,deltaCOL)
plot(distance,deltaSEN)
plot(distance,deltaPRO)
legend('deltaHD','deltaCOL','deltaSEN','deltaPRO')

% Figure that represents the values of each of the calculated errors mixing
% in the same curve the SEN and PRO errors.
figure;hold on;grid on
plot(distance,100*deltaHD)
plot(distance,100*deltaCOL)
plot(distance,100*deltaPRO+100*deltaSEN)
legend('deltaHD','deltaCOL','deltaPRO+SEN')

% Calculate and plot the PDR:
PDR = 1 - deltaSEN - deltaPRO - deltaCOL - deltaHD;

figure; hold on; grid on
%     plot(distance,100*PDR,'LineWidth',2);
plot(distance,PDR,'LineWidth',2);

toc

%% ###############################################

function [deltaCOL,p_aux,p_int] = get_COL(distance_int_to_tx, p_SIM, deltaPRO_pre, DEBUG_reuse);
% Function to get deltaCOL.
global distance noise Psen step_dB Pt B MCS

if DEBUG_reuse
    figure;
    subplot(3,1,1)
    plot(distance_int_to_tx,p_SIM,'k');hold on
end

deltaCOL_interferer = zeros(length(distance),size(p_SIM,1));    % Matrix to store the

for d=1:length(distance)
    
    [PL_E_R(d), std_dev_E_R(d)] = get_PL_SH(distance(d));       % Calculate the pathloss and shadowing for a given Tx-Rx distance following the Winner+ B1 propagation model
    
    % Calculate probability of collision for each interfering vehicle:
    distance_int_to_rx = distance_int_to_tx + distance(d);  % Distances from all the interfering vehicles to the transmitting vehicle.
    
    for i=1:length(distance_int_to_rx)
        [PL_I_R , std_dev_I_R] = get_PL_SH(distance_int_to_rx(i));   % Pathloss and shadowing for interf and rx
        
        Pi_dB = Pt - PL_I_R; % Average received interference
        if deltaPRO_pre(d) > 0.999
            p_int(d,i) = 0;  % If the proability of packet loss due to propagation is 1
        else
            [SINR, PDF_SINR] = get_SINRdistribution( Pt-PL_E_R(d) , Pi_dB , std_dev_E_R(d) , std_dev_I_R , noise , Psen , step_dB); % PDF of the SINR experienced by the receiving vehicle
            p_SINR(d,i) = get_BLER ( SINR , PDF_SINR , B , MCS , step_dB );     % Probability that the receiver receives a packet with error due to low SINR. Equation (17)
            p_int(d,i) = ( p_SINR(d,i) - deltaPRO_pre(d) ) / ( 1 - deltaPRO_pre(d) );  % Probability that the interference generated on the receiver provokes that the packet transmitted by vt cannot be correctly received at vr.  Equation (18)
        end
    end
    % For each multiple, calculate deltaCOL_interferer as the sum of the
    % elements obtained multiplying p_int and p_SIM:
    for selectedMultiple = 1:size(p_SIM,1)
        p_aux(d,:) = p_int(d,:).*p_SIM(selectedMultiple,:);
        deltaCOL_interferer(d,selectedMultiple) = min(sum(p_aux(d,:)),1);
    end
    
    % Calculate deltaCOL:
    deltaCOL(d) = 1 - prod(1 - deltaCOL_interferer(d,:));
    
    if DEBUG_reuse
        subplot(3,1,1)
        plot(distance_int_to_tx,p_int(d,:));hold on;grid on
        subplot(3,1,2)
        plot(distance_int_to_tx,p_aux(d,:));hold on;grid on
        subplot(3,1,3)
        plot(distance(d),deltaCOL(d),'*');hold on;grid on
    end
end

return
end
%% ###############################################

function [ PL , std_dev ] = get_PL_SH( distance );
% Function to calculate the pathloss and shadowing for a given set of Tx-Rx
% distances following the Winner+ B1 propagation model.

% Parameters of the radio propagation model:
fc = 5.91e9;                % Carrier frequency (Hz)
hBS = 1.5;                  % Transmitter antenna height (m)
hMS = 1.5;                  % Receiver antenna height (m)
environmentHeight = 0;      % Average environmental height (m)
distance = abs(distance);

c = 3e8;
dBP = 4 * (hBS-environmentHeight) * (hMS-environmentHeight) * fc / c; % breakpoint distance

% Avoid errors for very small distances:
i = find(distance < 3);
distance(i) = 3;

% Calculate pathloss for distances lower than the breakpoint distance:
i = find(distance < dBP);
PL(i) = 22.7*log10(distance(i)) + 27 + 20*log10(fc/1e9);
std_dev(i) = 3;    % Standard deviation

% Calculate pathloss for distances higher than the breakpoint distance:
i = find(distance >= dBP);
PL(i) = 40*log10(distance(i)) + 7.56 - 17.3*log10(hBS-environmentHeight) - 17.3*log10(hMS-environmentHeight) + 2.7*log10(fc/1e9);
std_dev(i) = 3;    % Standard deviation

% Compares obtained pathloss with free-space pathloss:
PLfree = 20*log10(distance) + 46.4 + 20*log10(fc*1e-9 / 5);
i = find(PLfree > PL);
PL(i) = PLfree(i);

return
end

%% #####################################################

function [ SINR , PDF_SINR ] =  get_SINRdistribution( Pr_dBm_avg , Pi_dBm_avg , std_dev_Pr , std_dev_Pi , noise_dBm , sensingThreshold , step_dB );
% Function to calculate the PDF of the SNR or SINR at the receiver based on
% average received power and interference levels, noise and other
% parameters.

x = -200:step_dB:200; % Wide range of values in dB to build the PDF

for i=1:length(Pr_dBm_avg)
    
    distrib_Pr = normpdf(x,Pr_dBm_avg(i),std_dev_Pr(i));    % PDF of the received signal
    ind = find( x < sensingThreshold);
    distrib_Pr(ind) = 0;                                    % Remove values below the sensing threshold
    distrib_Pr = distrib_Pr / sum(distrib_Pr) / step_dB;    % Normalize so that the integral between -inf and +inf is equal to 1
    
    if Pi_dBm_avg < -200
        % If there is no interference:
        distrib_Pi_noise = zeros(1,length(x));
        distrib_Pi_noise( round( (noise_dBm-x(1))/step_dB ) +1 ) = 1 / step_dB;
    else
        % If there is interference, compute the PDF of the SINR:
        aux = length( find(x <= noise_dBm) );   % Remove values below noise, since Pi+noise will never be lower than noise
        noise = 10^(noise_dBm/10);              % Noise power in linear units
        distrib_Pi_noise = 10.^(x(aux+1:end)/10)./((10.^(x(aux+1:end)/10)-noise)*std_dev_Pi*sqrt(2*pi)) .* exp( -(10*log10(10.^(x(aux+1:end)/10)-noise)-Pi_dBm_avg).^2/(2*std_dev_Pi^2) ); % PDF of interference + noise
        distrib_Pi_noise = [ zeros(1,aux) distrib_Pi_noise];                        % Null probability for values below the noise
        distrib_Pi_noise = distrib_Pi_noise / sum(distrib_Pi_noise) / step_dB;   	% Normalize so that the integral between -inf and +inf is equal to 1
    end
    
    % Calculate PDF of SINR through the cross-correlation of the PDF of the received power and the PDF of the interference+noise:
    PDF_SINR(i,:) = xcorr(distrib_Pr,distrib_Pi_noise);
    PDF_SINR(i,:) = PDF_SINR(i,:) / sum(PDF_SINR(i,:)) / step_dB;    % Normalize so that the integral between -inf and +inf is equal to 1
    
    
    % Adapt the range of the x axes to the values provided by the xcorr function:
    SINR(i,:) = min(x)*2:step_dB:max(x)*2;
end

return
end

%% #####################################################

function [ avg_BLER ] = get_BLER( SINR , PDF , B , MCS , step_dB);
% Function to calculate the average BLER experienced given the PDF of the
% SINR at the receiver and a given coding scheme.

vector_SNR_paper  =[-500 -4 -2 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 500];

if B == 190
    % 190 Bytes
    if MCS == 7
        % 14RB: QPSK 0.5 - NOT USED HERE
        vector_BLER_paper = [1 1 0.95 0.8 0.5 0.3 0.107 0.04 0.01 0.003 0.001 0.0003 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001] ;
    elseif MCS == 9
        % 10RB: QPSK 0.7
        vector_BLER_paper = [1 1 1.0 0.95 0.7 0.5 0.21 0.1 0.04 0.012 0.004 0.0012 0.0006 0.0004 0.0003 0.0003 0.0003 0.0003 0.0003 0.0003];
    end
elseif B == 300
    % 300 Bytes
    if MCS == 7
        % 20RB: QPSK 0.5
        vector_BLER_paper = [1 1 0.95 0.7 0.5 0.2 0.07 0.01 0.0015 0.0005 0.0002 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001] ;
    elseif MCS == 9
        % 16RB: QPSK 0.7 - NOT USED HERE
        vector_BLER_paper = [1 1 1.0000 0.95 0.8 0.6 0.3 0.12 0.04 0.01 0.003 0.001 0.0006 0.0004 0.0003 0.0003 0.0003 0.0003 0.0003 0.0003];
    end
end

for i=1:size(SINR,1)
    BLER_interp = interp1(vector_SNR_paper , vector_BLER_paper , SINR(i,:) , 'linear');   % BLER vector for a given SINR vector for a given Tx-Rx distance
    avg_BLER(i) =  PDF(i,:) * BLER_interp' * step_dB;   % Equation (13) or Equation (17) depending on whether the input is the PDF(SINR) or PDF(SNR)
end

return
end

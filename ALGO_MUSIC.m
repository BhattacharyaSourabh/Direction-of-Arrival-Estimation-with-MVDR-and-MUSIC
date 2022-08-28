%Sourabh Bhattacharya
%MUSIC algorithm
%21EE64R18
%Defining the constants
nat_freq=(5*(10^9));
wavelength=(3*(10^8)/nat_freq);
separation=wavelength/2;
datapoint=4500;
num_ele=10;
tot_ang=[-30 -50 -20]*pi/180; 
num_ang=length(tot_ang);
SNR=50;
%%
%Preparing the simulation 
%Preparing the noise eigen vectors and source eigen vectors
steer_vec=exp(-1i*2*pi*(separation/wavelength)*(0:num_ele-1)'*sin(tot_ang));
signal=randn(num_ang,datapoint); 
signal=exp(1).^(1i*2*pi*nat_freq*signal);
rand_noise=(randn(num_ele,datapoint)); 
X=steer_vec*signal+rand_noise; 
X_new=X*ctranspose(X); 
[eig_vect ,eig_val]=eig(X_new); 
[eig_val,idx]=sort(diag(eig_val),1,'descend'); 
eig_vect=eig_vect (:,idx);
sig_vect=eig_vect (:,1:num_ang); 
noise_vect=eig_vect(:,num_ang+1:num_ele); %Get the noise eigenvectors
%%
%Graph of noise, signal and eigen values
%Eigen values have an increasing tendency with index
figure;
imagesc(abs(X_new))
%grid on
title('AUTOCORRELATION GRAPH')
xlabel('NUMBER OF NODES')
ylabel('NUMBER OF NODES')
plot((0:datapoint-1),rand_noise)
grid on
title('GRAPH OF NOISE')
xlabel('DATAPOINTS')
ylabel('NOISE')
plot(idx,eig_val)
grid on
title('EIGEN VALUE VS INDEX')
xlabel('INDEX')
ylabel('EIGEN VALUE')
plot((0:datapoint-1),signal)
grid on
title('GRAPH OF SIGNAL')
xlabel('DATAPOINTS')
ylabel('SIGNALS')
%%
% MUSIC algorithm
% Define angles at which MUSIC spectrum will be computed
tot_ang=(-90:0.1:90);
freq_spectrum=zeros(length(tot_ang),1);
%Compute steering vectors corresponding values in angles
steer_vect_tot=exp(-1i*2*pi*(separation/wavelength)*(0:num_ele-1)'*sin(tot_ang*pi/180));
for iter=1:length(tot_ang)
freq_spectrum(iter)=1/(ctranspose(steer_vect_tot(:,iter))*(noise_vect*ctranspose(noise_vect))*steer_vect_tot(:,iter));
end
%%
%PLOTTING THE FREQUENCY SPECTRUM
%PLOTTING TO SHOW PEAKS
figure
plot(tot_ang,abs(freq_spectrum))
grid on
title('FREQUENCY SPECTRUM MUSIC ALGO')
xlabel('ANGLES IN DEGREES')
ylabel('SPECTRUM VALUES FOR SPECIFIC DOA')


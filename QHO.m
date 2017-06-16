%% first 10 eigenstates within zero potential well 
% degeneracy 
hbar=1;
m=1;
omega=1;
pos=1;

% Jacobian Laplacian function for derivatives:

% Dimensions of laplacian function
N=100;
x=0:N-1;

%Laplacian function NxN with non zero diagonals and off diagonals
Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1) ...
+ diag(ones((N-1),1),-1));%/(dx^2)
% Next modify Lap so that it is consistent with f(0) = f(L) = 0
Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; % So that f(0) = 0
Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0;% So that f(L) = 0

% Hamiltonian 
% Define Hamiltonian for the system. Simple QHO 
H= (-(hbar.^2)/(2.*m)).*(Lap.^2) + 1/2*m*(omega^2)*(pos^2);  
eig(H); % find eigensolutions to Hamiltonian matrix 

% eigenenergies and eigenvectors 
[V,E]=eig(H);
V; % eigenstates
E; % eigenenergies 

% E given as square matrix with values on diagonal 
% diagonalise for eigenenergies 
En=diag(E); 

% compare with normalised eigenstates 
sqrt(2/L)*sin((10.*pi.*x)/L)';
V(:,10);

figure
for n=1:10
    t=0;
    evo=cos(t*En(n)/hbar);
    wf=evo.*(V(:,n));
    hold on
    plot(wf)
    drawnow % not required 
    F=getframe;
    pause(0.1)
end

figure
movie(F)


% Time evolution of individual eigenstate n 
for n=[1 5]
figure 
for t=0:0.1:2 % first n eigenstates
    evo=cos(t*En(n)/hbar); % = 1 for t=0 
    wf=(evo.*(V(:,n))); % wavefunction(x,t) 
    hold off 
    % split graphs for n 
    plot(wf) % probability density 
    axis([0,100,-0.5,0.5])
    F=getframe;
    pause(0.1)
end
end


sumwf=zeros(N,1); 
% Probability Density Function 
% probability density for the first n eigenstates 
% initial state t=0 to t
for n=1:5 
figure
for t=0:0.1:2% first n eigenstates
    evo=cos(t*En(n)/hbar);
    wf=(evo.*(V(:,n)));
    sumwf=sumwf+wf;
    hold off
    %figure % split graphs for n 
    plot(sumwf.^2)
    axis([0,100,0,0.5])
    drawnow  % not required 
    F=getframe;
    pause(0.1)
end
end

% Probability Density Function
% over an extended time 
for n=1:5 
figure
for t=0:10:100% first n eigenstates
    evo=cos(t*En(n)/hbar);
    wf=(evo.*(V(:,n)));
    sumwf=sumwf+wf;
    hold off
    %figure % split graphs for n 
    plot(sumwf.^2)
    axis([0,100,0,0.5])
    drawnow  % not required 
    F=getframe;
    pause(0.1)
end
end


sumwf=zeros(N,1);

% increasing number of eigenstates on potential of initial state 
% entanglement effects 
figure 
for t=0
for n=1:10 % firt n eigenstates
    evo=cos(t*En(n)/hbar); % = 1 for t=0 
    wf=(evo.*(V(:,n))); % wavefunction(x,t) 
    sumwf=sumwf+wf; % sum of eigenstates 
    hold off 
    %figure % split graphs for n 
    plot(sumwf.^2) % probability density 
    axis([0,100,0,0.5])
    F=getframe;
    pause(0.1)
end
end

% close all
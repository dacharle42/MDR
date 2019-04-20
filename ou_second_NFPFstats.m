function cou = ou_second_NFPFstats ( tmax, ntsteps, ncells, threshold, seed )

%*****************************************************************************80
%
%% OU_EULER applies the second-derivative method to the Ornstein-Uhlenbeck SDE.
%
%  Original code by John Burkardt was modified.   
%
%  Reference:
%
%    Daniel T. Gillespie,
%    Exact numerical simulation of the Ornstein-Uhlenbeck process 
%    and its integral,
%    Physical review. E, 54(2):2084-2091 · September 1996.
%
%  Parameters:
%
%    Input, real TMAX, NTSTEPS, NCELLS, THRESHOLD, SEED, the value of problem parameters.
%
%
%    Input, real TMAX, the final time.
%
%    Input, integer NTSTEPS, the number of time steps.
%
%    Input, integer NCELLS, the maximum number of cells.
%
%    Input, real THRESHOLD, the threshold for survival.
%
%    Input, integer SEED, a seed for the random number generator.
%
% for example: ou_second_NFPFstats(150, 10000, 1000, 0.24, 'shuffle');

thetaPF = 0.4;
thetaNF = 2.5;
muPF = 0.16;
muNF = 0.16;
CVPF=1.0;
CVNF=0.5;
sigmaPF=CVPF*muPF;
sigmaNF=CVNF*muNF;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'OU_EULER:\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Use an Euler method to approximate the solution of\n' );
  fprintf ( 1, '  the Ornstein-Uhlenbeck stochastic differential equation:\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '    d x(t) = theta * ( mu - x(t) ) dt + sigma dW\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  with initial condition x(0) = x0.\n' );

 
  if ( nargin < 1 )
    tmax = 3.0;
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Using default final time TMAX = %g\n', tmax );
  end 

  if ( nargin < 2 )
    n = 10000;
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Using default number of timesteps N = %d\n', n );
  end
  
  if ( nargin < 3 )
    ncells = 1;
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Using default number of cells ncells = %d\n', ncells );
  end
  
  if ( nargin < 4 )
    threshold = mu;
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Using default survival threshold = %g\n', threshold );
  end

  if ( nargin < 5 )
    seed = 123456789;
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Using default value of random SEED = %d\n', seed );
  end
  
  NMAX=ncells*10;   %maximum number of cells
  growthRate=0.5;

  figure;
  xlabel ( 't (days)', 'FontSize', 20 )
  ylabel ( 'N(t)', 'FontSize', 20, 'HorizontalAlignment', 'right' )
  set(gca,'FontSize',20);
  tstring = sprintf('threshold = %g', threshold);
  title( tstring,  'FontSize', 20);
  hold on;
  
% simulate NF first
for rep=1:10
%
%  Initialize the random number generator.
%  The following way to initialize the random number generator 
%  may not be available in older versions of MATLAB.
%
  rng ( seed )
%
%  Set the discrete time stepsize.
%
  dt = tmax / ntsteps;
%
%  Compute the Brownian increments.
%
  dw = sqrt ( dt ) * randn ( NMAX, ntsteps );
%
%  Carry out the Euler approximate integration process for NF.
%
  t = linspace ( 0, tmax, ntsteps + 1 );
  x = NaN*zeros ( NMAX, ntsteps + 1 );
  NcurrentNF = zeros ( 1 , ntsteps + 1 );

  x(1:ncells,1) = muNF + sigmaNF * randn(ncells,1);
  NcurrentNF(1)=ncells;
  for j = 1 : ntsteps
    %x(:,j+1) = x(:,j) + dt * thetaNF * ( muNF - x(:,j) ) + sigmaNF * dw(:,j);
    x(:,j+1) = muNF + ( x(:,j) - muNF ) + ( -(x(:,j)-muNF)*dt*thetaNF + sigmaNF * dw(:,j) ) * (1 - dt*thetaNF/2.0);
    indDead=find(x(:,j+1)<threshold);   %find dead cells
    x(indDead,j+1)=NaN;
    indNaN=find(isnan(x(:,j+1)));       % find spaces available for growth
    indAlive=find(~isnan(x(:,j+1)));    % find cells that can divide
    indGrow=find((rand(1,length(indAlive))<growthRate*dt));
    %nGrow=sum((rand(1,length(indAlive))>growthRate*dt));
    nGrow=length(indGrow);
    if(nGrow>length(indNaN))    %ensure there is no overgrowth
        nGrow=length(indNaN);
        indGrow=indGrow(1:nGrow);
    end
    indNew=indNaN(1:nGrow);
    x(indNew,j+1)=x(indAlive(indGrow),j+1);
    NcurrentNF(j+1)=sum(isfinite(x(:,j+1)));
  end
    
%   figure(1);
%   plot ( t, x, 'b-', 'LineWidth', 1 )
%   xlabel ( 't', 'FontSize', 16 )
%   ylabel ( 'X(t)', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right' )
%   title ( 'Euler solution of Ornstein-Uhlenbeck SDE', 'FontSize', 16 )
%   grid on; hold on;
  
  plot ( t, NcurrentNF, 'b-', 'LineWidth', 2 );
  
% simulate PF next
%
%  Initialize the random number generator.
%  The following way to initialize the random number generator 
%  may not be available in older versions of MATLAB.
%
  rng ( seed )
%
%  Set the discrete time stepsize.
%
  dt = tmax / ntsteps;
%
%  Compute the Brownian increments.
%
  dw = sqrt ( dt ) * randn ( NMAX, ntsteps );
%
%  Carry out the Euler approximate integration process for PF.
%
  t = linspace ( 0, tmax, ntsteps + 1 );
  x = NaN*zeros ( NMAX, ntsteps + 1 );
  NcurrentPF = zeros ( 1 , ntsteps + 1 );

  x(1:ncells,1) = muPF + sigmaPF * randn(ncells,1);
  NcurrentPF(1)=ncells;
  for j = 1 : ntsteps
    %x(:,j+1) = x(:,j) + dt * thetaPF * ( muPF - x(:,j) ) + sigmaPF * dw(:,j);
    x(:,j+1) = muPF + ( x(:,j) - muPF ) + ( -(x(:,j)-muPF)*dt*thetaPF + sigmaPF * dw(:,j) ) * (1 - dt*thetaPF/2.0);
    indDead=find(x(:,j+1)<threshold);   %find dead cells
    x(indDead,j+1)=NaN;
    indNaN=find(isnan(x(:,j+1)));       % find spaces available for growth
    indAlive=find(~isnan(x(:,j+1)));    % find cells that can divide
    indGrow=find((rand(1,length(indAlive))<growthRate*dt));
    %nGrow=sum((rand(1,length(indAlive))>growthRate*dt));
    nGrow=length(indGrow);
    if(nGrow>length(indNaN))    %ensure there is no overgrowth
        nGrow=length(indNaN);
        indGrow=indGrow(1:nGrow);
    end
    indNew=indNaN(1:nGrow);
    x(indNew,j+1)=x(indAlive(indGrow),j+1);
    NcurrentPF(j+1)=sum(isfinite(x(:,j+1)));
  end
  
  plot ( t, NcurrentPF, 'r-', 'LineWidth', 2 );
  
end
  
  
%
%  Plot the approximate solutions for NF and PF.
%
%   figure(1);
%   plot ( t, x, 'r-', 'LineWidth', 1 )

% 
%   filename = 'ou_euler_cells.png';
%   print ( '-dpng', filename )
%   fprintf ( 1, '\n' );
%   fprintf ( 1, '  Plot saved as "%s"\n', filename );
  
  legend('NF','PF');

  cou =  x;
  return;
end
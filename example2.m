%%%%%
%%  Numerical experiment 2:
%%    "Large" length-scale with varying number of nodes
%%%%

  % Parameter alpha
  a = 1/sqrt(2);
  
  % Parameters for the experiment
  Ns = [2:20];
  ells = [1.5 2 3 4 5];
  
  % Storage
  W = cell(length(ells),length(Ns));
  Wa = cell(length(ells),length(Ns));
  wmins = zeros(length(ells),length(Ns));
  wamins = zeros(length(ells),length(Ns));  
  
  for i = 1:length(ells);
    l = ells(i);
    % Kernel for this length-scale
    k = @(x,y) exp(-(x-y)^2/(2*l^2));
    kmean = @(x) (l^2 / (1+l^2))^(1/2) * exp( -norm(x)^2 /(2*(1+l^2)) );
    for j = 1:length(Ns);
      N = Ns(j);
      % Construct the nodes and compute the weights
      [X, wa] = kq_approx(l,a,N);
      w = kqw_symm(X, k, kmean);
      % Save      
      W{i,j} = w;
      Wa{i,j} = wa;
      wmins(i,j) = min(w);
      wamins(i,j) = min(wa);      
    end
    
  end
  
  % Compute some errors
  esq = zeros(length(ells),length(Ns));   % Norm of the weight-wise relative error
  emax = zeros(length(ells),length(Ns));  % Maximal absolute relative error
  eabs = zeros(length(ells),length(Ns));  % Norm of the absolute weight-wise error
  
  for i = 1:length(ells)
    for j = 1:length(Ns)
      esq(i,j) = sqrt(sum(((W{i,j}-Wa{i,j}) ./ W{i,j}).^2));
      eabs(i,j) = sqrt(sum((W{i,j}-Wa{i,j}).^2));
      emax(i,j) = max(abs((W{i,j}-Wa{i,j}) ./ W{i,j}));
    end
  end
  
  semilogy(Ns,esq,'LineWidth',1.5)
  title('Norm of the relative weight error vector')
  xlabel('N')
  legendCell = cellstr(num2str(ells', 'l = %.1f'));
  legend(legendCell)

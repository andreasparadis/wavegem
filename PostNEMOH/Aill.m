function [] = Aill(A,B,w)
  rho = 1025; L = 23; Nw = length(w);

  [ii, jj] = meshgrid([1:6], [1:6]);

  % Aill: For illustration purposes - normalized Added mass matrix

  for i = 1:6
    for j = 1:6
      for k = 1:Nw
        if and(le(i,3),le(j,3))
          Aill(i,j,k) = A(i,j,k)/(rho*L^3);
        elseif and(ge(i,4),ge(j,4))
          Aill(i,j,k) = A(i,j,k)/(rho*L^5);
        else
          Aill(i,j,k) = A(i,j,k)/(rho*L^4);
        endif
      endfor
    endfor
  endfor

  figure(); hold on;
  colormap('turbo');
  scatter([0:7],[0:7],1,".");
  for i = 1:6
    c = Aill(i,:,1);
    scatter(i*ones(1,6),[1:6],2000,c, "filled", "square")
  endfor
  title('\textbf{Visualization of A(0) matrix}')
  xticks([1 2 3 4 5 6])
  yticks([1 2 3 4 5 6])
  colorbar; %caxis([0 1e5]);
  plot_properties;
  grid off; axis square; axis ij;

  figure(); hold on;
  colormap('turbo');
  scatter([0:7],[0:7],1,".");
  for i = 1:6
    c = Aill(i,:,end);
    scatter(i*ones(1,6),[1:6],2000,c, "filled", "square")
  endfor
  title('\textbf{Visualization of} $\mathbf{A(\infty)}$ \textbf{matrix}')
  xticks([1 2 3 4 5 6])
  yticks([1 2 3 4 5 6])
  colorbar; %caxis([0 1e5]);
  plot_properties;
  grid off; axis square; axis ij;
end


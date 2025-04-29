function hi = plotcc(rho, pvalue, sampleVarNames, interestNames)
%PLOTCC Plot correlation coefficient matrix.

%narginchk(2, 4);
[ninterest, nsampleVar] = size(rho);

% plot PRCC value
hi = imagesc(rho, [-1 1]);
ha = gca();
set(ha, 'XTick', 1:nsampleVar, 'YTick', 1:ninterest);
if nargin() > 2 && ~isempty(sampleVarNames)
  set(ha, 'XTickLabel', sampleVarNames);
end
if nargin() > 3 && ~isempty(interestNames)
  set(ha, 'YTickLabel', interestNames);
end

m = size(get(gcf, 'colormap'), 1);
first = [0 0 1];
middle = [1 1 1];
last = [1 0 0];
cmap = interp1([1; m / 2; m], [first; middle; last], (1:m).');
colormap(cmap);
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
    set(0,'DefaultAxesFontName', 'Times New Roman')
    set(0,'DefaultAxesFontSize', 25)
    set(0,'DefaultTextFontname', 'Times New Roman')
    set(0,'DefaultTextFontSize', 25)
% add markers
for l = 1:nsampleVar
  for k = 1:ninterest
    if pvalue(k, l) < 0.01
      label = '**';
    elseif pvalue(k, l) < 0.05
      label = '*';
    else
      label = '';
    end
    text(l, k, label, 'FontWeight', 'bold', ...
      'Color', 'k', ...
      'BackgroundColor', 'none', ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
  end
end
title('*: (p-value) < 0.05, **: (p-value) < 0.01');

end


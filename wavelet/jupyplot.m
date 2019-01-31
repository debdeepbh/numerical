% a helping file to plot in jupyter
% for octave version 4.2
% in 4.2 this is resolved
% e.g. in Ubuntu 18.10, this is resolved
function jupyplot(texthere)
if ~exist("texthere")
    texthere = '';
end
    legend(texthere)
function varargout = unpack(S,varargin)
% unpack extracts selected fields of S and puts them in new variables
% EXAMPLE
% par.a = 3.5;
% par.b = [4;1.2];
% par.c = 'ale';
% [b,c] = unpack(par,'b','c') % no cell array

varargout = varargin; % simpler preallocation :)
for k = 1:numel(varargin)
    varargout{k} = S.(varargin{k});
end
end %end function

% function varargout = unpack(par,varnames)
% % See matlab answers:
% % 
% nvar = numel(varnames);
% varargout = cell(nvar,1);
% 
% for i=1:nvar
%     varargout{i} = par.(varnames{i});
% end
% 
% end
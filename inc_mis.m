function [z inc_z] = inc_mis(varargin)
%funzione per il calcolo delle incertezze
%
% X= x +- inc_x
% Y= y +- inc_y
% X,Y indipendenti (no covarianza)
%
%esempi
%  X=2 +- 0.2 Y=3 +- 0.3 Z=X^2 * Y^3 =108 +- 38.94
% [z inc_z]=inc_mis('molt_pot',[2 0.2 3 0.3 2 3])  
%
%  x=[1;2;3;4]  inc_x=[0.1;0.2;0.3;0.4] n=5   
%  Q=X^5  => z=[1;32;243;1024] inc_z=[0.5;16;121.5;512]
% [z inc_z]=inc_mis('pot',[x inc_x (n.*ones(size(x)))])
%  
%options(varargin)
%   'add' addizione di grandezze Q=X+Y argomento [x inc_x y inc_y]
%   'sott' sottrazione di grandezze Q=X-Y argomento [x inc_x y inc_y]
%   'dp' diretta proporzionalità Q=bX argomento [x inc_x b]
%   'pot' elevamento a potenza Q=X^n argomento [x inc_x n]
%   'molt' moltiplicazione di grandezze Q=XY argomento [x inc_x y inc_y]
%   'div' divisione di grandezze Q=X/Y argomento [x inc_x y inc_y]
%   'esp' esponenziale di grandezza Q=e^X argomento [x inc_x]
%   'log' logaritmo di grandezza Q=log(X) argomento [x inc_x]
%   'molt_pot' funzione della forma Q=X^(a)*Y^(b) argomento [x inc_x y inc_y a b]
%
%   x, y, inc_x, inc_y, b, n, a, b vettori colonna 
%
%Libardi Gabriele 15/03/2019

if length(varargin)>0
    for jj=1:length(varargin)
        if strcmp(varargin{jj},'add')
            mis_add = varargin{jj+1};
        end
        if strcmp(varargin{jj},'sott')
            mis_sott = varargin{jj+1};
        end
        if strcmp(varargin{jj},'dp')
            mis_dp = varargin{jj+1};
        end
        if strcmp(varargin{jj},'pot')
            mis_pot = varargin{jj+1};
        end
        if strcmp(varargin{jj},'molt')
            mis_molt = varargin{jj+1};
        end
        if strcmp(varargin{jj},'div')
            mis_div = varargin{jj+1};
        end
        if strcmp(varargin{jj},'esp')
            mis_esp = varargin{jj+1};
        end
        if strcmp(varargin{jj},'log')
            mis_log = varargin{jj+1};
        end
        if strcmp(varargin{jj},'molt_pot')
            mis_molt_pot = varargin{jj+1};
        end
    end
end

if exist('mis_add')
    x = mis_add(:,1);
    inc_x = mis_add(:,2);
    y = mis_add(:,3);
    inc_y = mis_add(:,4);
    
    z = x+y;
    inc_z = sqrt((inc_x).^2+(inc_y).^2);
end

if exist('mis_sott')
    x = mis_sott(:,1);
    inc_x = mis_sott(:,2);
    y = mis_sott(:,3);
    inc_y = mis_sott(:,4);
    
    z = x-y;
    inc_z = sqrt((inc_x).^2+(inc_y).^2);
end

if exist('mis_dp')
    x = mis_dp(:,1);
    inc_x = mis_dp(:,2);
    b = mis_dp(:,3);
    
    z = b.*x;
    inc_z = abs(b) .* inc_x;
end

if exist('mis_pot')
    x = mis_pot(:,1);
    inc_x = mis_pot(:,2);
    n = mis_pot(:,3);
    
    z = x.^n;
    inc_z = n .* abs(x.^(n-1)) .* inc_x;
end

if exist('mis_molt')
    x = mis_molt(:,1);
    inc_x = mis_molt(:,2);
    y = mis_molt(:,3);
    inc_y = mis_molt(:,4);
    
    z = x.*y;
    inc_z = sqrt((y.^2).*(inc_x.^2)+(x.^2).*(inc_y.^2));
end

if exist('mis_div')
    x = mis_div(:,1);
    inc_x = mis_div(:,2);
    y = mis_div(:,3);
    inc_y = mis_div(:,4);
    
    z = x./y;
    inc_z = sqrt((1./(y.^2)) .* (inc_x.^2) + ((x.^2)./(y.^4)) .* (inc_y.^2));
end

if exist('mis_esp')
    x = mis_esp(:,1);
    inc_x = mis_esp(:,2);
    
    z = exp(1).^x;
    inc_z = exp(1).^(x) .* inc_x;
end

if exist('mis_log')
    x = mis_log(:,1);
    inc_x = mis_log(:,2);
    
    z = log(x);
    inc_z = inc_x ./x;
end

if exist('mis_molt_pot')
    x = mis_molt_pot(:,1);
    inc_x = mis_molt_pot(:,2);
    y = mis_molt_pot(:,3);
    inc_y = mis_molt_pot(:,4);
    a = mis_molt_pot(:,5);
    b = mis_molt_pot(:,6);
    
    z = x.^(a) .* y.^(b);
    inc_z = sqrt( a.^2 .* x.^(2.*a -2) .* y.^(2.*b) .* (inc_x.^2) + b.^2 .* x.^(2.*a) .* y.^(2.*b -2) .* (inc_y.^2) );
end
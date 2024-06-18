function p = bvnl( dh, dk, r )
%BVNL
%  A function for computing bivariate normal probabilities.
%  bvnl calculates the probability that x < dh and y < dk. 
%    parameters  
%      dh 1st upper integration limit
%      dk 2nd upper integration limit
%      r   correlation coefficient
%  Example:
%    p = bvnl( 3, 1, .35 )
%

%   Author
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%   This function is based on the method described by 
%        Drezner, Z and G.O. Wesolowsky, (1989),
%        On the computation of the bivariate normal inegral,
%        Journal of Statist. Comput. Simul. 35, pp. 101-107,
%    with major modifications for double precision, for |r| close to 1,
%    and for matlab by Alan Genz - last modifications 7/98.
%
      p = bvnu( -dh, -dk, r );
%
%   end bvnl
%
function p = bvnu( dh, dk, r )
%BVNU
%  A function for computing bivariate normal probabilities.
%  bvnu calculates the probability that x > dh and y > dk. 
%    parameters  
%      dh 1st lower integration limit
%      dk 2nd lower integration limit
%      r   correlation coefficient
%  Example: p = bvnu( -3, -1, .35 )
%  Note: to compute the probability that x < dh and y < dk, 
%        use bvnu( -dh, -dk, r ). 
%

%   Author
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%
%    This function is based on the method described by 
%        Drezner, Z and G.O. Wesolowsky, (1989),
%        On the computation of the bivariate normal inegral,
%        Journal of Statist. Comput. Simul. 35, pp. 101-107,
%    with major modifications for double precision, for |r| close to 1,
%    and for Matlab by Alan Genz. Minor bug modifications 7/98, 2/10.
%
  if dh ==  inf | dk ==  inf, p = 0;
  elseif dh == -inf, if dk == -inf, p = 1; else p = phid(-dk); end
  elseif dk == -inf, p = phid(-dh);
  elseif r == 0, p = phid(-dh)*phid(-dk);  
  else, tp = 2*pi; h = dh; k = dk; hk = h*k; bvn = 0; 
    if abs(r) < 0.3      % Gauss Legendre points and weights, n =  6
      w(1:3) = [0.1713244923791705 0.3607615730481384 0.4679139345726904];
      x(1:3) = [0.9324695142031522 0.6612093864662647 0.2386191860831970];
    elseif abs(r) < 0.75 % Gauss Legendre points and weights, n = 12
      w(1:3) = [.04717533638651177 0.1069393259953183 0.1600783285433464];
      w(4:6) = [0.2031674267230659 0.2334925365383547 0.2491470458134029];
      x(1:3) = [0.9815606342467191 0.9041172563704750 0.7699026741943050];
      x(4:6) = [0.5873179542866171 0.3678314989981802 0.1252334085114692];
    else,                % Gauss Legendre points and weights, n = 20
      w(1:3) = [.01761400713915212 .04060142980038694 .06267204833410906];
      w(4:6) = [.08327674157670475 0.1019301198172404 0.1181945319615184];
      w(7:9) = [0.1316886384491766 0.1420961093183821 0.1491729864726037];
      w(10) =   0.1527533871307259;
      x(1:3) = [0.9931285991850949 0.9639719272779138 0.9122344282513259];
      x(4:6) = [0.8391169718222188 0.7463319064601508 0.6360536807265150];
      x(7:9) = [0.5108670019508271 0.3737060887154196 0.2277858511416451];
      x(10) =   0.07652652113349733;
    end, w = [w  w]; x = [1-x 1+x]; 
    if abs(r) < 0.925, hs = ( h*h + k*k )/2; asr = asin(r)/2;  
      sn = sin(asr*x); bvn = exp((sn*hk-hs)./(1-sn.^2))*w';
      bvn = bvn*asr/tp + phid(-h)*phid(-k);  
    else, if r < 0, k = -k; hk = -hk; end
      if abs(r) < 1, as = 1-r^2; a = sqrt(as); bs = (h-k)^2;
        asr = -( bs/as + hk )/2; c = (4-hk)/8 ; d = (12-hk)/80; 
        if asr > -100, bvn = a*exp(asr)*(1-c*(bs-as)*(1-d*bs)/3+c*d*as^2); end
        if hk  > -100, b = sqrt(bs); sp = sqrt(tp)*phid(-b/a);
          bvn = bvn - exp(-hk/2)*sp*b*( 1 - c*bs*(1-d*bs)/3 );
        end, a = a/2; xs = (a*x).^2; asr = -( bs./xs + hk )/2; 
        ix = find( asr > -100 ); xs = xs(ix); sp = ( 1 + c*xs.*(1+5*d*xs) ); 
        rs = sqrt(1-xs); ep = exp( -(hk/2)*xs./(1+rs).^2 )./rs; 
        bvn = ( a*( (exp(asr(ix)).*(sp-ep))*w(ix)' ) - bvn )/tp; 
      end 
      if r > 0, bvn =  bvn + phid( -max( h, k ) ); 
      elseif h >= k, bvn = -bvn;
      else, if h < 0, L = phid(k)-phid(h); else, L = phid(-h)-phid(-k); end
        bvn =  L - bvn;
      end
    end, p = max( 0, min( 1, bvn ) );
  end
%
%   end bvnu
%
function p = phid(z), p = erfc( -z/sqrt(2) )/2; % Normal cdf
%
% end phid
%
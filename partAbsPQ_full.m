function [Q,Pb,Pp,Pm] = partAbsPQ_full(ka,dc,L1,L2,refbounds,startunif)
% get P and Q values for a state on a partially absorbing interval
% particle can start uniformly within an interval,or can start at a node
% between two edges
% ------------
% input values:
% ------------
% refbounds = [false,false] for both edges absorbing, [true,false] for first edge
% reflecting, [false,true] for second edge reflecting
% L1, L2 = length of edges attached to the node
% ka = 2x1 or 1x2 vector containing absorption rate on edges connected to
% the node
% dc = diffusion coefficient
% all other parameters could be vectors (but must be all of the same size)
% -------------
% output values:
% -------------
% Pm = splitting probability to leave 1st boundary (next to absorbing
% interval)
% Pp = splitting probability to leave last boundary (far from absorbing
% interval)
% Pb = splitting probability to bind
% Q = average time to leave the state
% -------------


a1 = sqrt(ka(1)/dc);
a2 = sqrt(ka(2)/dc);

if(startunif)
	if(~any(refbounds))
		% both boundaries absorbing
		
		if(a1~=0)
			Pm = tanh(a1.*L1/2)./(a1.*L1);
			Pp = Pm;
			Q = a1.^(-3).*dc.^(-1).*L1.^(-1).*(a1.*L1+(-2).*tanh((1/2).*a1.*L1));
			Pb = 1-(2.*tanh(a1.*L1/2))./(a1.*L1);
		else 
			Pm = 1/2;
			Pp = 1/2;
			Q = L1^2/12./dc;
			Pb = 0;
		end
		
	elseif(refbounds(1) && ~refbounds(2))
		% ref at minus end, abs at plus end
						
		if(a1~=0)
			Pp = tanh(a1.*L1)./(a1.*L1);
			Pm = 0; % no exit from reflecting bound
			Q = a1.^(-3).*dc.^(-1).*L1.^(-1).*(a1.*L1+(-1).*tanh(a1.*L1));
			Pb = 1-tanh(a1.*L1)./(a1.*L1);
		else
			Pp = 1;
			Pm = 0;
			Q = L1.^2./3./dc;
			Pb = 0;
		end
		
	elseif(~refbounds(1) && refbounds(2))
		% abs at minus end, ref at plus end
				
		if(ka(1)~=0)
			Pp = 0; % no exit from reflecting bound
			Pm = tanh(a1.*L1)./(a1.*L1); 
			Q = a1.^(-3).*dc.^(-1).*L1.^(-1).*(a1.*L1+(-1).*tanh(a1.*L1));
			Pb = 1-tanh(a1.*L1)./(a1.*L1);
		else
			Pp = 0;
			Pm = 1;
			Q = L1.^2./3./dc;
			Pb = 0;
		end
		
	else
		% both boundaries reflecting
		Pp = 0;
		Pm = 0;
		Q = 1/ka;
		Pb = 1;
	end
	
else
	if(~any(refbounds))
		% both boundaries absorbing
		
		if(a1~=0 && a2~=0) %k1 and k2 both nonzero
			Pp = a2.*(a2.*cosh(a2.*L2)+a1.*coth(a1.*L1).*sinh(a2.*L2)).^(-1);
			Pm = a1.*(a1.*cosh(a1.*L1)+a2.*coth(a2.*L2).*sinh(a1.*L1)).^(-1);
			Pb = a1.*((-1)+cosh(a1.*L1)).*(a1.*cosh(a1.*L1)+a2.*coth(a2.*L2).*sinh( ...
							  a1.*L1)).^(-1)+a2.*((-1)+cosh(a2.*L2)).*sinh(a1.*L1).*(a2.*cosh( ...
							  a2.*L2).*sinh(a1.*L1)+a1.*cosh(a1.*L1).*sinh(a2.*L2)).^(-1);
			Q = (a1.^2.*a2.*dc.*coth(a1.*L1)+a1.*a2.^2.*dc.*coth(a2.*L2)).^(-1).*( ...
					  a2.*tanh((1/2).*a1.*L1)+a1.*tanh((1/2).*a2.*L2));			
		elseif(a1==0 && a2~=0) %k1 = 0, k2 nonzero
			Pp = a2.*L1.*(a2.*L1.*cosh(a2.*L2)+sinh(a2.*L2)).^(-1);
			Pm = (1+a2.*L1.*coth(a2.*L2)).^(-1);
			Pb = a2.*L1.*((-1)+cosh(a2.*L2)).*(a2.*L1.*cosh(a2.*L2)+sinh(a2.*L2)).^(-1);
			Q = (1/2).*a2.^(-1).*dc.^(-1).*L1.*(a2.*L1.*cosh(a2.*L2)+sinh(a2.*L2)) ...
					  .^(-1).*((-2)+2.*cosh(a2.*L2)+a2.*L1.*sinh(a2.*L2));
		elseif(a1~=0 && a2==0) %k2=0, k1 nonzero
			Pp = (1+a1.*L2.*coth(a1.*L1)).^(-1);
			Pm = a1.*L2.*(a1.*L2.*cosh(a1.*L1)+sinh(a1.*L1)).^(-1);
			Pb = a1.*L2.*((-1)+cosh(a1.*L1)).*(a1.*L2.*cosh(a1.*L1)+sinh(a1.*L1)).^(-1);
			Q = (2.*a1.*dc+2.*a1.^2.*dc.*L2.*coth(a1.*L1)).^(-1).*(a1.*L2.^2+2.* ...
					L2.*tanh((1/2).*a1.*L1));
		elseif(a1==0 && a2==0) % no absorption
			Pp = 1 - L2/(L1+L2);
			Pm = 1 - L1/(L1+L2);
			Pb = 0;
			Q = L1.*L2./2./dc;
		end
		
	elseif(refbounds(1) && ~refbounds(2))
		% ref at minus end, abs at plus end
		Pm = 0;
		if(a1~=0 && a2~=0) %k1 and k2 both nonzero
			Pp = a2.*(a2.*cosh(a2.*L2)+a1.*sinh(a2.*L2).*tanh(a1.*L1)).^(-1);
			Pb = a1.*(a1+a2.*coth(a1.*L1).*coth(a2.*L2)).^(-1)+a2.*cosh(a1.*L1).*(( ...
					-1)+cosh(a2.*L2)).*(a2.*cosh(a1.*L1).*cosh(a2.*L2)+a1.*sinh(a1.*L1).*sinh(a2.*L2)).^(-1);
			Q = a1.^(-1).*a2.^(-1).*dc.^(-1).*(a1+a2.*coth(a1.*L1).*coth(a2.*L2)) ...
				  .^(-1).*(a2.*cosh(a2.*L2).*coth(a1.*L1)+a1.*sinh(a2.*L2)).*(a2.* ...
				  cosh(a1.*L1).*cosh(a2.*L2)+a1.*sinh(a1.*L1).*sinh(a2.*L2)).^(-1).* ...
				  (a2.*sinh(a1.*L1)+a1.*cosh(a1.*L1).*tanh((1/2).*a2.*L2));
		elseif(a1==0 && a2~=0) %k1 = 0, k2 nonzero
			Pp = sech(a2.*L2);
			Pb = 1+(-1).*sech(a2.*L2);
			Q = a2.^(-2).*dc.^(-1).*(1+(-1).*sech(a2.*L2)+a2.*L1.*tanh(a2.*L2));
		elseif(a1~=0 && a2==0) %k2=0, k1 nonzero
			Pp = (1+a1.*L2.*tanh(a1.*L1)).^(-1);
			Pb = a1.*L2.*(a1.*L2+coth(a1.*L1)).^(-1);
			Q = (a1.*L2.^2.*cosh(a1.*L1)+2.*L2.*sinh(a1.*L1)).*(2.*a1.*dc.*cosh( ...
						a1.*L1)+2.*a1.^2.*dc.*L2.*sinh(a1.*L1)).^(-1);
		elseif(a1==0 && a2==0) % no absorption
			Pp = 1;
			Pb = 0;
			Q = (1/2).*dc.^(-1).*L2.*(2.*L1+L2);
		end

	elseif(~refbounds(1) && refbounds(2))
		error('Function not set up for this boundary condition yet');
		% abs at minus end, ref at plus end

	else
		% both boundaries reflecting
		Pm = 0;
		Pp = 0;
		if(a1~=0 && a2~=0)
			Pb = 1;
			Q = (a2+a1.*coth(a1.*L1).*tanh(a2.*L2)).*(a1.^2.*a2.*dc+a1.*a2.^2.*...
						dc.*coth(a1.*L1).*tanh(a2.*L2)).^(-1);
		elseif(a1==0 && a2~=0)
			Pb =  a2.*(a2+a1.*coth(a2.*L2).*tanh(a1.*L1)).^(-1);
			Q = a2.^(-2).*dc.^(-1).*(1+a2.*L1.*coth(a2.*L2));
		elseif(a1~=0 && a2==0)
			Pb = a1.*(a1+a2.*coth(a1.*L1).*tanh(a2.*L2)).^(-1);
			Q = a1.^(-2).*dc.^(-1).*(1+a1.*L2.*coth(a1.*L1));
		elseif(a1==0 && a2==0)
			Pb = 0;
			Q = inf;
		end
	end		
end



end

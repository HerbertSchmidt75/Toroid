package IMOC25
  package Force
  end Force;

  package Leakage
    model OuterHalfHollowTorus
      extends BaseClasses.Force;
      SI.Length l=s;
      parameter SI.Length w=0.1;
      parameter SI.Radius r=0.01;
    
      Real eta, root, arccot, ln;
      SI.Length g=l;
      SI.Radius ri=g/2, ro=t+ri;
      parameter SI.Radius R=w/(2*pi);
      parameter SI.Length t=r;
      parameter SI.Permeance Gm0=pi*mu_0*t;
    equation
      eta=R/t*ln(ro/ri);
      if eta<1 then
        root = sqrt(1-eta^2);
        ln = log((1+root)/(1-root));
        G_m = Gm0*2*root/ln;
        dGmBydx = -Gm0*R/(ri*ro)*(2/eta-eta/root*ln)/ln^2*dlBydx;
      elseif eta>1 then
        root = sqrt(eta^2-1);
        arccot = atan(1/root);
        G_m = Gm0*root/(pi/2-arccot);
        dGmBydx = -Gm0*R/(ri*ro)*eta/(2*root)*(pi/2-arccot-root/eta^2)/(pi/2- arccot)^2*dlBydx;
      else
        G_m = Gm0;
        dGmBydx = -Gm0*R/(3*ri*ro)*dlBydx; 
      end if;
    end OuterHalfHollowTorus;
  end Leakage;
end IMOC25;

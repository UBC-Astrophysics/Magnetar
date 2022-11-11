#include <math.h>

void emissivity(    
            double *thetak,
		    double *phik,
		    double *ene,
		    int N,
		    double mag_strength,
		    double thetab,
		    double dens,
		    double Z,
		    double A,
		    int fixedion,
		    double *emisX,
		    double *emisO)  {
	
  // quantities that are independent of the photon angles and energy (ene)
  
  double B13=mag_strength/1e13;
  double ece=115.77*B13;               // electron cyclotron energy (Eq. 1)
  double epe=0.0288*sqrt(dens*Z/A);    // plasma frequency (Eq. 3)
  double eci=0.0635*(Z/A)*B13;         // ion cyclotron energy (Eq. 8)
  double ecfx=epe*epe/ece;	       // max energy of L wave 
  double ec=eci+ecfx;                  // characteristic energy (Eq. 9) : max L-wave including ion cyclotron
  double n0=sqrt(1+epe*epe/2/ece/eci); // index of refraction at eci (Eq. 10) after Eq. 19
  double twon0overoneplusn0sqr=2*n0/(1+n0)/(1+n0);
  double stb=sin(thetab);
  double ctb=cos(thetab);
  double fourthrootstb=sqrt(sqrt(stb));
  double p=0.1*(1+stb)/(1+B13);
  double a1factor=1/(1+0.6*B13*ctb*ctb);
  double j1ectilde=0.5+0.05/(1+B13)+stb*0.25;			  // Eq. 26
  double j0=4/(sqrt(ec/eci)+1)/(sqrt(eci/ec)+1);                  // after Eq. 16 

    
  for (int i=0;i<N;i++) {
    double stk, ctk, sphik, cphik;
    double calphai, calphar, dum;
    double alphai, alphar, alpha, salphar;
    double salpha, calpha;
    double epetilde,ecfxtilde,ectilde;
    double ephoton=ene[i];
    double ypre1r,xpre1r,ypre2r,xpre2r;
    double emis1, emis2, emisXloc, emisOloc;
    double ja, jb, jeci;
    double a1;
    double jectilde, logectildeoeci;
    double eratio;
	
    stk=sin(a1=thetak[i]);
    ctk=cos(a1);
    sphik=sin(a1=phik[i]);
    cphik=cos(a1);
      
//    sincos(thetak[i],&stk,&ctk);
//    sincos(phik[i],&sphik,&cphik);

    ctk=(ctk<0 ? 0 : (ctk>0.999 ? 0.999 : ctk));
    stk=(stk<0.01414178207 ? 0.01414178207 : (stk>1 ? 1 : stk));

    dum=ctb*ctk;
    calphai=(stb*stk*cphik);
    calphar=calphai+dum;			// cosine of angle between relfected photon k and field Eq. 5
    calphai-=dum;                               // cosine of angle between incident photon k and field  Eq. 5

    calphar = (calphar < -0.9999 ? -0.9999 : (calphar>0.9999 ? 0.9999: calphar));
    calphai = (calphai < -0.9999 ? -0.9999 : (calphai>0.9999 ? 0.9999: calphai));
    
    alphai=acos(calphai);    
    alphar=acos(calphar);
    salphar=sqrt(1-calphar*calphar);

    epetilde=epe*sqrt(3.-2.*ctk);	       // Eq. 14
    ecfxtilde=epetilde*epetilde/ece; 
    ectilde=eci+ecfxtilde;                     // Eq. 12 -- max L wave including ion cyclotron
    logectildeoeci=log(ectilde/eci);	  
	  
    if (alphai<alphar) {                       // minimum of the two alphas Eq. 16
      alpha=alphai;
      salpha=sqrt(1-calphai*calphai);
      calpha=calphai;
    } else {
      alpha=alphar;
      salpha=salphar;
      calpha=calphar;
    }
	
    a1=(1-ctb*ctb*cphik-stb*stb*calpha);                           // Eq. 25
	  
    if (ephoton<eci) {
      double nplus,nminus,Rplus,Rminus;
      nplus=sqrt(1+epe*epe/(ece*(ephoton+eci)));			// Eq. 10
      nminus=sqrt(1-epe*epe/(ece*(ephoton-eci)));
      Rplus=(nplus-1)/(nplus+1);
      Rplus*=Rplus;
      Rminus=(nminus-1)/(nminus+1);
      Rminus*=Rminus;
#define AE ((1-ctb)/2/sqrt(1+B13)+stk*stk*stk*stk*(0.7-0.45/j0)*(1-calpha)) // Eq. 16
      ja=(1-AE)*(1-0.5*(Rminus+Rplus));                            // Eq. 15
    } else {
      ja=0;
    }
    jectilde=0.5+0.05/(1+B13)*(1+ctb*sphik)-0.15*(1-ctb)*salpha;   // Eq. 18
    eratio=ephoton/ectilde;
    if (fixedion) {
      jb=jectilde/(1-p+p*pow(eratio,-0.6));             // Eq. 22
    } else { 
      jeci=twon0overoneplusn0sqr*(1+(ctb-cphik)/(2*(1+B13)));		       // Eq. 19
      jb=pow(eratio,log(jectilde/jeci)/logectildeoeci)*jectilde;  // Eq. 17
    } 
    
    if ((ephoton<eci) && (ja>jb)) {
      emis1=(1-a1*a1factor)*ja;                                              // Eq. 25
      emis2=2*ja-emis1;
    } else {
      double eL=(epe*(1+1.2*(1-ctk)*sqrt(1-ctk))*(1-0.3333333333*stb*stb));
      double wL=0.8*pow(ectilde/epe,0.2)*sqrt(sin(alpha*0.5))*(1+stb*stb);
      double bigX=(ephoton-eL)/(2*epe*wL*(1-ctk));
      double dd=ephoton/epe;
      double l=(0.17*epe/(ec*(1+bigX*bigX*bigX*bigX))+0.21*exp(-dd*dd))*stk*stk*wL;    
      double rl=fourthrootstb*(2-salpha*salpha*salpha*salpha)*(l/(1+l));
		   
      double ntilde=1-epetilde*epetilde/(ece*(ephoton-eci));
      double jc=1+ntilde;
      double jb1;

      jc=4*ntilde/jc/jc;
	
      if (fixedion) {
	jb1=j1ectilde/(0.1+0.9*pow(eratio,-0.4));                      // Eq. 27		
      } else {
	double j1eci;
	j1eci=(1-a1)*jeci;                                              // Eq. 26   
	jb1=pow(eratio,log(j1ectilde/j1eci)/logectildeoeci)*j1ectilde;  // Eq. 26
      }

      emis1=jb1*(1-jc)+jc*(1-rl);
      emis2=2*(jb*(1-jc)+jc/(1+l))-emis1;
    }
  
  // transformations from mode 1 and 2 to X and O (Eq. B.9 - B.12)	  
  ypre1r = (ctb*stk-stb*ctk*cphik)/salphar;
  ypre1r *=  ypre1r;
  xpre1r = stb*sphik/salphar;
  xpre1r *= xpre1r;
  ypre2r = xpre1r;
  xpre2r = ypre1r;

  emisXloc=xpre1r*emis2+xpre2r*emis1;
  emisOloc=ypre1r*emis2+ypre2r*emis1;
  emisX[i]=(emisXloc<0 ? 0 : (emisXloc>1 ? 1 : emisXloc));
  emisO[i]=(emisOloc<0 ? 0 : (emisOloc>1 ? 1 : emisOloc));
       
  }
}

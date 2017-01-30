#include <iostream>
#include <fstream>

#include "particles.h"

#include "TMath.h"

using namespace std;
using namespace TMath;

TPseudoParticle::TPseudoParticle(
                                 const PARTICLETYPE& hd,
                                 const double& A_,
                                 const double& rf_,
                                 const double& phif_,
                                 const double& thetaf_,
                                 const double& qf_,
                                 const double& tf_,
                                 const double& af_,
                                 const double& alphaf_,
                                 const double& charge_,
                                 double* bparam) :
Anumber(A_),
rf(rf_),
phif(phif_),
thetaf(thetaf_),
qf(qf_),
tf(tf_),
af(af_),
alphaf(alphaf_),
iE(qf_),
charge(charge_)
{
    if (hd == NUCLEUS) M = 0.938;
    else if (hd == LEPTON) M = 0.000511;
    
    //rn = TRandom3(0);
    bf = Bfield(bparam[0],bparam[1],bparam[2],bparam[3],bparam[4], bparam[5], bparam[6], bparam[7], bparam[8], bparam[9], bparam[10], bparam[11]);
    dW = vector<double>(4,0);
    vd = vector<double>(3,0);
    Ktensor = vector<double>(4,0);
    Kderiv = vector<double>(5,0);
    
}


void TPseudoParticle::Evolve(const TEnvironment* env, TRandom3* rn) {
    
    //TFile rootfile("p0.root","RECREATE");
    // TNtupleD nt("nt","Particle path", "rf:thetaf:phif:qf:af:alphaf:tf");
    
    const double Rmax = env->GetRmax();
    const double dt0 = env->Getdt();
    //  const double sqrtdt = sqrt(dt);
  
    //cout << endl << endl;   

    //ofstream outfile("log.dat", ios::out);
 
    while (rf < Rmax && rf > SunRadius) {
        
	//cout << " --- rf = " << rf << "; Rmax = " << Rmax << endl;       

        bf.Set(rf, thetaf, phif, qf, beta_velocity(qf), momentum(qf), charge);
        bf.GetVd(vd);
        bf.GetKTensors(Ktensor, Kderiv);
        double Vsw = bf.GetVsw();
        
        double dt = dt0*pow(rf,2)/Ktensor[0];
        double sqrtdt = sqrt(dt);
        
        tf += dt;
        
        for (int j = 0; j < 2; j++) {
            double R1 = sqrt(-2.0*log(rn->Rndm()));
            double R2 = TMath::TwoPi()*(rn->Rndm());
            dW[2*j] = R1*cos(R2)*sqrtdt;
            dW[2*j+1] = R1*sin(R2)*sqrtdt;
        }
        
        
        double Ar = pow(rf,-2)*Kderiv[0] + 1.0/(rf*sin(thetaf))*Kderiv[1] - Vsw - vd[0];
        
        double Aphi = pow(rf*sin(thetaf),-2)*Kderiv[3] + 1.0/(pow(rf,2)*sin(thetaf))*Kderiv[4] - vd[1]/(rf*sin(thetaf));
        
        double Atheta = 1.0/(pow(rf,2)*sin(thetaf))*Kderiv[2] - vd[2]/rf;
        
        double Gamma = (qf+2.0*M)/(qf+M);
        double AE = pow(rf,-2)/3.0*2.0*(rf*Vsw)*Gamma*qf;
        
        double secterm = sqrt(2.0/Ktensor[3])*Ktensor[1];
        double Wt = sqrt(2.0*Ktensor[2])/rf;
        double Wphi = sqrt(2.0*Ktensor[3])/(rf*sin(thetaf));

	/*
	cout << " pow(rf,-2)*Kderiv[0]*dt = " << pow(rf,-2)*Kderiv[0]*dt << endl;
	cout << " 1.0/(rf*sin(thetaf))*Kderiv[1]*dt = " << 1.0/(rf*sin(thetaf))*Kderiv[1]*dt << endl;
	cout << " Vsw*dt = " << Vsw*dt << endl;
	cout << " vd*dt (drift term, charge dependent)  = " << sqrt(vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2])*dt << endl;
	cout << " sqrt(2.0*Ktensor[0]-pow(secterm,2))*dW[0] (stochastic) " <<   sqrt(2.0*Ktensor[0]-pow(secterm,2))*dW[0] << endl;
	cout << " secterm*dW[1] (stochastic) " << secterm*dW[1] << endl << endl;     

	outfile << rf << "\t" << pow(rf,-2)*Kderiv[0]*dt << "\t" << 1.0/(rf*sin(thetaf))*Kderiv[1]*dt << "\t" << Vsw*dt << "\t" << sqrt(vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2])*dt << "\t" << sqrt(2.0*Ktensor[0]-pow(secterm,2))*dW[0] << " " << secterm*dW[1] << endl;
	*/

        rf += Ar*dt + sqrt(2.0*Ktensor[0]-pow(secterm,2))*dW[0] + secterm*dW[1];
        
        thetaf += (Atheta*dt + Wt*dW[2]);
        
        phif += (Aphi*dt + Wphi*dW[1]);
        
        qf += AE*dt;
        //cout << AE*dt << endl;
        
        // Consistency conditions
        
        thetaf = fmod(thetaf,TMath::TwoPi());
        
        if (thetaf < 0) {
            thetaf *= -1.0;
            phif += TMath::Pi();
        }
        // while (thetaf > TMath::TwoPi()) thetaf -= TMath::TwoPi();
        
        
        if (thetaf > TMath::Pi()) {
            thetaf = TMath::TwoPi() - thetaf;
            phif += TMath::Pi();
        }
        
        if (rf<0) {
            phif += TMath::Pi();
            thetaf = TMath::Pi()-thetaf;
            rf *= -1.0;
        }
        
        phif = fmod(phif,TMath::TwoPi());
        if (phif < 0) phif += TMath::TwoPi();
        //if (phif > TMath::TwoPi()) phif -= TMath::TwoPi();
        //while (phif > TMath::TwoPi()) phif -= TMath::TwoPi();
        
        
        // qf does not evolve
        //        LF -= dt*(loss->Eval(yf) + k0*(yf/Ymax)/Ymax);
        
        /*        if ( fabs(tf-tsource) < 0.05 && fabs(xf-xs) < ds && fabs(yf-ys) < ds && fabs(zf-zs) < ds) {
         alphaf = exp(LF);
         af += alphaf*source->Eval(xf,yf,zf);
         }
         */
    }
    // alphaf = exp(LF);
    af += boundarycondition(qf,tf,rf);
    
    //outfile.close();	
    //nt.Write();
    //rootfile.Close();
    
    return ;
    
}



# FourTopProject

//Single Lepton Channel Code

#define Events_cxx
#include "Events.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
void Events::Loop()
{
if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes=0,nb=0;
TH1F*MyHist = new TH1F("MyHist","Events vs H_t",10,2,12);
Long64_t n = 101; 
float ne=0; // event counter
Int_t counter = 0;
float CSVV2[960869] = {}; // these arrays will contain the b-jet discrim value for each jet, size is the number of jets in 
                           // data set
for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	Long64_t elec =0;
	Long64_t muon = 0;
	Long64_t u=0;
	float SumTrans=0; // scalar sum of transverse momentum counter
	Long64_t p=0;
	for(Long64_t j=0;j<nJet;++j){         // this small loop calculates scalar sum of transverse momentum
	SumTrans = SumTrans+Jet_pt[j];
	}
	float x_t=0.988;   // these are the tight, medium and loose working points. i denotes the users intermediate value.
	float x_m=0.93;     // these are the working points for the CSVV2 algorithm
	float x_l=0.67;     // THe b denotes that these are the working points for the CSVV2 algorithm

	float x_tb = 0.925;
	float x_mb = 0.69;
	float x_lb = 0.23;
	float pt_misselecx=0; // These missing transverse momentum variables of electron etc... are used to calculate the 
	float pt_misselecy=0; // sum of missing transverse momentum later
        float pt_missmuonx=0;	
	float pt_missmuony=0;
	float pt_missx=0;
	float pt_missy =0;
	float pt_miss=0;
	float pt_misstotx=0;
	float pt_misstoty=0;
	float lep_pt =0;  // This is used to calculate the sum of lepton momentum  in an event
	Long64_t repeat =0;
	float Discrim[40]={};
	float First_Discrim=0; // This is the discriminatory value of the b-jet with the highest discriminatory value in an event
	float Second_Discrim=0; // second highest etc...
	float Third_Discrim = 0;
	Long64_t T=0; // counters for how many b-jets in an event have a discrim val greater than T etc...
	Long64_t M=0;
	Long64_t L=0;
	Long64_t mine=0;
	Long64_t Di_elec =-10;   
	float second_elec=-10;
	float second_muon=-10;
	Long64_t Electron_number =0;
	Long64_t Muon_number =0;
for (Long64_t Jet =0;Jet<nJet;++Jet) {     // This loop finds out if electrons in an event obey the isolation condition 
		if (Jet_nElectrons[Jet] > 0) {          // and momentum condition 
		float I_rel = 0;
		float pT = 0;
		for (Long64_t i=0;i<Jet_nElectrons[Jet];++i) {
			for (Long64_t NumPart = 0;NumPart<nGenPart;++NumPart){
				float deltaphi = Jet_phi[Jet]-GenPart_phi[NumPart];
				if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
				float deltaeta = Jet_eta[Jet]-GenPart_eta[NumPart];
				float deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
				if (deltaR<0.4){
				deltaphi = Electron_phi[Electron_number]-GenPart_phi[NumPart];
				if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
				deltaeta = Electron_eta[Electron_number]-GenPart_eta[NumPart];
				deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
				if (deltaR<0.4){
				pT += GenPart_pt[NumPart];
				}
				}
			}
		I_rel = (pT)/Electron_pt[Electron_number];
		if (I_rel < 0.15 && Electron_pt[Electron_number] >35 && abs(Electron_eta[Electron_number])<2.1){
		if (Di_elec ==-10) Di_elec = Electron_number;
		else second_elec = Electron_number;
		lep_pt+=Electron_pt[Electron_number];
		elec=elec+1;
		pt_misselecx = Electron_pt[Electron_number]*cos(Electron_phi[Electron_number]);
		pt_misselecy = Electron_pt[Electron_number]*sin(Electron_phi[Electron_number]);
}
Electron_number += 1;
}
}
}

	Long64_t Di_muon =-10;
for (Long64_t Jet =0;Jet<nJet;++Jet) { // this loop finds out if muons in an events obey isolation and momentum conditions
	if (Jet_nMuons[Jet] > 0) {
	float I_rel = 0;
	float pT = 0;
	for (Long64_t i=0;i<Jet_nMuons[Jet];++i) {
		for (Long64_t NumPart = 0;NumPart<nGenPart;++NumPart){
			float deltaphi = Jet_phi[Jet]-GenPart_phi[NumPart];
			if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
			float deltaeta = Jet_eta[Jet]-GenPart_eta[NumPart];
			float deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
			if (deltaR<0.4){
			deltaphi = Muon_phi[Muon_number]-GenPart_phi[NumPart];
			if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
			deltaeta = Muon_eta[Muon_number]-GenPart_eta[NumPart];
			deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
			if( deltaR<0.4){
			pT += GenPart_pt[NumPart];
			}
			}
		}
	I_rel = (pT)/Muon_pt[Muon_number];
	if (I_rel < 0.15 && Muon_pt[Muon_number] >26 && abs(Muon_eta[Muon_number])<2.1){
	if (Di_muon ==-10) Di_muon = Muon_number;
	else second_muon = Muon_number;
	lep_pt=Muon_pt[Muon_number];
	muon=muon+1;
	pt_missmuonx = Muon_pt[Muon_number]*cos(Muon_phi[Muon_number]);
	pt_missmuony = Muon_pt[Muon_number]*sin(Muon_phi[Muon_number]);
	}
	Muon_number += 1;
	}
	}
}

		if (elec+muon==1){  // this requires exactly one electron or muon therefore in single lepton channel
		Long64_t jetnumber =0;  // This condition can be trivially changed if one only wants to look at single electron or single
		if (elec==1) jetnumber =8;   // muon channel.
		if (muon==1) jetnumber =7;
		for (Long64_t l=0;l<nJet;++l){
		pt_misstotx += Jet_pt[l]*cos(Jet_phi[l]); 
		pt_misstoty += Jet_pt[l]*sin(Jet_phi[l]);
		}
		pt_missx = pow(pt_misstotx +pt_misselecx+pt_missmuonx,2);
		pt_missy = pow(pt_misstoty +pt_misselecy+pt_missmuony,2);      
		pt_miss = sqrt(pt_missx+pt_missy);  // final caluation of missing transverse momentum 
		if (nJet>= jetnumber){
		for (Long64_t NumJet= 0;NumJet<nJet;++NumJet) {  // this loop stores discriminatory value of a jet and counts how many
                        counter = counter + 1;    // jets pass the jet requirements for your baseline selection
                        CSVV2[counter] = Jet_btagCSVV2[NumJet];  
                        if (abs(Jet_eta[NumJet])<2.5 && Jet_pt[NumJet] > 30){
                        p+=1;
                        Discrim[NumJet] = CSVV2[counter];
                        }  
                }
                for(Long64_t i = 0;i < 40; ++i){   // this loop counts how many b-jets pass tight, medium and loose cuts etc..
                        if (Discrim[i] > x_l) L+=1;
                        if (Discrim[i] > x_m) M+=1;
                        if (Discrim[i] > x_t) T+=1;
			if (Discrim[i] >x_i) mine +=1;
                }
                for(Long64_t i = 0;i < 40; ++i){    // this loop stores which jet has the highest, second and third highest
                        if(Discrim[39] < Discrim[i]) Discrim[39] = Discrim[i]; // discrim val. It can be extended further
                }
                repeat =0;
                for(Long64_t i = 0;i < 40; ++i){
                        if(Discrim[i]== Discrim[39]) repeat += 1;
                        if(repeat==3) Discrim[38]=Discrim[39];
                        if(Discrim[38] < Discrim[i] && Discrim[i]<Discrim[39]) Discrim[38] = Discrim[i];
                }
                repeat =0;
                if (Discrim[39]==Discrim[38]) repeat =-2;
                for(Long64_t i = 0;i < 40; ++i){
                        if( Discrim[i]==Discrim[38]) repeat +=1;
                        if( repeat ==3) Discrim[37]=Discrim[38];
                        if(Discrim[37] < Discrim[i] && Discrim[i]<Discrim[38]) Discrim[37] = Discrim[i];
                }
                First_Discrim = Discrim[39]; 
                Second_Discrim = Discrim[38];
                Third_Discrim = Discrim[37]; 
                }


else {
for (Long64_t NumJet= 0;NumJet<nJet;++NumJet) {
counter=counter+1;
CSVV2[counter]=-10;
CMVA[counter]=-10;
}}
\\ the below condition bascially takes all the variables that have been calculated for this event and checks if they pass
\\ the baseline selection 
if ((muon ==1 && T>=1 && L >=2 && M>=2 && SumTrans>500 && p>=7 && pt_miss>50) or (elec ==1 && T>=1 && M>=2 && L>=2 && SumTrans>500 && p>=8 && pt_miss>50) ){
 ne=ne+1;
 MyHist->Fill(M);} \\ fills a histogram with whatever variable that has been calculated that you want to look at 
}
if (ientry < 0) break;
nb=fChain->GetEntry(jentry);nbytes+=nb;
   
}
TCanvas *c1 = new TCanvas("c1");
c1->SetLogy(); 
MyHist->Draw(); // draws the histogram

TFile myGraphFile("yourhist.root","RECREATE"); \\ saves the histogram to your directory
myGraphFile.cd();
MyHist->Write();
myGraphFile.Close();

}

Opposite Sign Dilepton Channel Code

#define Events_cxx
#include "Events.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
void Events::Loop()
{
if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes=0,nb=0;
TH1F*MyHist = new TH1F("MyHist","Events vs H_t",10,2,12);
Long64_t n = 101; 
float ne=0;
Int_t counter = 0;
float CSVV2[960869] = {};
float CMVA[960869]={};	
float H_t[65555]={};
        Long64_t electroncount=0;
Long64_t muoncount=0;
for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	Long64_t elec =0;
	Long64_t muon = 0;
	Long64_t u=0;
	float SumTrans=0;
	Long64_t p=0;
	for(Long64_t j=0;j<nJet;++j){
	SumTrans = SumTrans+Jet_pt[j];
	}
	H_t[jentry]=SumTrans;
	float Discrim[40]={};
	float x_t=0.988;
	float x_m=0.93;
	float x_l=0.67;
	float x_i=0.855;
	float x_ib = 0.53;
	float x_tb =0.925;
	float x_mb = 0.69;
	float x_lb =0.23;
	Long64_t mine =0;
	Long64_t repeat =0;
	float pt_misselecx=0;
	float pt_misselecy=0;
        float pt_missmuonx=0;	
	float pt_missmuony=0;
	float pt_missx=0;
	float pt_missy =0;
	float pt_miss=0;
	float pt_misstotx=0;
	float pt_misstoty=0;
	float lep_pt =0;
	float Discrim1 =0;
	float First_Discrim=0;
	float Second_Discrim=0;
	float Third_Discrim = 0;
	Long64_t T=0;
	Long64_t M=0;
	Long64_t L=0;
        Long64_t  eleccharge[8]={};
        float elecmass[8]={};
        float  elecpt[8]={};
	float elecphi[8]={};
	float eleceta[8]={};
	float Di_elec =-10;
	float second_elec=-10;
	float second_muon=-10;
	Long64_t Electron_number =0;
	Long64_t Muon_number =0;
for (Long64_t Jet =0;Jet<nJet;++Jet) {
		if (Jet_nElectrons[Jet] > 0) {
		float I_rel = 0;
		float pT = 0;
		for (Long64_t i=0;i<Jet_nElectrons[Jet];++i) {
			for (Long64_t NumPart = 0;NumPart<nGenPart;++NumPart){
				float deltaphi = Jet_phi[Jet]-GenPart_phi[NumPart];
				if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
				float deltaeta = Jet_eta[Jet]-GenPart_eta[NumPart];
				float deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
				if (deltaR<0.4){
				deltaphi = Electron_phi[Electron_number]-GenPart_phi[NumPart];
				if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
				deltaeta = Electron_eta[Electron_number]-GenPart_eta[NumPart];
				deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
				if (deltaR<0.4){
				pT += GenPart_pt[NumPart];
				}
				}
			}
		I_rel = (pT)/Electron_pt[Electron_number];
		if (I_rel < 0.15 && Electron_pt[Electron_number] >20 && abs(Electron_eta[Electron_number])<2.4){
		if (Di_elec ==-10) Di_elec = Electron_number;
		else second_elec = Electron_number;
		lep_pt += Electron_pt[Electron_number];
		eleccharge[Electron_number] =Electron_charge[Electron_number];
		elecpt[Electron_number] = Electron_pt[Electron_number];
		eleceta[Electron_number]= Electron_eta[Electron_number];
		elecphi[Electron_number]= Electron_phi[Electron_number];
		elec=elec+1;
		pt_misselecx = Electron_pt[Electron_number]*cos(Electron_phi[Electron_number]);
		pt_misselecy = Electron_pt[Electron_number]*sin(Electron_phi[Electron_number]);
}
Electron_number += 1;
}
}
}

	Long64_t muoncharge[8]={};
	float muonmass[8]={};
	float muonpt[8]={};
	float muonphi[8]={};
	float muoneta[8]={};
	float Di_muon =-10;
for (Long64_t Jet =0;Jet<nJet;++Jet) {
	if (Jet_nMuons[Jet] > 0) {
	float I_rel = 0;
	float pT = 0;
	for (Long64_t i=0;i<Jet_nMuons[Jet];++i) {
		for (Long64_t NumPart = 0;NumPart<nGenPart;++NumPart){
			float deltaphi = Jet_phi[Jet]-GenPart_phi[NumPart];
			if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
			float deltaeta = Jet_eta[Jet]-GenPart_eta[NumPart];
			float deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
			if (deltaR<0.4){
			deltaphi = Muon_phi[Muon_number]-GenPart_phi[NumPart];
			if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
			deltaeta = Muon_eta[Muon_number]-GenPart_eta[NumPart];
			deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
			if (deltaR<0.4){
			pT += GenPart_pt[NumPart];
			}
			}
		}
	I_rel = (pT)/Muon_pt[Muon_number];
	if (I_rel < 0.15 && Muon_pt[Muon_number] >20 && abs(Muon_eta[Muon_number])<2.4){
	if (Di_muon ==-10) Di_muon = Muon_number;
	else second_muon = Muon_number;
	lep_pt += Muon_pt[Muon_number];
	muoncharge[Muon_number] =Muon_charge[Muon_number];
	muonpt[Muon_number] = Muon_pt[Muon_number];
	muoneta[Muon_number]=Muon_eta[Muon_number];
	muonphi[Muon_number]=Muon_phi[Muon_number];
	muon=muon+1;
	pt_missmuonx = Muon_pt[Muon_number]*cos(Muon_phi[Muon_number]);
	pt_missmuony = Muon_pt[Muon_number]*sin(Muon_phi[Muon_number]);
	}
	Muon_number += 1;
	}
	}
}
		Long64_t leadposmuon =0;  \\ this section basically catergorises every electron and muon that have pt>7 for electron and
		Long64_t subposmuon =0;   \\ pt>5 for muon and I_rel>0.15 which is the standard isolation condition 
		Long64_t leadnegmuon=0;   \\ it catergorises them by their charge and pt 
		Long64_t subnegmuon=0;
		Long64_t leadposelec =0;
		Long64_t subposelec =0;
		Long64_t leadnegelec=0;
		Long64_t subnegelec =0;
		float lpm_pt=0; float lpm_eta=0; float lpm_phi=0; float lpm_pt1=0; float lpm_eta1=0; float lpm_phi1=0;
		float spm_pt=0; float spm_eta=0; float spm_phi=0; float spm_pt1=0; float spm_eta1=0; float spm_phi1=0;
		float lnm_pt=0; float lnm_eta=0; float lnm_phi=0; float lnm_pt1=0; float lnm_eta1=0; float lnm_phi1=0;
		float snm_pt=0; float snm_eta=0; float snm_phi=0; float snm_pt1=0; float snm_eta1=0; float snm_phi1=0;		
		float lpe_pt=0; float lpe_eta=0; float lpe_phi=0; float lpe_pt1=0; float lpe_eta1=0; float lpe_phi1=0;
		float lne_pt=0; float lne_eta=0; float lne_phi=0; float lne_pt1=0; float lne_eta1=0; float lne_phi1=0;
		float spe_pt=0; float spe_eta=0; float spe_phi=0; float spe_pt1=0; float spe_eta1=0; float spe_phi1=0;
		float sne_pt=0; float sne_eta=0; float sne_phi=0; float sne_pt1=0; float sne_eta1=0; float sne_phi1=0;
		for(Long64_t k=0;k<8;++k){
			if(muoncharge[k]==1 && muonpt[k]>25){
			leadposmuon +=1;
			if(leadposmuon ==1){
			lpm_pt = muonpt[k];
			lpm_eta = muoneta[k];
			lpm_phi =muonphi[k];
			}
			if(leadposmuon==2){
			lpm_pt1 = muonpt[k];
			lpm_eta1 = muoneta[k];
			lpm_phi1 = muonphi[k];
			}}
			if(muoncharge[k]==1 && muonpt[k]>20 && muonpt[k]<=25){
			subposmuon +=1;
			if (subposmuon ==1){
			spm_pt = muonpt[k];
			spm_eta = muoneta[k];
			spm_phi =muoneta[k];
			}
			if(subposmuon ==2){
                        spm_pt1 = muonpt[k];
                        spm_eta1 = muoneta[k];
                        spm_phi1 =muoneta[k];
                        }}
			if(muoncharge[k]==-1 && muonpt[k]>25){
			leadnegmuon +=1;
			if(leadnegmuon ==1){
			lnm_pt = muonpt[k];
			lnm_eta = muoneta[k];
			lnm_phi = muonphi[k];
			}
			if(leadnegmuon ==2){
                        lnm_pt1 = muonpt[k];
                        lnm_eta1 = muoneta[k];
                        lnm_phi1 = muonphi[k];
                        }}
			if(muoncharge[k]==-1 && muonpt[k]>20 && muonpt[k]<=25){
			subnegmuon +=1;
			if(subnegmuon ==1){
			snm_pt = muonpt[k];
			snm_eta= muoneta[k];
			snm_phi= muonphi[k];
			}
                        if(subnegmuon ==2){
                        snm_pt1 = muonpt[k];
                        snm_eta1= muoneta[k];
                        snm_phi1= muonphi[k];
                        }}}
		for(Long64_t k=0;k<8;++k){
			if(eleccharge[k]==1 && elecpt[k]>25){
			leadposelec +=1;
			if(leadposelec ==1){
			lpe_pt = elecpt[k];
			lpe_eta = eleceta[k];
			lpe_phi =elecphi[k];
			}
                        if(leadposelec ==2){
                        lpe_pt1 = elecpt[k];
                        lpe_eta1 = eleceta[k];
                        lpe_phi1 =elecphi[k];
                        }}
			if(eleccharge[k]==1 && elecpt[k]>20 && elecpt[k]<=25){
			subposelec +=1;
			if( subposelec ==1){
			spe_pt = elecpt[k];
			spe_eta = eleceta[k];
			spe_phi = elecphi[k];
			}
			if(subposelec ==2){
			spe_pt1 =elecpt[k];
			spe_eta1=eleceta[k];
			spe_phi1=elecphi[k];
			}}
			if(eleccharge[k]==-1 && elecpt[k]>=25){
			leadnegelec +=1;
			if(leadnegelec ==1){
			lne_pt=elecpt[k];
			lne_eta = eleceta[k];
			lne_phi = elecphi[k];
			}
                        if(leadnegelec ==2){
                        lne_pt1 = elecpt[k];
                        lne_eta1 = eleceta[k];
                        lne_phi1 = elecphi[k];
                        }}
                        if(eleccharge[k]==-1 && elecpt[k]>20 && elecpt[k]<=25){
			subnegelec +=1;
			if(subnegelec ==1){
                        sne_pt=elecpt[k];
                        sne_eta = eleceta[k];
                        sne_phi = elecphi[k];
                        }
                        if(subnegelec ==2){
                        sne_pt1 = elecpt[k];
                        sne_eta1 = eleceta[k];
                        sne_phi1 = elecphi[k];
                        }}
		}
/// this section calculates the invariant mass of different potential pairs of leptons 

		float M_lpm_lnm = sqrt(2*lpm_pt*lnm_pt*(cosh(lpm_eta-lnm_eta)-cos(lpm_phi-lnm_phi)));
		float M_lpm_snm = sqrt(2*lpm_pt*snm_pt*(cosh(lpm_eta-snm_eta)-cos(lpm_phi-snm_phi)));
		float M_spm_lnm = sqrt(2*spm_pt*lnm_pt*(cosh(spm_eta-lnm_eta)-cos(spm_phi-lnm_phi)));


///
		float M_lpm_sne = sqrt( 2*lpm_pt*sne_pt*(cosh(lpm_eta-sne_eta)-cos(lpm_phi-sne_phi)));
		float M_lnm_spe = sqrt(2*lnm_pt*spe_pt*(cosh(lnm_eta-spe_eta)-cos(lnm_phi-sne_phi)));
		float M_lpm_lne = sqrt(2*lpm_pt*lne_pt*(cosh(lpm_eta-lne_eta)-cos(lpm_phi-lne_phi)));
		float M_lnm_lpe = sqrt(2*lnm_pt*lpe_pt*(cosh(lnm_eta-lpe_eta)-cos(lnm_phi-lpe_phi)));
		float M_lne_spm = sqrt(2*lne_pt*spm_pt*(cosh(lne_eta-spm_eta)-cos(lne_phi-spm_phi)));
		float M_lpe_snm = sqrt(2*lpe_pt*snm_pt*(cosh(lpe_eta-snm_eta)-cos(lpe_phi-snm_phi)));
///
		float M_lpe_lne = sqrt(2*lpe_pt*lne_pt*(cosh(lpe_eta-lne_eta)-cos(lpe_phi-lne_phi)));
		float M_lpe_sne = sqrt(2*lpe_pt*sne_pt*(cosh(lpe_eta-sne_eta)-cos(lpe_phi-sne_phi)));
		float M_spe_lne = sqrt(2*spe_pt*lne_pt*(cosh(spe_eta-lne_eta)-cos(spe_phi-lne_phi)));

// this is an accept/reject condition bases on the different pairs one wants and the range in which their invariant
// mass falls
		if ((leadnegelec==1 && leadposelec==1 && M_lpe_lne >20 && (M_lpe_lne >105 or M_lpe_lne < 75) && Q==2) or
		   (leadposelec==1 && subnegelec==1 && M_lpe_sne>20 && (M_lpe_sne > 105 or M_lpe_sne <75) && Q==2) or
		   (leadposmuon==1 && leadnegmuon==1 && M_lpm_lnm>20 && (M_lpm_lnm > 105 or M_lpm_lnm <75) && Q==2) or
		   (leadposmuon==1 && subnegmuon==1 && M_lpm_snm>20 && (M_lpm_snm > 105 or M_lpm_snm <75) && Q==2) or
		   (subposmuon==1 && leadnegmuon==1 && M_spm_lnm>20 && (M_spm_lnm > 105 or M_spm_lnm <75) && Q==2) or
		   (leadposmuon==1 && subnegelec==1 && Q==2) or
		   (leadnegmuon==1 && subposelec==1 && Q==2) or
		   (leadposmuon==1 && leadnegelec==1 && Q==2) or
		   (leadnegmuon==1 && leadposelec==1 && Q==2) or
		   (leadnegelec==1 && subposmuon==1 && Q==2) or
		   (leadposelec==1 && subnegmuon==1 && Q==2) or
		   (Q==2 && subposelec ==1 && leadnegelec==1 && M_spe_lne >20 && ( M_spe_lne<75 or M_spe_lne>105)) ) {
		Long64_t jetnumber =4;
		for (Long64_t l=0;l<nJet;++l){
		pt_misstotx += Jet_pt[l]*cos(Jet_phi[l]); 
		pt_misstoty += Jet_pt[l]*sin(Jet_phi[l]);
		}
		pt_missx = pow(pt_misstotx +pt_misselecx+pt_missmuonx,2);
		pt_missy = pow(pt_misstoty +pt_misselecy+pt_missmuony,2);      
		pt_miss = sqrt(pt_missx+pt_missy);
		if (nJet>= jetnumber){
		for (Long64_t NumJet= 0;NumJet<nJet;++NumJet){
			counter = counter + 1;
			CSVV2[counter] = Jet_btagCSVV2[NumJet];
			if (abs(Jet_eta[NumJet])<2.4 && Jet_pt[NumJet] > 30) p+=1;
			if (abs(Jet_eta[NumJet])<2.4 && Jet_pt[NumJet] > 25){
			Discrim[NumJet] = CSVV2[counter];
			}
		}
		for(Long64_t i = 0;i < 40; ++i){
			if (Discrim[i] > x_l) L+=1;
			if (Discrim[i] > x_m) M+=1;
			if (Discrim[i] > x_t) T+=1;
			if (Discrim[i] > x_i) mine+=1;
		}
		for(Long64_t i = 0;i < 40; ++i){
			if(Discrim[39] < Discrim[i]) Discrim[39] = Discrim[i];
		}
		repeat =0;
		for(Long64_t i = 0;i < 40; ++i){
			if(Discrim[i]== Discrim[39]) repeat += 1;
			if(repeat==3) Discrim[38]=Discrim[39];
			if(Discrim[38] < Discrim[i] && Discrim[i]<Discrim[39]) Discrim[38] = Discrim[i];
		}
		repeat =0;
		if (Discrim[39]==Discrim[38]) repeat =-2;
		for(Long64_t i = 0;i < 40; ++i){
			if( Discrim[i]==Discrim[38]) repeat +=1;
			if( repeat ==3) Discrim[37]=Discrim[38];
			if(Discrim[37] < Discrim[i] && Discrim[i]<Discrim[38]) Discrim[37] = Discrim[i];
		}
		First_Discrim = Discrim[39];
		Second_Discrim = Discrim[38];
		Third_Discrim = Discrim[37];
		}

else {
for (Long64_t NumJet= 0;NumJet<nJet;++NumJet) {
counter=counter+1;
CSVV2[counter]=-10;
CMVA[counter]=-10;
}}
if ( T>=0 && M>=2  && L>=2 && SumTrans>500 && p>=4 && pt_miss>50){
 ne=ne+1;
 MyHist->Fill(M);}
}
if (ientry < 0) break;
nb=fChain->GetEntry(jentry);nbytes+=nb;
   
}
TCanvas *c1 = new TCanvas("c1");
c1->SetLogy();
MyHist->Draw();

TFile myGraphFile("yourhist.root","RECREATE");
myGraphFile.cd();
MyHist->Write();
myGraphFile.Close();

}

Code For Combining Two Histograms

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

float error(float Evtot,float EvTag){   \\ This is a montecarlo simulation of errors on each data point and is used later
Int_t newb=0;
float x=1;
float q = 6;
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
        for (Long64_t inc =1;inc<(1000);++inc){
                if ((EvTag-inc)<=0) break;
                float D = 0;
                float Prob = (EvTag-inc)/Evtot;
                float NoB[1000]={};
                std::binomial_distribution<int> distribution (Evtot,Prob);
//              cout<< distribution(generator)<< endl;
                for(Long64_t i=0; i<1000;++i){
                        NoB[i] = distribution(generator);
//                      cout<< "NoB="<< NoB[i]<<"    B'="<<B+inc<<"             D="<< D<<endl;
                        if (NoB[i] > EvTag)
                        D=D+1;
                }

                if (D/1000<=(x/q)){
                newb =( EvTag -inc); \\ gives value of errors
                break;
                }



        }
return (abs(newb-EvTag));

}

TH1F *graph_41t;
TH1F *graph_21t;
TH1F *MyHist4= new TH1F("MyHist4","4Top",3,2,5); // Define your histograms and their size
TH1F *MyHist2= new TH1F("MyHist2","2Top",3,2,5);
void HistGrapher()
{
TFile file_21t("LAST2.root");   
TFile file_41t("LAST4.root");
graph_21t = (TH1F*) file_21t.Get("MyHist");  // import the graphs you want to combine
graph_41t = (TH1F*) file_41t.Get("MyHist");
if (graph_21t == NULL)
{
std::cout << "2tGraph not found. Exiting!" << std::endl;
return;
}
if (graph_41t == NULL)
{
std::cout << "4tGraph not found. Exiting!" << std::endl;
return;
}

Long64_t NumBins4 = 1+graph_41t->GetNbinsX(); \\ the +1 makes sure every bin is counted
Long64_t NumBins2 = 1+graph_21t->GetNbinsX();

for(Long64_t i=0;i<NumBins4;++i){
float Content=graph_41t->GetBinContent(i); 
float centre =graph_41t->GetBinCenter(i);
for (Long64_t j=0;j<Content;++j){ 
MyHist4->Fill(centre,0.7232);
}
Long64_t Total = MyHist4->GetEntries();
Long64_t Bin_Content = MyHist4->GetBinContent(i);
float Bin_Error = error(Total,Content);
MyHist4->SetBinError(i,Bin_Error);
}

for(Long64_t i=0;i<NumBins2;++i){
Long64_t Content=graph_21t->GetBinContent(i);
float centre =graph_21t->GetBinCenter(i);
for (Long64_t j=0;j<Content;++j){
MyHist2->Fill(centre,469.3);
}
Long64_t Total = MyHist2->GetEntries();
Long64_t Bin_Content = MyHist2->GetBinContent(i);
float Bin_Error = error(Total,Content);
MyHist2->SetBinError(i,Bin_Error);
}


MyHist4->SetFillColor(kOrange);
MyHist2->SetFillColor(kAzure);
MyHist3->SetFillColor(kBlack);
TCanvas *c1 = new TCanvas("c1");
//c1->SetLogy();
//c1->DrawFrame(4,0,8,200000);
THStack *hs = new THStack("hs","CSL");
hs->Add(MyHist4);
hs->Add(MyHist2);
//hs->Add(MyHist3);
hs->SetMinimum(0.4);
gPad->SetLogy();
hs->Draw(" HIST E2");
hs->GetXaxis()->SetTitle("N^{m}_{b}");
hs->GetYaxis()->SetTitle("Events");
gPad->BuildLegend(0.75,0.75,0.9,0.9,"");

//TFile myGraphFile("CMS2019HistoPaperCopyCSVV2STACKEDdisc.root","RECREATE");
//myGraphFile.cd();
//hs->Write();
//myGraphFile.Close();
}

Code For Combining Two Histograms into a Signal vs sqrt(Background) graph

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

float error(float Evtot,float EvTag){
Int_t newb=0;
float x=1;
float q = 6;
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
        for (Long64_t inc =1;inc<(1000);++inc){
                if ((EvTag-inc)<=0) break;
                float D = 0;
                float Prob = (EvTag-inc)/Evtot;
                float NoB[1000]={};
                std::binomial_distribution<int> distribution (Evtot,Prob);
//              cout<< distribution(generator)<< endl;
                for(Long64_t i=0; i<1000;++i){
                        NoB[i] = distribution(generator);
//                      cout<< "NoB="<< NoB[i]<<"    B'="<<B+inc<<"             D="<< D<<endl;
                        if (NoB[i] > EvTag)
                        D=D+1;
                }

                if (D/1000<=(x/q)){
                newb =( EvTag -inc);
                break;
                }



        }
return (abs(newb-EvTag));

}

TH1F *graph_41t;
TH1F *graph_21t;
TH1F *MyHist4= new TH1F("MyHist4","4Top",24,40,280);
TH1F *MyHist2= new TH1F("MyHist2","2Top",24,40,280);
void SignalBackgroundGrapher()
{
TFile file_21t("TMM2Toplep_ptOS.root");
TFile file_41t("TMM4Toplep_ptOS.root");
graph_21t = (TH1F*) file_21t.Get("MyHist");
graph_41t = (TH1F*) file_41t.Get("MyHist");
if (graph_21t == NULL)
{
std::cout << "2tGraph not found. Exiting!" << std::endl;
return;
}
if (graph_41t == NULL)
{
std::cout << "4tGraph not found. Exiting!" << std::endl;
return;
}

float SB[300]={}; float SB60[300]={}; float SB100[300]={};
float Con2[300]={}; float Con260[300]={};float Con2100[300]={};
float Con4[300]={}; float Con460[300]={};float Con4100[300]={};
Long64_t i=300;
float x[300]={};
float Fill4[300]={}; float Fill460[300]={}; float Fill4100[300]={};
float Fill2[300]={}; float Fill260[300]={}; float Fill2100[300]={};
float w4 = 0.7232; float w2= 436.9;
float w460 = 1.21; float w260 =732.2;
float w4100 = 2.01; float w2100=1220.3;
float Fill2error[300]={}; float Fill260error[300]={}; float Fill2100error[300]={};
float Fill4error[300]={}; float Fill460error[300]={}; float Fill4100error[300]={};
float SBerror[300]={}; float SBerror60[300]={}; float SBerror100[300]={};
float tot2 =65555;
float tot4= 595;
for(Long64_t j=0;j<i;++j){
float Content4 = graph_41t->GetBinContent(j);
float Content2 = graph_21t->GetBinContent(j);
Con4[j] =w4*Content4;
Con460[j] = w460*Content4;
Con4100[j] = w4100*Content4;
Con2[j] =w2 * Content2;
Con260[j] = w260*Content2;
Con2100[j] = w2100*Content2;
x[j]= 40+ j*10; 
}
for(Long64_t k=0;k<i;++k){
for(Long64_t l=k;l<i;++l){
Fill4[k] +=Con4[l]; Fill460[k] += Con460[l]; Fill4100[k] += Con4100[l];
Fill2[k] +=Con2[l]; Fill260[k] += Con260[l]; Fill2100[k] += Con2100[l];
}
Fill2error[k]= error(Fill2[0],Fill2[k]); Fill260error[k]= error(Fill260[0],Fill260[k]); Fill2100error[k]= error(Fill2100[0],Fill2100[k]);
Fill4error[k]= error(Fill4[0],Fill4[k]); Fill460error[k]= error(Fill460[0],Fill460[k]); Fill4100error[k]= error(Fill4100[0],Fill4100[k]);
if(Fill2[k] != 0){
SB[k]=Fill4[k]*pow(Fill2[k],-0.5);
SBerror[k]=SB[k]*sqrt(pow((Fill4error[k])/Fill4[k],2)+0.5*pow(Fill2error[k]/Fill2[k],2));
}
else SB[k]=0;
if(Fill260[k] != 0){
SB60[k]=Fill460[k]*pow(Fill260[k],-0.5);
SBerror60[k]=SB60[k]*sqrt(pow((Fill460error[k])/Fill460[k],2)+0.5*pow(Fill260error[k]/Fill260[k],2));
//cout<<SBerror[k]<<endl;
}
else SB60[k]=0;
if(Fill2100[k] != 0){
SB100[k]=Fill4100[k]*pow(Fill2100[k],-0.5);
SBerror100[k]=SB100[k]*sqrt(pow((Fill4100error[k])/Fill4100[k],2)+0.5*pow(Fill2100error[k]/Fill2100[k],2));
//cout<<SBerror[k]<<endl;
}
else SB100[k]=0;
cout<< "35 = " << SB[k] << "	60 = " << SB60[k] << " 	SB100 = " << SB100[k] << endl;
}

TCanvas *c1 = new TCanvas("c1","S/sqrt(B) vs no. extra B-jets at loose constraint",200,10,700,500);
c1->DrawFrame(40,0,160,2,"OS; Lepton p_{T};S/#sqrt{B}");
TGraph *gr1 = new TGraphErrors(i,x,SB,0,SBerror);
TGraph *gr2 = new TGraphErrors(i,x,SB60,0,SBerror60);
TGraph *gr3 = new TGraphErrors(i,x,SB100,0,SBerror100);
gr1 -> SetLineColorAlpha(kBlue,0.6);
gr2 -> SetLineColorAlpha(kGreen,0.6);
gr3 -> SetLineColorAlpha(kBlack,0.6);
 TLegend *legend = new TLegend(0.75,0.75,0.9,0.9); 				// option "C" allows to center the header
   legend->AddEntry(gr1,"35.9fb^{-1}","l");
   legend->AddEntry(gr2,"60fb^{-1}","l");
   legend->AddEntry(gr3,"100fb^{-1}","l");
gr1->Draw();
gr2->Draw();
gr3->Draw();

legend->Draw();
//TFile myGraphFile("CMS2019HistoPaperCopyCSVV2STACKEDdisc.root","RECREATE");
//myGraphFile.cd();
//hs->Write();
//myGraphFile.Close();
}


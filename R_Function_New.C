#include "TH1.h"
#include "TFile.h"
#include "TCut.h"
#include "TCutG.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TROOT.h"
#include "math.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TLine.h"
#include <iostream>
Bool_t verbose=kFALSE;
#define nplanes 6
using namespace std;

class R_function {
private:
  TString plane_name[nplanes];
  Double_t r_plane[nplanes];

public:
  R_function();
  R_function(TString input_file); 
  //Double_t ph[100],dp[100],th[100],y[100];
  TCutG *mycut[nplanes];
  enum planes_enum{ph_dp,ph_th,ph_y,th_dp,th_y,y_dp};
  Double_t getR2D(int select);
  Double_t R_function_2D(Double_t x_ev, Double_t y_ev, int select);
  Double_t Global_R_function(Double_t ph, Double_t dp, Double_t th, Double_t y);  
  Double_t minimumdistance(Double_t x_ev, Double_t y_ev, Double_t xi, Double_t yi, Double_t xj, Double_t yj);
  Double_t return_min_r(Double_t previous_R, Double_t thi_R);
  TString GetPlaneName(Int_t j);
};

Double_t R_function::getR2D(int select)
{
  return r_plane[select];
}

TString R_function::GetPlaneName(Int_t j)
{
  TString out;
  if(j<nplanes)
    out=plane_name[j];
  else
    out=Form("This plane number %d doesn't exist",j);
  
  return out;
}

Double_t R_function::return_min_r(Double_t previous_R, Double_t this_R)
{
  Double_t selectR=-99999999;
  if(previous_R*this_R>=0) // they are both of the same sign
    if(fabs(previous_R)<fabs(this_R)) // choose the one with the smallest absolute value
      selectR=previous_R;
    else
      selectR=this_R;
  else // both R have different signs
    if(previous_R<this_R) // choose the smallest one (with sign) 
      selectR=previous_R;
    else
      selectR=this_R;
 
  // if(verbose)cout<<" chosen R is ="<<selectR<<"\n =======\n";
  
  return selectR;  
}
  
R_function()	
{
  plane_name[0]= "ph_dp";
  plane_name[1]=  "ph_th";
  plane_name[2]=  "ph_y";
  plane_name[3]=  "th_dp";
  plane_name[4]=  "th_y";
  plane_name[5]=  "y_dp";
  TString inputfile= "run12518cuts.root";
  TFile *myfile = TFile::Open(inputfile);
  for(int i=0;i<nplanes;i++)
    {
      mycut[i]= (TCutG*) myfile->Get("mycut_"+plane_name[i]);
    }
  myfile->Close();
  delete myfile;
}
  
 
R_function::R_function(TString inputfile)
{
  plane_name[0]= "ph_dp";
  plane_name[1]=  "ph_th";
  plane_name[2]=  "ph_y";
  plane_name[3]=  "th_dp";
  plane_name[4]=  "th_y";
  plane_name[5]=  "y_dp";
  TFile *myfile = TFile::Open(inputfile);
  for(int i=0;i<nplanes;i++)
    {
      mycut[i]= (TCutG*) myfile->Get("mycut_"+plane_name[i]);
    }
  myfile->Close();
  delete myfile;
}

Double_t R_function::R_function_2D(Double_t x_ev, Double_t y_ev, int select)
{
  
  Double_t nseg=mycut[select]->GetN()-1;
  Double_t *x_cut=mycut[select]->GetX();
  Double_t *y_cut=mycut[select]->GetY();
  

  /*
  Double_t xmax=-999999.;
  Double_t xmin=999999.;
  Double_t ymax=-999999.;
  Double_t ymin=999999.;

  for(int i=0;i<nseg;i++)
    {
    // Double_t x_cut=mycut[select]->GetX();
    // Double_t y_cut=mycut[select]->GetY();
 
    
      if (x_cut>xmax) xmax=x_cut;
      if (x_cut>xmin) xmin=x_cut;
      if (y_cut>ymax) ymax=y_cut;
      if (y_cut>ymin) ymin=y_cut;
    }
    cout<< xmax<<endl;
    cout<< ymax<<endl;
    
  for(int i=0;i<nseg;i++)
    {
      x_cut[i]=x_cut[i]/(xmax-xmin);
      y_cut[i]=y_cut[i]/(ymax-ymin);

      x_ev=x_ev/(xmax-xmin);
      y_ev=y_ev/(ymax-ymin);
    }
    */
    
 ///////////////////////////
   
  
  
  Double_t min2D=9999999;
  Double_t min;
  verbose=kTRUE;
  if(verbose)
    {
     // cout<<"select="<<select<<endl;r_
     // cout<<"x,y="<<x_ev<<"  "<<y_ev<<endl;
    }
  for(int i=0;i<nseg;i++)
    {
      min = minimumdistance(x_ev,y_ev,x_cut[i],y_cut[i],x_cut[i+1],y_cut[i+1]); 
      if(min<min2D) min2D = min;
      if(verbose)
	{
  // cout<<" i="<<i<<" min distance="<<min<<"  x_cut[i],y_cut[i],x_cut[i+1],y_cut[i+1]"<<x_cut[i]<<"  "<<y_cut[i]<<"  "<<x_cut[i+1]<<"  "<<y_cut[i+1]<<endl;
	}
    }
   // min2D=min2D/sqrt(mycut[select]->Area());//((xmax-xmin)*(ymax-ymin)));
  
  Int_t inout=mycut[select]->IsInside(x_ev,y_ev); 
  if(inout==0)inout=-1;
  return inout*min2D;
}


 
Double_t R_function::minimumdistance(Double_t x_ev, Double_t y_ev,
				     Double_t xi, Double_t yi,
				     Double_t xj, Double_t yj)
{
  Double_t mindistance=9999999;
  Double_t ABAE=(xj-xi)*(x_ev-xi)+(yj-yi)*(y_ev-yi); //cosine of the  angle between AB (the segment) and AE
  Double_t BABE=(xi-xj)*(x_ev-xj)+(yi-yj)*(y_ev-yj); //cosine of the  angle between BA (the segment) and BE
  Double_t AB=sqrt(pow(xj-xi,2)+pow(yj-yi,2));
  Double_t AE=sqrt(pow(x_ev-xi,2)+pow(y_ev-yi,2));
  Double_t BE=sqrt(pow(x_ev-xj,2)+pow(y_ev-yj,2));
  Double_t dist;
  
  if(ABAE<0 && BABE<0)
    {
      cout<<"this is trouble\n";
    }
  else if(ABAE<=0)// if the cosine is negative the angle is larger than 90deg
    mindistance=AE;       
  else if(BABE<=0)
    mindistance=BE;
  else
    {
      Double_t s=(AB+AE+BE)*0.5;
      Double_t area=sqrt(s*(s-AE)*(s-BE)*(s-AB));
      mindistance=2*area/AB;	  
    }
 
  return mindistance; 
}									      										      
Double_t  R_function::Global_R_function(Double_t ph, Double_t dp, Double_t th, Double_t y)
{
  Double_t global_R=9999999;
  Double_t thisR=R_function_2D(dp,ph,ph_dp);
  global_R=return_min_r(thisR,global_R);
  r_plane[ph_dp]=thisR;
  thisR=R_function_2D(th,ph,ph_th);
  global_R=return_min_r(thisR,global_R);
  r_plane[ph_th]=thisR;
  thisR=R_function_2D(y,ph,ph_y);
  global_R=return_min_r(thisR,global_R);
  r_plane[ph_y]=thisR;
  thisR=R_function_2D(dp,th,th_dp);
  global_R=return_min_r(thisR,global_R);
  r_plane[th_dp]=thisR;
  thisR=R_function_2D(y,th,th_y);
  global_R=return_min_r(thisR,global_R);
  r_plane[th_y]=thisR;
  thisR=R_function_2D(dp,y,y_dp);
  global_R=return_min_r(thisR,global_R);
  r_plane[y_dp]=thisR;
  
  return global_R;
}

/////////////////////////////////////////////////////////////////////////
int r_function1(Bool_t bool_w=kTRUE) {
 R_function *rfunc= new R_function("run12518cuts.root");
  Double_t ph,dp,th,y;
  Double_t ph1,dp1,th1,y1;
  Double_t R=0;
  Double_t weight;
  Int_t events=100;
  gRandom->SetSeed();
  Double_t r[nplanes];
  TFile *f = new TFile("test_new_code.root","RECREATE");
  TH2F *phvdp = new TH2F("phvdp","phvdp",1000,-0.07,0.07,1000,-0.07,0.07);
  TH2F *phvth = new TH2F("phvth","phvth",1000,-0.07,0.07,1000,-0.07,0.07);
  TH2F *phvy = new TH2F("phvy","phvy",1000,-0.07,0.07,1000,-0.07,0.07);
  TH2F *thvdp = new TH2F("thvdp","thvdp",1000,-0.07,0.07,1000,-0.07,0.07);
  TH2F *thvy = new TH2F("thvdy","thvdy",1000,-0.07,0.07,1000,-0.07,0.07);
  TH2F *yvdp = new TH2F("yvdp","yvdp",1000,-0.07,0.07,1000,-0.07,0.07);
  TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","dp:th:ph:y:R:rphdp:rphth:rphy:rthdp:rthy:rydp:weight");
  TH1F *ph33 = new TH1F("ph33", "ph", 100, -0.1, 0.1);
  TH1F *dp33 = new TH1F("dp33", "dp", 100, -0.1, 0.1);
  TH1F *th33= new TH1F("th33", "th", 100, -0.1, 0.1);
  TH1F *y33 = new TH1F("y33", "y", 100, -0.1, 0.1);

    //events=1;
  for (int i=0;i<events;i++)
    {
     // if(i%1000)cout<<"so far "<<i<<" out of "<<events<<endl; 
     dp= gRandom->Uniform(-0.05,0.07);
     th= gRandom->Uniform(-0.07,0.07);
     ph= gRandom->Uniform(-0.05,0.05);
     y= gRandom->Uniform(-0.05,0.05);
     if(bool_w)
     weight=1+((ph+0.1)/0.1*2.);
      //weight=1+(sin(ph/0.035*3.14));
      //weight=1-(sin(ph/0.1*3.14));
       
     else
       weight=1.;
      
      
      R = rfunc->Global_R_function( ph,  dp,  th,  y);      
      
      for(int j=0;j<nplanes;j++)
	{
	  r[j]=rfunc->getR2D(j);
	}
      //verbose=kFALSE;
     //if(verbose) {
	// cout<<"rphdp=r[0]="<<r[0]<<endl;
	// cout<<"rphth=r[1]="<<r[1]<<endl;
	// cout<<"rphy=r[2]="<<r[2]<<endl;
	// cout<<"rthdp=r[3]="<<r[3]<<endl;
	// cout<<"rthy=r[4]="<<r[4]<<endl;
	// cout<<"rydp=r[5]="<<r[5]<<endl;
	 cout<<"R="<<R<<endl;
      //}
      if(R > 0.0)
	{
	  phvdp->Fill(dp,ph);
	  phvth->Fill(th,ph);
	  phvy->Fill(y,ph);
	  thvdp->Fill(dp,th);
	  thvy->Fill(y,th);
	  yvdp->Fill(dp,y);
	}			   
     float ntu[]=
	    {
	   dp,th,ph,y,R,r[0],r[1],r[2],r[3],r[4],r[5],weight
	    };
	    
	    ntuple->Fill(ntu);
	    
      
 
  ph33->Fill(ph,weight);
  dp33->Fill(dp,weight);
  th33->Fill(th,weight);
  y33->Fill(y,weight);
 	  
      
   }
  f->Write();

  return 0;
}

  
///////////////////////////////////R_function from Simulation///////////////////////////////////////////////////////
  
int r_function_sim(TString arg_inputfile="arg_siminputfile481.txt")//inputfile is different for each kinematic
{

 
  TString cut_inputfile;
  TString simulation_inputfile;
  Double_t th_HRS;
  Double_t p_HRS;
  TString outputfile;
  Int_t kin;
  Double_t xBn;
  Double_t q2n;
  Double_t En;
  Double_t E_scattern;
  Double_t R_test;
  Double_t n_gen;
  Double_t z_mzx;
  Double_t z_min;  
  //// read the input file
  // TString inputfilename=Form("%s.txt",arg_inputfile.Data());
  FILE *finput= fopen(arg_inputfile.Data(),"r");
  TString tline;
  Int_t counter=0;
  
  while(tline.Gets(finput)  )
    {
      TString localtstring=(tline.Strip(TString::kBoth,' '));
      localtstring.Remove(TString::kBoth,'\n');
      localtstring.Remove(TString::kBoth,'\r');
      if(localtstring(0,2)=="//")
	continue;
      // at this point only reading lines that do not have two leading
      //, that is not comment-lines
      if(counter==0)
	simulation_inputfile=Form("%s",localtstring.Data());
      else if(counter==1)
	{
	  th_HRS=atof(localtstring.Data());
	  th_HRS = th_HRS*TMath::DegToRad();
	}
      else if(counter==2)
	p_HRS=atof(localtstring.Data());
      else if(counter==3)
	cut_inputfile=Form("%s",localtstring.Data());
      else if(counter==4)
	outputfile=Form("%s",localtstring.Data());
      else if(counter==5)
	kin=atoi(localtstring.Data());
	 else if(counter==6)
	xBn=atof(localtstring.Data());
	 else if(counter==7)
	q2n=atof(localtstring.Data());
	 else if(counter==8)
	En=atof(localtstring.Data());
	 else if(counter==9)
	E_scattern=atof(localtstring.Data());
	 else if(counter==10)
	R_test=atof(localtstring.Data());
	 else if(counter==11)
	n_gen=atof(localtstring.Data());
	else if(counter==12)
	z_max=atof(localtstring.Data());
	else if(counter==13)
	z_min=atof(localtstring.Data());
      counter+=1;
    }
    cout<<"sim = "<<simulation_inputfile<<endl;
    cout<<"theta = "<<th_HRS<<endl;
    cout<<"momentum = "<<p_HRS<<endl;
    cout<<"output = "<<outputfile<<endl;
    cout<<"kinematic = "<<kin<<endl;
    cout<<"cut input file ="<<cut_inputfile<<endl;
    cout<<"xBn ="<<xBn<<endl;
    cout<<"q2n ="<<q2n<<endl;
    cout<<"En ="<<En<<endl;
    cout<<"E_sscattern ="<<E_scattern<<endl;
    cout<<"R_test ="<<R_test<<endl;
    cout<<"n_gen ="<<n_gen<<endl;
    cout<<"z_max ="<<z_max<<endl;
    cout<<"z_min ="<<z_min<<endl;
    

  R_function *rfunc= new R_function(cut_inputfile);
  Double_t ph2[100],dp2[100],th2[100],y2[100];
  Double_t R=0;
  Double_t x_v = 0;
  ///////////////////////
  
  Double_t kx;
  Double_t ky;
  Double_t kz;
  Double_t theta;
  Double_t psf_dis;
  Double_t q2;
  Double_t theta;
  
    /////////////////////////////
    
  Double_t qx_v;
  Double_t qy_v;
  Double_t qz_v;
  Double_t qx_tg;
  Double_t qy_tg;
  Double_t qz_tg;
  Double_t z_react;
  Double_t y_tg;
  Double_t x_tg;
  Double_t phi_tg;
  Double_t th_tg;
  Double_t delta;
  Double_t p; 
  Double_t xB;
  Double_t r[nplanes];
  Double_t alpha=1./137.;
  Double_t sum=0;
  Double_t sum2=0;
  Int_t count=0;
  Int_t npassrcut=0;
  //////////////////////////variables forcalculating DIS_cross_section//////////////////////////////////  
  // Double_t deadtime_corr= 0.985;// different foe each kinematics 
  // Double_t n_s2m= ;// different for each kinematics
   //Double_t n_pr= ;// different for each kinematics
  ///////////////////////////////////////////////////////////////////////////////// 
  TFile *f = new TFile(Form("%s",outputfile.Data(),kin),"RECREATE");
  TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","dp:th:ph:y:R:rphdp:rphth:rphy:rthdp:rthy:rydp:z_react:DIS:psf_dis:weight");
  TFile*MCData = new TFile(Form("/adaqfs/home/a-onl/dvcs12/onlana/gula/%s_%d.root",simulation_inputfile.Data(),kin));
  TTree *T = (TTree*)MCData->Get("tr1");
  tr1->SetBranchAddress("psf",&psf_dis); 
  tr1->SetBranchAddress("kpx_v",&qx_v);
  tr1->SetBranchAddress("kpy_v",&qy_v);
  tr1->SetBranchAddress("kpz_v",&qz_v);
  tr1->SetBranchAddress("kx_v",&kx);
  tr1->SetBranchAddress("ky_v",&ky);
  tr1->SetBranchAddress("kz_v",&ky);
  tr1->SetBranchAddress("vz_v",&z_react);
  tr1->SetBranchAddress("q2_v",&q2);
  tr1->SetBranchAddress("xB_v",&xB);
   
  //------------- a list of slopes and intercepts for xb = 0.35,0.45,0.55,0.65
  //       xb           slope           intercept
  //
  //
  //
  //     0.08          0.03619           0.290
  //     0.125        -0.00285           0.35681
  //     0.175         0.01286           0.32135
  //     0.25          0.00339           0.3194
  //     0.35         -0.00437           0.28132
  //     0.45         -0.00631           0.22217
  //     0.55         -0.00771           0.16741
  //     0.65         -0.00424           0.09744
  //     0.75         -0.00152           0.04254
  //     0.85         -0.00028           0.01039
  //-------------------------------------------------------------------
  //
  //  want F2(x,Q2) = m(x)*Q2 + b(x) when x is not one of the above, but close
  //
  //*************************************************************
  //   m(x) = "ms" * x + "mi"
  //   find "ms" and "mi" for the different intervals in listed xb:
  //
  //-----------------------------------------------------------------
  //
  //     interval             "ms"              "mi"
  //       
  //    0.35->0.45          -0.0194           0.00242 
  //    0.45->0.55          -0.014           -0.00001
  //    0.55->0.65           0.0347          -0.026795
  //    0.65->0.75           0.0272          -0.02192
  //    0.75->0.85           0.0124          -0.01082
  //
  //*************************************************************
  //   b(x) = "bs" * x + "i"
  //   find "bs" and "bi" for the different intervals in listed xb:
  //
  //-----------------------------------------------------------------
  //
  //     interval             "bs"              "bi"
  //       
  //    0.35->0.45          -0.5915           0.4888 
  //    0.45->0.55          -0.5476           0.4686
  //    0.55->0.65          -0.6997           0.5522
  //    0.65->0.75          -0.549            0.45429
  //    0.75->0.85          -0.3215           0.283665
  // 
  //************************************************************

    Double_t x_test[10] = {0.08,0.125,0.175,0.25,0.35,0.45,0.55,0.65,0.75,0.85};
  Double_t ms_array[10] = {-0.86755,0.3142,-0.12626,-0.0776,-0.0194,-0.014,0.0347,0.0272,0.0124,0.0124};
  Double_t mi_array[10] = {0.1055,-0.042125,0.03495,0.02279,0.00242,-0.00001,-0.026795,-0.02192,-0.01082,-0.01082};
  Double_t bs_array[10] = {1.4847,-0.7092,-0.026,-0.3808,-0.5915,-0.5476,-0.6997,-0.549,-0.3215,-0.3215};
  Double_t bi_array[10] = {0.17122,0.44546,0.3259,0.4146,0.4888,0.4686,0.5522,0.45429,0.283665,0.283665};
  Double_t ms;
  Double_t mi;
  Double_t bs;
  Double_t bi; 
 

  
  ////////////////////////////////////

 Int_t nentries = tr1->GetEntries();
 // Int_t nentries = 2;
  cout<<"nentries="<<nentries<<endl;
 
  for(Int_t n=0;n<nentries;n++)
    {
      tr1->GetEntry(n);
      qx_tg = -qy_v;
      qy_tg = TMath::Cos(th_HRS)*qx_v - TMath::Sin(th_HRS)*qz_v;
      qz_tg = TMath::Cos(th_HRS)*qz_v + TMath::Sin(th_HRS)*qx_v;
      y_tg = (TMath::Cos(th_HRS)*x_v - TMath::Sin(th_HRS)*z_react)*0.01;
      th_tg = qx_tg/qz_tg;
      phi_tg = qy_tg/qz_tg;
      p = TMath::Sqrt((TMath::Power(qx_v,2)+TMath::Power(qy_v,2)+TMath::Power(qz_v,2)));
      delta = (p-p_HRS)/p_HRS;
      
      //////weight//////
      Double_t E_scatter=sqrt(pow(qx_v,2) + pow(qy_v,2) + pow(qz_v,2));
      Double_t E=sqrt(pow(kx,2)+pow(ky,2)+pow(kz,2));
      Double_t theta2 = atan(sqrt(pow(qx_v,2) + pow(qy_v,2))/qz_v);
      /////////////////////////////different for each kinematic//////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////// 
        Double_t nu = q2/(2*xB);
        Double_t  y_dis = nu/E;
        Double_t F2;
   /////////////// nominal kinematic_different for each kinematic///////////////////////  
   /////////////////////////////////////////////////////////// 
   
      Double_t nu_n = q2n/(2*xBn);
      Double_t  y_disn = nu_n/En;
      Double_t F2n;


      for(Int_t k=0; k<10; k++){
	      if(k<9){
		if(xBn<=x_test[0])
		  F2n = 0.036*xBn+0.29;
		if(xBn <= x_test[k+1] && xBn >= x_test[k]){
		  //cout<<"here     "<<j<<endl;
		  ms = ms_array[k];
		  mi = mi_array[k];
		  bs = bs_array[k];
		  bi = bi_array[k];
		  F2n = (ms*xBn + mi)*q2n + (bs*xBn + bi);
		}
	      }
	      if(k==9){
		if(xBn >= x_test[k]){
		  //cout<<"here     "<<j<<endl;
		  F2n= -0.00028*q2n + 0.01039;
		}	    
	      }
	    }
	    
//////////////////////////////////////////////////////////////	    
////////////////////////////////////////////////////////////
      for(Int_t j=0; j<10; j++){
	      if(j<9){
		if(xB<=x_test[0])
		  F2 = 0.036*xB+0.29;
		if(xB <= x_test[j+1] && xB >= x_test[j]){
		  //cout<<"here     "<<j<<endl;
		  ms = ms_array[j];
		  mi = mi_array[j];
		  bs = bs_array[j];
		  bi = bi_array[j];
		  F2 = (ms*xB + mi)*q2 + (bs*xB + bi);
		}
	      }
	      if(j==9){
		if(xB >= x_test[j]){
		  //cout<<"here     "<<j<<endl;
		  F2 = -0.00028*q2 + 0.01039;
		}	    
	      }
	    }
	   

	    Double_t DIS = ((2*TMath::Pi()*alpha*alpha)/(q2*q2*xB))*F2*(2-(2*y_dis) + (y_dis*y_dis));  //DIS = (2*alpha*alpha/Q2xb)/(Q2*Q2)*F2*(2-2*y_dis + y_dis*y_dis); //arxiv nuber 9703012
	    Double_t weight = DIS*10E3;
	               // cout<<weight<<endl;	
     
                 
       
	
      y2[0]=y_tg;
      ph2[0]=phi_tg;
      th2[0]=th_tg;
      dp2[0]=delta;

      verbose=kFALSE;
      if ( z_react>z_min && z_react<z_max )
      {
      R = rfunc->Global_R_function( ph2[0],  dp2[0],  th2[0], y2[0]); 
      for(int j=0;j<nplanes;j++)
	{
	  r[j]=rfunc->getR2D( j );
	} 
	//cout <<"R="<<R<<endl;
	
	
	
	   /////////////////////////////nominal kinematic/////////////////////////// 
	   /////////////////////////////////////////////////////////////////////////////
	   
	   if(R>R_test)// different for each kinematics
	   { 
         sum+=psf_dis;
        sum2+=DIS*psf_dis;
         npassrcut+=1; 
       }
      ntuple->Fill(dp2[0],th2[0],ph2[0],y2[0],R,r[0],r[1],r[2],r[3],r[4],r[5],z_react,DIS,psf_dis,weight);
      }
     
	  }
	  /////////////////calculating the DIS_cross_section//////////////////////////////////
	  
   Double_t DIS_psf=(sum/n_gen);
   cout <<" DIS_psf = "<< DIS_psf<<endl; 
   cout <<" DIS_MC=  "<<sum<<endl;
   Double_t DIS_nominal = ((2.*TMath::Pi()*alpha*alpha)/(q2n*q2n*xBn))*F2*(2.-(2.*y_disn) + (y_disn*y_disn));
   cout <<" DIS_nominal= "<< DIS_nominal<<endl;
   Double_t alpha_coeff= sum2/(sum* DIS_nominal);
   cout << "alpha_coeff = " << alpha_coeff << endl; 
   cout<<" saved events= "<< nentries<<"       "<<"npassrcut= "<< npassrcut<<endl; 
   Double_t me=0.5110E-03;//GeV
   Double_t logterm = TMath::Log(q2n/(me*me));
   Double_t pie=TMath::Pi();
   Double_t delta_0_R=(alpha/(pie))*(-0.95 -pie*pie/3.+(logterm*logterm)/2.);// +(alpha_fs/(pie))*TMath::Log((deltaE*deltaE)/(k*kprime))*(logterm-1.);
   cout <<"delta_0_R= "<< delta_0_R<<endl;                                  
   Double_t delta_ver=(alpha/(pie))*((3./2.)*logterm-2.+(pie*pie/6.)-(logterm*logterm)/2.);
   cout <<"delta_vertex= "<< delta_ver<<endl;
   Double_t delta_vac= (1.*alpha/(3.*pie))*(logterm-(5./3.));
   cout << "delta_vacccum= "<<delta_vac<<endl;
   Double_t n_exp=1.;        //= d_deadtime*n_cer*n_s0*n_s2m*n_pr; 
   Double_t n_virt = exp(delta_0_R + delta_ver)/(1.- delta_vac)**2.;
   cout <<"n_virt= "<< n_virt<<endl;

	 ////////////////////////////////////////////////////////////////////////////////////
  
  f->Write();
  f->Close();
  MCData->Close();
  return 0;
}
  

/////////////////////////////////R_function from data////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int r_function_from_data(TString arg_datainputfile="arg_data481.txt")
{
  TString cut_inputfile;
  TString basic;
  Int_t runi;
  Int_t kin;     
  //// read the input file
  FILE *finput= fopen(arg_datainputfile.Data(),"r");
  TString tline;
  Int_t counter=0;
  while(tline.Gets(finput)  )
    {
      TString localtstring=(tline.Strip(TString::kBoth,' '));
      localtstring.Remove(TString::kBoth,'\n');
      localtstring.Remove(TString::kBoth,'\r');
      if(localtstring(0,2)=="//")
	continue;
      // at this point only reading lines that do not have two leading //, that is not comment-lines
      if(counter==0)
      {
      cut_inputfile=Form("%s",localtstring.Data());
      cout <<"COUNTER = : "<< counter<<" cut_inputfile= "<< cut_inputfile<< endl;
      }
      else if(counter==1) 
      {
      basic=Form("%s",localtstring.Data());
       cout <<"COUNTER = : "<< counter<<" basic= "<< basic<< endl;
       }
      else if(counter==2)
      {
     runi=atoi(localtstring.Data());
      cout <<"COUNTER = : "<< counter<<" run_numbe = "<< runi<< endl;
      }
      else if(counter==3)
      {
       kin=atoi(localtstring.Data());
        cout<<"COUNTER = : "<< counter<<" kin= "<< kin<< endl;
        }
        
        counter+=1;
        }	
 R_function *rfunc= new R_function(Form("%s.root",cut_inputfile.Data()));
 Double_t ph2[100],dp2[100],th2[100],y2[100];
 Double_t R=0;
 Double_t r[nplanes];
 TFile *f = new TFile(Form("%s_%d.root",basic.Data(),kin),"RECREATE");
 TH2F *phvdp = new TH2F("phvdp","phvdp",1000,-0.07,0.07,1000,-0.07,0.07);
 TH2F *phvth = new TH2F("phvth","phvth",1000,-0.07,0.07,1000,-0.07,0.07);
 TH2F *phvy = new TH2F("phvy","phvy",1000,-0.07,0.07,1000,-0.07,0.07);
 TH2F *thvdp = new TH2F("thvdp","thvdp",1000,-0.07,0.07,1000,-0.07,0.07);
 TH2F *thvy = new TH2F("thvy","thvy",1000,-0.07,0.07,1000,-0.07,0.07);
 TH2F *yvdp = new TH2F("yvdp","yvdp",1000,-0.07,0.07,1000,-0.07,0.07);
 ///////////////////////////////////////////////////////////  
 
 

   
  
 TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","dp:th:ph:y:R:rphdp:rphth:rphy:rthdp:rthy:rydp:ver:cher:track:t3coin:prsum:pr1val:prl1:prl2:u1clst:u2clst:v1clst:v2clst:dis:triggerPatternWord:missing:single:edtm",300);
 TChain *tree=new TChain("T");
 tree->Add(Form("/adaqfs/home/a-onl/dvcs12/onlana/Rootfiles/left_dvcs_%d*.root",runi));
 
 Double_t cher,track;
 //TCut c="L.cer.asum_c>500 &&L.tr.n==1";
 Double_t nnevent= tree->GetEntries();
 Double_t cher,track,vertex,prl1,prl2,u1clst,u2clst,v1clst,v2clst;
 UInt_t triggerPatternWord, edtm;
 Double_t w1;
 Double_t w2;
 Double_t t3coin[50];
 Double_t wa1[9]= {2.47,2.68,3.33,1.33,3.41,2.51,2.90,2.97,2.64}; // 36_1, 36_2, 36_3, 48_1, 48_2, 48_3, 48_4,60_1,60_3
 Double_t wa2[9]={3.35,3.72,4.51,2.13,4.61,3.48,4.01,4.08,3.63} ;
 Int_t count=0;
 Int_t count1=0;
 Double_t prsum;
 Double_t pr1val;
 
if(runi>=10553 && runi<=10760)
    {
      w1= wa1[0];
      w2=wa2[0];
    }
  else if (runi>=14150 && runi<=14260)
    {
      w1= wa1[1];
      w2=wa2[1];
    }
  else if (runi>=14476 && runi<=14525)
    {
      w1= wa1[2];
      w2=wa2[2];
    }
  else if (runi>=12508 && runi<=12647)
    {
      w1= wa1[3];
      w2=wa2[3];
    }
  else if (runi>=13000 && runi<=13015 ||runi>=13183 && runi<=13236 )
    {
      w1= wa1[4];
      w2=wa2[4];
    }
  else if (runi>=12838 && runi<=12991)
    {
      w1= wa1[5];
      w2=wa2[5];
    }
  else if (runi>=13100 && runi<=13162 ||runi>=13279 && runi<=13418 )
    {
      w1= wa1[6];
      w2=wa2[6];
    }
    else if (runi==15017)
    {
      w1= wa1[7];
      w2=wa2[7];
    }
    else if (runi==14628 )
    {
      w1= wa1[8];
      w2=wa2[8];
    }
  else
    cout<<"Your run may not be the good one please double check it"<<endl;
  cout<<"w1 = "<<w1<<"   w2 = "<<w2<<endl;
   
  

  cout << "Entries in inputfile " << nnevent << endl;
 tree->SetBranchAddress("L.tr.tg_th",&th2);
 tree->SetBranchAddress("L.tr.tg_ph",&ph2);
 tree->SetBranchAddress("L.tr.tg_dp",&dp2);
 tree->SetBranchAddress("L.tr.tg_y",&y2);
 tree->SetBranchAddress("L.cer.asum_c",&cher);
 tree->SetBranchAddress("L.tr.n",&track);
 tree->SetBranchAddress("rpl.z",&vertex);
 tree->SetBranchAddress("L.prl1.asum_c",&prl1);
 tree->SetBranchAddress("L.prl2.asum_c",&prl2);
 tree->SetBranchAddress("L.vdc.u1.nclust",&u1clst);
 tree->SetBranchAddress("L.vdc.u2.nclust",&u2clst);
 tree->SetBranchAddress("L.vdc.v1.nclust",&v1clst);
 tree->SetBranchAddress("L.vdc.v2.nclust",&v2clst);
 tree->SetBranchAddress("dvcs_scaler_EDTM",&edtm);
 tree->SetBranchAddress("triggerPatternWord",&triggerPatternWord);
 tree->SetBranchAddress("DL.t3",&t3coin);
 
 tree->SetBranchStatus("*",0);
 tree->SetBranchStatus("L.tr.tg_th",1);
 tree->SetBranchStatus("L.tr.tg_ph",1);
 tree->SetBranchStatus("L.tr.tg_dp",1);
 tree->SetBranchStatus("L.tr.tg_y",1);
 tree->SetBranchStatus("L.cer.asum_c",1);
 tree->SetBranchStatus("L.tr.n",1);
 tree->SetBranchStatus("rpl.z",1);
 tree->SetBranchStatus("L.prl1.asum_c",1);
 tree->SetBranchStatus("L.prl2.asum_c",1);
 tree->SetBranchStatus("L.vdc.u1.nclust",1);
 tree->SetBranchStatus("L.vdc.u2.nclust",1);
 tree->SetBranchStatus("L.vdc.v1.nclust",1);
 tree->SetBranchStatus("L.vdc.v2.nclust",1);
 tree->SetBranchStatus("dvcs_scaler_EDTM",1);
 tree->SetBranchStatus("triggerPatternWord",1);
  tree->SetBranchStatus("DL.t3",1);
     // for(Int_t i=0;i<1000;i++) 
        for(Int_t i=0;i<nnevent;i++) 
      
	{
	  tree->GetEntry(i);
	  UInt_t dis = (triggerPatternWord &0x00080)==128;
	  // UInt_t dis = (triggerPatternWord &0x00100)==256; is only for kin60
      UInt_t missing  =(triggerPatternWord & 0xffffffc0)==2048;
      UInt_t single = (triggerPatternWord & 0x3f);
      prsum= (prl1/w1 + prl2/w2); // for pion rejector sum 
      pr1val=prl1/w1; // for pion rejector layer first
	  if (cher>150 && vertex>-0.06 && vertex<0.07 && prsum>600 && pr1val>200 && edtm==0 && track ==1 && u1clst==1 && dis==1 && fabs(dp2[0])<0.1 && fabs(ph2[0])<0.1 && fabs(th2[0])<0.15 && fabs(y2[0])<0.15 ) 
	
	  
	  {           
	  R = rfunc->Global_R_function( ph2[0],  dp2[0],  th2[0], y2[0]);
	  for(int j=0;j<nplanes;j++)
	    {
	     r[j]=rfunc->getR2D( j );
		}
	
	   for(int j=0;j<nplanes;j++) 
	   if(R > 0.004)
	    {
	    phvdp->Fill(dp2[0],ph2[0]);
	    phvth->Fill(th2[0],ph2[0]);
	    phvy->Fill(y2[0],ph2[0]);
	    thvdp->Fill(dp2[0],th2[0]);
	    thvy->Fill(y2[0],th2[0]);
	    yvdp->Fill(dp2[0],y2[0]);
	     }
	     		  
	    verbose=kFALSE;
	    float dis1[]=
	    {
	    dp2[0],th2[0],ph2[0],y2[0],R,r[0],r[1],r[2],r[3],r[4],r[5],vertex,cher,track, t3coin[0], prsum, pr1val,prl1, prl2,u1clst,u2clst,v1clst,v2clst,dis,triggerPatternWord,missing,single,edtm
	    };
	    
	    ntuple->Fill(dis1); 
     Bool_t dis_events_selection = kFALSE;
      // This is the condition for the DIS electron selection
    if(cher>150 && track==1 && vertex>-0.06 && vertex<0.07 && prsum>600 && pr1val>200 && edtm==0 && track ==1 && u1clst==1 && u2clst==1 && v1clst==1 && v2clst==1 && R> 0.004 ) dis_events_selection= kTRUE;
     //  if(R>0.0) dis_events_selection= kTRUE;
    if (dis_events_selection== kTRUE && dis ==1)// after with prescale and dead time  correction add with missing, if any
 	 {
 	count++;
     }  
  if (dis_events_selection== kTRUE && missing ==1 && single > 50 &&t3coin[0]>1300 && t3coin[0]<1600 )//selecting missing dis events
     {
    count1++; // no any prescale correction is needed 
     }	
 	 }
     } 
    cout<<"total number of events= "<< nnevent<<"       "<<"DIS events that pass R-cut= "<< count<<endl; 
    cout <<"missed = "<<count1<<endl; 	
    f->Write();
    return 0;
}
 


 
 

  
 

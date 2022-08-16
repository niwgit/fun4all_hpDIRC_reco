#include "/work/eic3/users/nwickjlb/hpDIRC_reco_test/fun4all_hpDIRC_reco/geometric_reco/geo_reco_f4a.h"
#include <iostream>
#include "/work/eic3/users/nwickjlb/hpDIRC_reco_test/fun4all_hpDIRC_reco/prttools/prttools.cpp"

using namespace std;

TGraph gg_gr;

void geo_reco_f4a(TString infile, TString lutfile, TString filedir, int verbose)
{
  fVerbose = verbose;

  //gROOT->ProcessLine(".L PrtLutNode.cxx+");
  gROOT->ProcessLine(".L /work/eic3/users/nwickjlb/hpDIRC_reco_test/fun4all_hpDIRC_reco/geometric_reco/PrtLutNode.cxx+");

  TH1F*  fHist1 = new TH1F("Time1","1", 1000,0,20);
  TH1F*  fHist2 = new TH1F("Time2","2", 1000,-10,10);
  TH1F*  fHistDiff[3];
  TH2F*  fDiff = new TH2F("diff",";measured time [ns];t_{measured}-t_{calc} [ns]", 500,0,100,150,-5,5);
  TH2F*  fHist4 = new TH2F("Time4","4", 200,-1,1, 200,-1,1);
  TH2F*  fHist5 = new TH2F("Time5","5", 200,-1,1, 200,-1,1);
  TH1F*  fFindTime = new TH1F("ft",";t_{measured}-t_{calculated} [ns];entries [#]",2000,-100,100);
  TH1F*  fFindTimeA[20];
  TH1F*  fFindTimeRes = new TH1F("ftr","ftr",100,-2,2);
  TH2F*  fdtt = new TH2F("dtt",";t_{measured}-t_{calculated} [ns];#theta_{l} [deg]", 1000,-2,2, 1000,0,90);
  TH2F*  fdtl = new TH2F("dtl",";t_{measured}-t_{calculated} [ns];path length [m]", 1000,-2,2, 1000,0,15);
  TH2F*  fdtp = new TH2F("dtp",";#theta_{l} [deg];path length [m]", 1000,0,90, 1000,0,15);
  TH2F*  fhChrom = new TH2F("chrom",";t_{measured}-t_{calculated} [ns];#Delta#theta_{C} [mrad]", 100,-2,2, 100,-30,30);
  TH2F*  fhChromL = new TH2F("chroml",";(t_{measured}-t_{calculated})/L_{path};#Delta#theta_{C} [mrad]", 100,-0.0002,0.0002, 100,-30,30);
  TH1F*  fHistMcp[28];
  double fCorr[28];

  TH1F* Hist_lut_theta = new TH1F("lut_theta","lut theta (rad)", 100, -TMath::Pi(), TMath::Pi());
  //TF1* chrom_func = new TF1("chrom_func","0.02619 + (-0.000592973)*x + (4.05258e-06)*x*x + (1.15741e-09)*x*x*x",20,150);

  int gg_i(0);
  
  fChain = new TChain("mG4EvtTree");
  fChain->Add(infile);

  fFile = new TFile(lutfile);
  fTree=(TTree *) fFile->Get("prtlut") ;

  for(Int_t b=0; b < 10; b++){
    fLut[b] = new TClonesArray("PrtLutNode");
    fTree->SetBranchAddress(Form("LUT_%d",b),&fLut[b]); 
  }
  
  fTree->GetEntry(0);

  fFit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +[3]",0.35,0.9);
  fSpect = new TSpectrum(10);
  fMethod = 2;
  fCriticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm) 

  int col[]={kRed+1,kBlue+1,kBlack};
  for(int i=0; i<3; i++){
    fHistDiff[i] = new TH1F(Form("TimeDiff_%d",i),";t_{measured}-t_{calculated} [ns];entries [#]", 500,-10,10);
    fHistDiff[i]->SetLineColor(col[i]);
  }  
  for(int h=0; h<5; h++){
    hthetac[h] = new TH1F(Form("thetac_%d",h),";#theta_{C} [rad];entries [#]", 200,0.75,0.9);
    hthetac[h]->SetLineColor(prt_color[h]);
    hthetacd[h] = new TH1F(Form("thetacd_%d",h),";#Delta#theta_{C} [mrad];entries [#]", 200,-60,60);
    hthetacd[h]->SetLineColor(prt_color[h]);
    hnph[h] = new TH1F(Form("nph_%d",h),";detected photons [#];entries [#]", 200,0,200);
    hnph[h]->SetLineColor(prt_color[h]);
    fFunc[h] = new TF1(Form("gaus_%d",h),"[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);    
    fFunc[h]->SetLineColor(prt_color[h]);
    }

    for(int i=0; i<20; i++){
    fFindTimeA[i] = new TH1F(Form("fta_%d",i),";t_{measured}-t_{calculated} [ns];entries [#]",1000,-10,10);
  }

  for(int i=0; i<28; i++){
    fHistMcp[i] = new TH1F(Form("fHistMcp_%d",i),Form("fHistMcp_%d;#theta_{C} [rad];entries [#]",i), 50,-0.05,0.05);
    }

  // read corrections  
  
  fCorrFile = infile+"_corr.root";
  for(int i=0; i<prt_nmcp; i++) fCorr[i]=0;
  if(!gSystem->AccessPathName(fCorrFile)){
    std::cout<<"------- reading  "<<fCorrFile <<std::endl;
    int pmt;
    double corr,cspr[5];
    TChain ch; ch.SetName("corr"); ch.Add(fCorrFile);
    ch.SetBranchAddress("pmt",&pmt);    
    ch.SetBranchAddress("corr",&corr);
    ch.SetBranchAddress("cspr",&cspr);
    for(int i=0; i<ch.GetEntries(); i++){
      ch.GetEvent(i);
      fCorr[pmt] = (fabs(corr)<0.007)? corr: 0.00001;
      fSigma = 0.001*cspr[2]*0.9;
      std::cout<<"pmt "<<pmt<<"  "<<fCorr[pmt]<< " spr = "<<fSigma<<std::endl;    
    }
  }else{
    std::cout<<"------- corr file not found  "<<fCorrFile <<std::endl;
  }
  

  TVector3 dird, dir, momInBar(0,0,1), posInBar;
  double mom, tangle, tdiff, evtime, bartime, lenz, dirz, luttheta, hitTime;
  int tofPid(0);
  bool reflected = kFALSE;
  gStyle->SetOptFit(111);
  
  TVector3 fnX1 = TVector3(1,0,0);   
  TVector3 fnY1 = TVector3(0,1,0);
  int nsHits(0),nsEvents(0);

  double theta(0),phi(0),cangle[5] = {0}, spr[5] = {0}, nph[5] = {0}, sigma_nph[5] = {0}, timeRes(0), test1(0), sep(0), sep_err(0);

  TGaxis::SetMaxDigits(3);
  prt_setRootPalette(1);
  prt_initDigi(2);
  
  TFile fRootFile("reco.root","recreate");

  TTree tree("reco","reco");
  tree.Branch("mom", &mom,"mom/D");
  tree.Branch("tofPid", &tofPid,"tofPid/I");
  tree.Branch("spr", &spr,"spr[5]/D");
  tree.Branch("nph",&nph,"nph[5]/D");
  tree.Branch("sigma_nph",&sigma_nph,"sigma_nph[5]/D");
  tree.Branch("cangle",&cangle,"cangle[5]/D");
  tree.Branch("sep",&sep,"sep/D");
  tree.Branch("sep_err",&sep_err,"sep_err/D");
  tree.Branch("timeres",&timeRes,"timeRes/D");
  tree.Branch("test1",&test1,"test1/D");
  tree.Branch("theta",&theta,"theta/D");
  tree.Branch("phi",&phi,"phi/D");
  
  timeRes = 0.2;
  test1 = 0.0005;
  //test1 = 0.;

  int nEvents = fChain->GetEntries();
  int end = nEvents;
  //std::cout << "No of events = " << nEvents << std::endl;
  //if(end==0) end = nEvents;

  const int arr_size = 500;

  int hit_size = 0;
  int Particle_id = 0;

  fChain->SetBranchAddress("nhits", &hit_size);

  int mcp_num[arr_size], pixel_id[arr_size], bar_id[arr_size];
  Double_t px, py, pz, theta_ang;
  double hit_pos[arr_size][3];
  double hit_mom[arr_size][3];
  double lead_time[arr_size];
  Long_t hit_pathId[arr_size];
  int nrefl[arr_size];

  fChain->SetBranchAddress("theta", &theta_ang);
  fChain->SetBranchAddress("px", &px);
  fChain->SetBranchAddress("py", &py);
  fChain->SetBranchAddress("pz", &pz);
  fChain->SetBranchAddress("bar_id", &bar_id);
  fChain->SetBranchAddress("mcp_id", &mcp_num);
  fChain->SetBranchAddress("pixel_id", &pixel_id);
  fChain->SetBranchAddress("pid", &Particle_id);
  fChain->SetBranchAddress("hit_pos", &hit_pos);
  fChain->SetBranchAddress("hit_mom", &hit_mom);
  fChain->SetBranchAddress("lead_time", &lead_time);
  fChain->SetBranchAddress("hit_pathId", &hit_pathId);
  fChain->SetBranchAddress("nrefl", &nrefl);

  double timeCut = 0.5;
  
  //std::cout<<"Run started for ["<<start<<","<<end <<"]"<<std::endl;

  int fEvId = 2030; 

  { // identify particles involved 
    fp1=0;
    fp2=0;
    for(int i=0; i<nEvents; i++){
      fChain->GetEntry(i);
      int pid = prt_get_pid(Particle_id);
      if(fp1==0) fp1=pid;
      if(fp2==0 && fp1!=pid) {
	fp2=pid;
	break;
      }
    }
    if(fabs(fp2)<fabs(fp1)){
      int t = fp1;
      fp1=fp2;
      fp2=t;      
      }
  }  

  int ntotal=0;  


  for (int ievent=0; ievent < nEvents; ievent++){
    fChain->GetEntry(ievent);
    int nHits = hit_size;  
    if(ievent%1000==0) std::cout<<"Event # "<< ievent << " has "<< nHits <<" hits"<<std::endl;
    
    theta = theta_ang*TMath::RadToDeg();
    ntotal+=nHits;
    prt_theta = theta;
    mom = TMath::Sqrt(px*px + py*py + pz*pz);
    //double corr_coeff = chrom_func->Eval(theta);
    
    tofPid = Particle_id;
    int pid = prt_get_pid(tofPid);
    int tnph[5]={0};

    double minChangle = 0.35;
    double maxChangle = 0.9;
    TVector3 mom_vec = TVector3(px,py,pz);
    //int barbox_number = (int) (mom_vec.phi()*TMath::RadToDeg()+15.0)/30.0;
    //mom_vec.RotateZ(-barbox_number*30*TMath::DegToRad());
    TVector3 rotatedmom = mom_vec.Unit();
    double sum1(0),sum2(0),noise(0.2);
    
    // track smearing
    TVector3 init = rotatedmom;
    rotatedmom.RotateY(prt_rand.Gaus(0,test1));
    rotatedmom.Rotate(TMath::Pi(),init);
 
    if(fSigma<0.003) fSigma=0.007;  

    if(ievent==1){
      double range = 160;
      if(mom>1.5) range = 100;
      if(mom>2) range = 60;
      if(mom<1) range = 300;
      if(mom<0.7) range = 400;
      if(mom<0.6) range = 500;

      if(fp2==3 && theta<120) range = 80;
      if(fp2==3 && mom<4) range = 500;
      
      for(int h=0; h<5; h++){
	fLnDiff[h] = new TH1F(Form("LnDiff_%d",h),  ";ln L(K) - ln L(#pi);entries [#]",100,-range,range);
	//fLnDiff[h] = new TH1F(Form("LnDiff_%d",h),  ";ln L(K) - ln L(#pi);entries [#]",100,-60,60);
	fLnDiff[h]->SetLineColor(prt_color[h]);
	}
    }
        
    for(int i=0; i<5; i++){
      fAngle[i] = acos(sqrt(mom*mom + prt_mass[i]*prt_mass[i])/mom/1.4738); //1.4738 = 370 = 3.35
      fFunc[i]->SetParameter(0,1);
      fFunc[i]->SetParameter(1,fAngle[i]);
      fFunc[i]->SetParameter(2,fSigma);
    }

    // hits loop
    for(int h=0; h < nHits; h++)
      {      
	hitTime = lead_time[h] + prt_rand.Gaus(0,0.1);
	lenz = 2555 + hit_pos[h][2]; // ECCE z-shift
	dirz = hit_mom[h][2]; 

	int barId = bar_id[h];
	int mcp = mcp_num[h];
	int pix = pixel_id[h];	
	int ch =  256*mcp + pix;
      
	TVector3 hit_mom_vec = TVector3(hit_mom[h][0], hit_mom[h][1], hit_mom[h][2]);
	TVector3 dir0 = hit_mom_vec.Unit();      
	TVector3 cz = TVector3(-rotatedmom.X(),rotatedmom.Y(),rotatedmom.Z());
	TVector3 cd = TVector3(-dir0.X(),dir0.Y(),dir0.Z());    
	TVector3 unitdir1 = rotatedmom.Unit();
	TVector3 unitdir2 = rotatedmom.Unit();
	cz.RotateUz(unitdir1);
	cd.RotateUz(unitdir2);
      
	double phi0 =  cd.Phi();
	if(dirz>0){
	  reflected = true;
	  lenz = 2*4235 - lenz; // ECCE
	}
	else{
	  reflected = false;
	  // continue;
	}
	
	double theta0 = rotatedmom.Angle(dir0);
	
	fHist5->Fill(theta0*TMath::Sin(phi0),theta0*TMath::Cos(phi0));

	PrtLutNode *node = (PrtLutNode*) fLut[barId]->At(ch);
	int size = node->Entries();
	bool isGoodHit(false);
      	
	int nreflections = nrefl[h];

	Long_t hpath = hit_pathId[h]; 
	//std::cout << "path id = " << hpath << std::endl;
	TString spath = Form("%ld",hpath);	
	
	for(Int_t i=0; i<size; i++){
	  dird = node->GetEntry(i);
	  Long_t lpath = node->GetPathId(i);
	  TString slpath = Form("%ld",lpath);
	  bool ipath=0;
	  if(hpath==lpath) ipath=1;
	  
	  evtime = node->GetTime(i);
	  for(int u=0; u<4; u++){
	    if(u == 0) dir = dird;
	    if(u == 1) dir.SetXYZ( dird.X(),-dird.Y(), dird.Z());
	    if(u == 2) dir.SetXYZ(-dird.X(), dird.Y(), dird.Z());
	    if(u == 3) dir.SetXYZ(-dird.X(),-dird.Y(), dird.Z());
	    if(reflected) dir.SetXYZ( dir.X(), dir.Y(),-dir.Z());  
	    	    
	    if(dir.Angle(fnX1) < fCriticalAngle || dir.Angle(fnY1) < fCriticalAngle) continue;

	    luttheta = dir.Theta();

	    if(luttheta > TMath::PiOver2()) luttheta = TMath::Pi()-luttheta;
	    Hist_lut_theta->Fill(luttheta);
      
	    bartime = (lenz/cos(luttheta)/198.5); //198

	    fHist1->Fill(hitTime);
	    double luttime = bartime+evtime;
	    tdiff = hitTime-luttime;
	    fHistDiff[reflected]->Fill(tdiff);
	    if(ipath) fHistDiff[2]->Fill(tdiff);

	    tangle = rotatedmom.Angle(dir)+fCorr[mcp];//45;
	    //tangle = rotatedmom.Angle(dir);
	    //if(tangle>TMath::PiOver2()) tangle = TMath::Pi()-tangle;

	    if(fabs(tdiff)<2) tangle -= 0.008*tdiff; // chromatic correction
	    //if(fabs(tdiff)<2) tangle -= 0.005*tdiff;
	    //if(fabs(tdiff)<2) tangle -= corr_coeff*tdiff;

	    if(fabs(tdiff)>timeCut+luttime*0.035) continue;
	    fDiff->Fill(hitTime,tdiff);
	      	      
	    hthetac[pid]->Fill(tangle);
	    hthetacd[pid]->Fill((tangle-fAngle[pid])*1000);    
	    fHistMcp[mcp]->Fill(tangle-fAngle[pid]);
	    fhChrom->Fill(tdiff,(tangle-fAngle[pid])*1000);
	    fhChromL->Fill(tdiff/(lenz/cos(luttheta)),(tangle-fAngle[pid])*1000);
	      
	    if(fabs(tangle-fAngle[fp2]) > 0.05 && fabs(tangle-fAngle[fp1]) > 0.05) continue;

	    if(tangle > minChangle && tangle < maxChangle){
	      TVector3 rdir = TVector3(-dir.X(),dir.Y(),dir.Z());
	      TVector3 unitdir3 = rotatedmom.Unit();
	      rdir.RotateUz(unitdir3);
	      double cphi =  rdir.Phi();	      
	      fHist4->Fill(tangle*TMath::Sin(cphi),tangle*TMath::Cos(cphi));
	      gg_gr.SetPoint(gg_i,tangle*TMath::Sin(cphi),tangle*TMath::Cos(cphi));
	      gg_i++;
	    }

	    isGoodHit=true;
	      
	    sum1 += -TMath::Log(fFunc[fp1]->Eval(tangle)+noise);
	    sum2 += -TMath::Log(fFunc[fp2]->Eval(tangle)+noise);
	  }
	}
      
	if(isGoodHit){
	  nsHits++;
	  tnph[pid]++;
	  if(pid == 2) prt_hdigi[mcp]->Fill(pix%16, pix/16);
	}
      }

    double sum = sum1-sum2;

    if(sum!=0) fLnDiff[pid]->Fill(sum);
    if(tnph[pid]>1) hnph[pid]->Fill(tnph[pid]);    
  
  if(fVerbose==1){
      prt_canvasAdd("ff",800,400);
      if(hthetac[fp1]->GetMaximum()>0) hthetac[fp1]->Scale(1/hthetac[fp1]->GetMaximum());
      hthetac[fp1]->Draw("hist");
      fFunc[fp1]->Draw("same");
      fFunc[fp2]->Draw("same");      
      prt_waitPrimitive("ff");
      prt_canvasDel("ff");

      FindPeak(cangle,spr);
      theta = theta_ang*TMath::RadToDeg();
      hthetac[fp1]->Reset();
    }
    
    if(++nsEvents>=end) break;
    
  }

  
    FindPeak(cangle,spr);

    for(int h=0; h<5; h++){
      if(hnph[h]->GetEntries()<20) continue;
      hnph[h]->Fit("gaus","","S",5,180);
      auto f = hnph[h]->GetFunction("gaus");
      if(f){
	nph[h] = f->GetParameter(1);
	sigma_nph[h] = f->GetParameter(2);
      }

    
      theta = theta_ang*TMath::RadToDeg();

    TF1 *ff;
    double m1=0,m2=0,s1=0,s2=0,dm1=0,dm2=0,ds1=0,ds2=0; 
    if(fLnDiff[fp2]->GetEntries()>10){
      fLnDiff[fp2]->Fit("gaus","S");
      ff = fLnDiff[fp2]->GetFunction("gaus");
      m1=ff->GetParameter(1);
      s1=ff->GetParameter(2);
      dm1=ff->GetParError(1);
      ds1=ff->GetParError(2);

      if(fp1==0 && mom <1.5){ //handle tails      
	fLnDiff[fp2]->Fit("gaus","S","",m1-1.8*s1,500);
	ff = fLnDiff[fp2]->GetFunction("gaus");
	m1=ff->GetParameter(1);
	s1=ff->GetParameter(2);
	dm1=ff->GetParError(1);
	ds1=ff->GetParError(2);
      }
      
    }
    if(fLnDiff[fp1]->GetEntries()>10){
      fLnDiff[fp1]->Fit("gaus","S");
      ff = fLnDiff[fp1]->GetFunction("gaus"); 
      m2=ff->GetParameter(1);
      s2=ff->GetParameter(2);
      dm2=ff->GetParError(1);
      ds2=ff->GetParError(2);

      if(fp1==0 && mom <1.5){ ///handle tails
	fLnDiff[fp1]->Fit("gaus","S","",-500,m2+1.8*s2);
	ff = fLnDiff[fp1]->GetFunction("gaus");      
	m2=ff->GetParameter(1);
	s2=ff->GetParameter(2);
	dm2=ff->GetParError(1);
	ds2=ff->GetParError(2);
      }
    }
    sep = (fabs(m1-m2))/(0.5*(s1+s2));
     
    double e1,e2,e3,e4;
    e1 =  2/(s1 + s2)*dm1;
    e2 = -2/(s1 + s2)*dm2;
    e3 = -2*(m1 - m2)/((s1 + s2)*(s1 + s2))*ds1;
    e4 = -2*(m1 - m2)/((s1 + s2)*(s1 + s2))*ds2;
    sep_err = sqrt(e1*e1+e2*e2+e3*e3+e4*e4);    
    
    std::cout<<Form("%3d : SPR = %2.2f N = %2.2f +/- %2.2f",prt_pdg[fp1],spr[fp1],nph[fp1],sigma_nph[fp1])<<std::endl;
    std::cout<<Form("%3d : SPR = %2.2f N = %2.2f +/- %2.2f",prt_pdg[fp2],spr[fp2],nph[fp2],sigma_nph[fp2])<<std::endl;
    std::cout<<Form("SEP = %2.2f +/- %2.2f ",sep,sep_err)<<std::endl;    
    }

  if(!fVerbose) gROOT->SetBatch(1);

  
  tree.Fill();
  tree.Write();  
  Hist_lut_theta->Write();
  
  if(fVerbose>1){
    TString nid = Form("_%1.2f_%1.4f_%1.2f",prt_theta,test1,mom);

    { // cherenkov angle
      prt_canvasAdd("tangle"+nid,800,400);
      prt_normalize(hthetac, 5);
          
      hthetac[fp1]->SetTitle(Form("theta %1.2f", prt_theta));
      hthetac[fp1]->Draw("");
      hthetac[fp2]->Draw("same");
      drawTheoryLines(6);

      prt_canvasAdd("tangled"+nid,800,400);
      prt_normalize(hthetacd, 5);    
      hthetacd[fp2]->SetTitle(Form("theta %1.2f", prt_theta));
      hthetacd[fp2]->Draw("");
      hthetacd[fp1]->Draw("same");
    }      

    { // nph
      prt_canvasAdd("nph"+nid,800,400);      
      prt_normalize(hnph, 5);    
      hnph[fp1]->SetStats(0);
      hnph[fp1]->Draw();
      hnph[fp2]->Draw("same");
    }
    
    { // sep
      fLnDiff[fp1]->SetStats(0);
      fLnDiff[fp2]->SetStats(0);
          
      prt_canvasAdd("lh"+nid,800,400);
      prt_normalize(fLnDiff, 5);
      fLnDiff[fp2]->SetName(Form("s_%2.2f",sep));
      fLnDiff[fp2]->SetTitle(Form("separation = %2.2f #pm %2.2f s.d.",sep, sep_err));
      fLnDiff[fp2]->GetXaxis()->SetTitle("ln L("+prt_lname[fp2]+") - ln L("+prt_lname[fp1]+")");      
      fLnDiff[fp2]->Draw();
      fLnDiff[fp1]->Draw("same");
    }
    
    { // chromatic corrections
      prt_canvasAdd("chrom"+nid,800,400);
      fhChrom->SetStats(0);
      fhChrom->Draw("colz");
      prt_canvasAdd("chroml"+nid,800,400);
      fhChromL->SetStats(0);
      fhChromL->Draw("colz");
    }
      
    { // hp
      auto cdigi = prt_drawDigi(fEvId); //2030
      cdigi->SetName("hp"+nid);
      prt_canvasAdd(cdigi);
    }
      
    { // cherenkov ring
      prt_canvasAdd("ring"+nid,500,500);
      
      fHist4->SetStats(0);
      fHist4->GetXaxis()->SetTitle("#theta_{c}sin(#varphi_{c})");
      fHist4->GetYaxis()->SetTitle("#theta_{c}cos(#varphi_{c})");
      fHist4->SetTitle(Form("Calculated from LUT, #theta = %1.2f#circ", prt_theta));
      fHist4->Draw("colz");
      double x0(0), y0(0);
      FitRing(x0,y0,fAngle[2]);
      TVector3 corr(x0,y0,1-TMath::Sqrt(x0*x0+y0*y0));
      //std::cout<<"Tcorr "<< corr.Theta()*1000<< "  Pcorr "<< corr.Phi() <<std::endl;

      TLegend *leg = new TLegend(0.32,0.42,0.67,0.59);
      leg->SetFillStyle(0); 
      leg->SetBorderSize(0);
      leg->AddEntry((TObject*)0,Form("Entries %0.0f",fHist4->GetEntries()),"");
      leg->AddEntry((TObject*)0,Form("#Delta#theta_{c} %f [mrad]",corr.Theta()*1000),"");
      leg->AddEntry((TObject*)0,Form("#Delta#varphi_{c} %f [rad]",corr.Phi()),"");
      leg->Draw();

      TArc *arc = new TArc(x0,y0,fAngle[2]);
      arc->SetLineColor(kRed);
      arc->SetLineWidth(1);
      arc->SetFillStyle(0);
      arc->Draw();
      gg_i=0;
      gg_gr.Set(0);
    }

    { // corrections
      if(fabs(fCorr[0])<0.00000001 && fabs(fCorr[10])<0.00000001 && fabs(fCorr[20])<0.00000001){
	std::cout<<"Writing "<<fCorrFile<<std::endl;
	  
        TFile fc(fCorrFile,"recreate");
        TTree *tc = new TTree("corr","corr");
        int pmt;
        double corr;
        tc->Branch("pmt",&pmt,"pmt/I");
        tc->Branch("corr",&corr,"corr/D");
	tc->Branch("cspr",&spr,"cspr[5]/D");
       
        for(int i=0; i<prt_nmcp; i++){  
	  if(fHistMcp[i]->GetEntries()<20) continue;  
	  prt_canvasAdd(Form("r_tangle_%d",i),800,400);
	    
	  corr= -prt_fit(fHistMcp[i],0.008,20,0.01).X();
	  //fHistMcp[i]->Fit("gaus","MQ","",-0.01,0.01);//.X();
	  //auto f = fHistMcp[i]->GetFunction("gaus");
	  //if(f)
	  {
	    pmt = i;
	    //corr= -f->GetParameter(1);
	    tc->Fill();
	    std::cout<<"if(mcpid=="<< i<<") tangle += "<<corr<<";" <<std::endl;  
	    fHistMcp[i]->Draw();
	  }
        }
	  
        tc->Write();
        fc.Write();
        fc.Close();
      }
    }
      
    { // time
      prt_canvasAdd("tdiff"+nid,800,400);
      prt_normalize(fHistDiff, 3);
      for(int i=0; i<3; i++){
	if(fHistDiff[i]->GetEntries()>100){
	  fHistDiff[i]->SetStats(0);
	  fHistDiff[i]->SetTitle(Form("theta %1.2f", prt_theta));
	  fHistDiff[i]->Draw((i==0)?"h":"hsame");
	}
      }
      
      prt_canvasAdd("diff"+nid,800,400);
      fDiff->SetStats(0);
      fDiff->Draw("colz");
    }
    
    /*TString filedir = fCorrFile;
    if(filedir.Contains("/")) {
      filedir.Remove(filedir.Last('/'));
      filedir += "/";
    }
    else filedir = ""; 
    */

    //prt_canvasSave(filedir+"reco",0,0,0);
    
    prt_canvasSave(filedir + "/reco_plots/" + Form("%1.2f_deg",prt_theta), 0,0,0);

    if(fVerbose>2) prt_waitPrimitive("lh"+nid,"none");
    
  }
}


  void FindPeak(double (&cangle)[5], double (&spr)[5]){
  for(int h=0; h<5; h++){    
    spr[h]=0;
    cangle[h]=0;
      
    if(hthetac[h]->GetEntries()>20 ){
      gROOT->SetBatch(1);
      int nfound = fSpect->Search(hthetac[h],1,"",0.9); //0.6
      if(nfound>0) cangle[h] = fSpect->GetPositionX()[0];
      else cangle[h] =  hthetac[h]->GetXaxis()->GetBinCenter(hthetac[h]->GetMaximumBin());

      fFit->SetParameters(100,cangle[h],0.005,10);   // peak
      fFit->SetParameter(2,0.005); // width
      fFit->FixParameter(2,0.008); // width
      hthetac[h]->Fit("fgaus","Q","",cangle[h]-3.5*fSigma,cangle[h]+3.5*fSigma);
      fFit->ReleaseParameter(2); // width
      hthetac[h]->Fit("fgaus","MQ","",cangle[h]-3.5*fSigma,cangle[h]+3.5*fSigma);
      cangle[h] = fFit->GetParameter(1);
      spr[h] = fFit->GetParameter(2)*1000; 
      if(fVerbose>2) gROOT->SetBatch(0);
    }
  }
  }

void circleFcn(int &, double *, double &f, double *par, int) {
  int np = gg_gr.GetN();
  f = 0;
  double *x = gg_gr.GetX();
  double *y = gg_gr.GetY();
  for (int i=0;i<np;i++) {
    double u = x[i] + par[0];
    double v = y[i] + par[1];
    double dr = par[2] - TMath::Sqrt(u*u+v*v);
    f += dr*dr;  
  }
}

void circleFcn2(int &, double *, double &f, double *par, int) {
  int np = gg_gr.GetN();
  f = 0;
  double *x = gg_gr.GetX();
  double *y = gg_gr.GetY();
  for (int i=0;i<np;i++) {
    double u = x[i] + par[0];
    double v = y[i] + par[1];
    double dr = par[2] - TMath::Sqrt(u*u+v*v);
    if(dr>0.07) f += dr*dr; 
    else f += fabs(dr);
  }
}

void FitRing(double& x0, double& y0, double& theta){

  TGraph ff_gr;
  int ff_i(0);
  int np = gg_gr.GetN();
  double *x = gg_gr.GetX();
  double *y = gg_gr.GetY();
  for (int i=0;i<np;i++) {
    if( fabs(theta - TMath::Sqrt(x[i]*x[i]+y[i]*y[i]))<0.05) {
      ff_gr.SetPoint(ff_i,x[i],y[i]);
      ff_i++;
    }
  }
  gg_gr = ff_gr;

  //Fit a circle to the graph points
  TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
  fitter->SetPrecision(0.00000001);
  fitter->SetMaxIterations(1000);
  double arglist[] = {-1};
  fitter->ExecuteCommand("SET PRINT", arglist,1);
  
  fitter->SetFCN(circleFcn);
  fitter->SetParameter(0, "x0",   0.03, 0.01, -0.05,0.05);
  fitter->SetParameter(1, "y0",   0, 0.01, -0.05,0.05);
  fitter->SetParameter(2, "R",    theta, 0.01, theta-0.05,theta+0.05);

  //fitter->FixParameter(0);
  //fitter->FixParameter(1);

  fitter->FixParameter(2);
  fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  // fitter->SetFCN(circleFcn2);
  // fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  x0 = fitter->GetParameter(0);
  y0 = fitter->GetParameter(1);
  theta = fitter->GetParameter(2);
}

int FindPdg(double mom, double cangle){
  double tdiff, diff=100;
  int minid=0;
  for(int i=0; i<5; i++){
    tdiff = fabs(cangle - acos(sqrt(mom*mom + prt_mass[i]*prt_mass[i])/mom/1.46907)); //1.46907 - fused silica
    if(tdiff<diff){
      diff = tdiff;
      minid = i;
    }
  }
  return prt_pdg[minid]; 
}
	
void drawTheoryLines(double mom){
  gPad->Update();
  
  for(int i=0; i<5; i++){
    fAngle[i] = acos(sqrt(mom*mom + prt_mass[i]*prt_mass[i])/mom/1.4738); //1.4738 = 370 = 3.35
  }
  
  TLine *line = new TLine(0,0,0,1000);
  line->SetX1(fAngle[fp2]);
  line->SetX2(fAngle[fp2]);
  line->SetY1(gPad->GetUymin());
  line->SetY2(gPad->GetUymax());
  line->SetLineColor(prt_color[fp2]);
  line->Draw();

  TLine *line1 = new TLine(0,0,0,1000);
  line1->SetX1(fAngle[fp1]);
  line1->SetX2(fAngle[fp1]);
  line1->SetY1(gPad->GetUymin());
  line1->SetY2(gPad->GetUymax());
  line1->SetLineColor(prt_color[fp1]);
  line1->Draw();
}

      
      

  
		   
			   

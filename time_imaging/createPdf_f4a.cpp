#include "../prttools/prttools.cpp"
//#include "/work/eic3/users/nwickjlb/hpDIRC_reco_test/fun4all_hpDIRC_reco/prttools/prttools.cpp"

#include <iostream>

const int nch(24*256);
TH1F *hlef[nch], *hles[nch];

void createPdf_f4a(TString in="G4DIRCTree.root", int pid=321)
{
  for(int i=0; i<nch; i++){
    hlef[i] = new TH1F(Form("lef_%d",i),";LE time [ns]; entries/N [#]", 2000,0,100);
    hles[i] = new TH1F(Form("les_%d",i),";LE time [ns]; entries/N [#]", 2000,0,100);
  }

  prt_ch = new TChain("mG4EvtTree");
  prt_ch->Add(in);
  int nEvents = prt_ch->GetEntries();
  std::cout << "Entries = " << nEvents << std::endl;

  const int arr_size = 500;

  int hit_size = 0;
  int Particle_id = 0;

  prt_ch->SetBranchAddress("nhits", &hit_size);
  prt_ch->SetBranchAddress("pid", &Particle_id);

  int mcp_num[arr_size], pixel_id[arr_size];
  Double_t lead_time[arr_size];
  
  prt_ch->SetBranchAddress("mcp_id", &mcp_num);
  prt_ch->SetBranchAddress("pixel_id", &pixel_id);
  prt_ch->SetBranchAddress("lead_time", &lead_time);
  
  int printstep = 2000;
  double time = 0;
 
  int pdg(0), totalf(0),totals(0), ch(0);
 
 for (int ievent=0; ievent < nEvents; ievent++)
   { 
     prt_ch->GetEntry(ievent);
     pdg = Particle_id;
  
     //if((pdg == 211 && ievent < 500) || (pdg == 321 && ievent < 1500)) continue;
     if((pdg == 211 && ievent < 5000) || (pdg == 321 && ievent < 25000)) continue; 

     int nHits = hit_size;
     if(ievent%printstep==0 && ievent!=0) std::cout << "Event # " << ievent << " # hits "<< nHits << std::endl;


    if(nHits < 5) continue;
    for(int i=0; i < nHits; i++)
      {
	ch = mcp_num[i]*256 + pixel_id[i];
	time = lead_time[i] + gRandom->Gaus(0,0.1);
      
	if(pdg==pid){
	  totalf++;
	  hlef[ch]->Fill(time);
	}
	if(pdg==211){
	  totals++;
	  hles[ch]->Fill(time);
	}
      }
   }

 std::cout<<"#1 "<< totalf <<"  #2 "<<totals <<std::endl;
  
 if(totalf>=0 && totals>0) 
   {
    in.ReplaceAll(".root",".pdf.root");
    TFile efile(in,"RECREATE");
    
    for(int i=0; i<nch; i++){
      hlef[i]->Scale(1/(double)totalf);
      hles[i]->Scale(1/(double)totals);
      
      hlef[i]->SetName(Form("hf_%d",i));
      hles[i]->SetName(Form("hs_%d",i));
      hlef[i]->Write();
      hles[i]->Write();
      
      if(0){
	int nrow=4, ncol=6, p = i/256;
	int np =p%ncol*nrow + p/ncol;
	hles[i]->SetName(Form("mcp%dpix%d",np,i%256));
	prt_canvasAdd(Form("pdf_mcp%dpix%d",np,i%256),800,400);
      	prt_normalize(hlef[i],hles[i]);
      	hlef[i]->SetLineColor(2);
	hles[i]->SetLineColor(4);
	hles[i]->GetXaxis()->SetRangeUser(5, 20);
      	hles[i]->Draw("hist");	
      	hlef[i]->Draw("hist same");
	prt_canvasSave("data/pdfs_7",0,0,1);
      }
    }
    
    efile.Write();
    efile.Close();
   }

}

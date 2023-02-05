#include <iostream>
//#include "/work/eic3/users/nwickjlb/hpDIRC_reco_test/fun4all_hpDIRC_reco/prttools/prttools.cpp"
#include "../prttools/prttools_32x8.cpp"

void draw_hp(TString infile)
{
  TChain *fChain = new TChain("mG4EvtTree");
  fChain->Add(infile);

  int nEvents = fChain->GetEntries();

  TGaxis::SetMaxDigits(3);
  prt_setRootPalette(1);
  prt_initDigi(2);

  TH1F* h_pip_nph = new TH1F("pip_nph",";multiplicity [#]; entries [#]", 200,0,200);

  const int arr_size = 500;

  int hit_size = 0;

  fChain->SetBranchAddress("nhits", &hit_size);

  int mcp_num[arr_size], pixel_id[arr_size];

  fChain->SetBranchAddress("mcp_id", &mcp_num);
  fChain->SetBranchAddress("pixel_id", &pixel_id);

  double mean_nph_pip, mean_nph_pip_err, sigma_nph_pip, sigma_nph_pip_err;
  double nhits_cut;

  for (int ievent=0; ievent < nEvents/2; ievent++)
    {
      fChain->GetEntry(ievent);
    
      int numHits = hit_size;
      int nph_pip = 0;
      
      for(int i=0; i < numHits; i++) 
	{
	  nph_pip++;
	}
      
      h_pip_nph->Fill(nph_pip);
    }

  TFitResultPtr ptr = h_pip_nph->Fit("gaus","SQ","",0,200);
  mean_nph_pip = ptr->Parameter(1);
  mean_nph_pip_err = ptr->ParError(1);
  sigma_nph_pip = ptr->Parameter(2);
  sigma_nph_pip_err = ptr->ParError(2);
    
  nhits_cut = mean_nph_pip - 3*sigma_nph_pip;
  double nhits_cut_val = nhits_cut;

  std::cout << "number of hits cut at " << nhits_cut << std::endl;

  for (int ievent=0; ievent < nEvents/2; ievent++){
  //for (int ievent=0; ievent < 2000; ievent++){
    fChain->GetEntry(ievent);
    int nHits = hit_size;  
    if(nHits < nhits_cut) continue;
    if(ievent%1000==0) std::cout<<"Event # "<< ievent << " has "<< nHits <<" hits"<<std::endl;
    
    for(int h=0; h < nHits; h++)
      {      
	int mcp = mcp_num[h];
	int pix = pixel_id[h];

	prt_hdigi[mcp]->Fill(pix%8, pix/8);
      }
  }

  { // hp
    auto cdigi = prt_drawDigi(2030); 
    cdigi->SetName("hp");
    prt_canvasAdd(cdigi);
  }

  prt_canvasSave(".", 0,0,0);

}



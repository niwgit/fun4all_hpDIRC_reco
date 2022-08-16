#include "PrtLutNode.h"
#include <iostream>
#include "../prttools/prttools.cpp"

void lutmean(TString baseFile = "../data/lut"){  

  gROOT->ProcessLine(".L PrtLutNode.cxx+");

  TString inFile =baseFile+".root";
  TString outFile =baseFile+"_avr.root";

  TFile* f = new TFile(inFile);
  TTree *t=(TTree *) f->Get("prtlut") ;
  TClonesArray* fLut[10];
  for(Int_t l=0; l < 10; l++){
    fLut[l] = new TClonesArray("PrtLutNode");
    t->SetBranchAddress(Form("LUT_%d",l),&fLut[l]);
  }
  
  TFile *fFileNew = TFile::Open(outFile, "RECREATE");
  TClonesArray *fLutNew[10];

  TTree *fTreeNew = new TTree("prtlut","Look-up table for DIRC. Averaged");

  Int_t Nnodes = 24*256;

  for(Int_t l=0; l < 10; l++){
    fLutNew[l] = new TClonesArray("PrtLutNode");
    fTreeNew->Branch(Form("LUT_%d",l),&fLutNew[l],256000,2);
    TClonesArray &fLutaNew = *fLutNew[l];
    for (Long64_t n=0; n<Nnodes; n++) {
      new((fLutaNew)[n]) PrtLutNode(-1);
    }
  }
  
  TCanvas* c = new TCanvas("c","c",0,0,800,1200); c->Divide(1,2);
  TH1F * histNode = new TH1F("LutNode","Node vs Multiplicity",30000,0,150000);
  TH1F * hTime = new TH1F("hTime","Time",5000,0,10);
  TH1F * hDir = new TH1F("hDir","X component",1000,-1,1);


  std::vector<TVector3> vArray[1000];
  std::vector<Double_t> tArray[1000];
  std::vector<Long_t> pArray;
  std::vector<int> rArray;
  
  TVector3 dir, dir2, sum;
  Double_t angle,minangle,time,sumt;
  Long_t pathid;
  int nrefl;
  PrtLutNode *node;

  for(Int_t l=0; l< 10; l++){
    t->GetEntry(l);
    for (Int_t inode=0; inode<fLut[l]->GetEntriesFast(); inode++){
      if(inode%1000==0) cout<<"Node # "<<inode << "  L "<<l<<endl;
      node= (PrtLutNode*) fLut[l]->At(inode);
      Int_t size = node->Entries();    
      if(size<1) continue;
      
      for(int i=0; i<size; i++){
	dir = node->GetEntry(i).Unit();
	time = node->GetTime(i);
	//if(node->GetNRefl(i) > 10) continue; 
	pathid = node->GetPathId(i);
	nrefl = node->GetNRefl(i);
      
	bool newid = true;
	for(int j=0; j<pArray.size(); j++){
	  if(pathid==pArray[j]){
	    vArray[j].push_back(dir);
	    tArray[j].push_back(time);
	    newid= false;
	  }
	}
	if(newid) {
	  vArray[pArray.size()].push_back(dir);
	  tArray[pArray.size()].push_back(time);
	  pArray.push_back(pathid);
	  rArray.push_back(nrefl);
	}
      }
    
      for(int j=0; j<pArray.size(); j++){
	sum = TVector3(0,0,0);
	sumt=0;
	hDir->Reset();
	hTime->Reset();
      
	for(int v=0; v<vArray[j].size(); v++) {
	  sum += vArray[j][v]; 
	  sumt += tArray[j][v]; 

	  hDir->Fill(vArray[j][v].Y());
	  hTime->Fill(tArray[j][v]);
	}

	if(vArray[j].size()<2) continue;
	if(rArray[j]>15) continue;
      
	if(hDir->GetStdDev()>0.02)
	  {
	    //continue;
	    std::cout<<inode<<" "<<pArray[j]<<" hDir->GetStdDev() ================  "<<hDir->GetStdDev()<<std::endl;
	    std::cout << "varray size = " << vArray[j].size() << "\t" << "tarray size = " << tArray[j].size() << std::endl;

	    c->cd(1);
	    hTime->Draw();
	    c->cd(2);
	    hDir->Draw();
	    c->Update();
	    c->WaitPrimitive();
	  }
      
	sum *= 1/(Double_t)vArray[j].size();
	sumt *= 1./(Double_t)tArray[j].size();

      
	((PrtLutNode*)(fLutNew[l]->At(inode)))->AddEntry(inode, sum,pArray[j],rArray[j],sumt, node->GetDigiPos(),node->GetDigiPos());
      
      }

      for(int i=0; i<1000; i++) 
	{
	  vArray[i].clear();
	  tArray[i].clear();
	}
      pArray.clear();
      rArray.clear();
    }
  }

  fTreeNew->Fill();
  fTreeNew->Write();
  fFileNew->Write();

}

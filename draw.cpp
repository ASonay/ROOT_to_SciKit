#include "TH1D.h"

#include "inc/MVAExtra.hh"

using namespace std;

double GetIntegral(TH1 *h){

  double sum=0;
  double dx = (double)(h->GetXaxis()->GetXmax()-h->GetXaxis()->GetXmin())/(double)h->GetNbinsX();
  
  for (int i=0;i<h->GetNbinsX();i++)
    sum+=h->GetBinContent(i+1)*dx;

  return sum;

}

void draw()
{
  gStyle->SetOptStat(0);
  gStyle->SetHatchesLineWidth(3);
  //TGaxis::SetMaxDigits(2);

  TCut Sig = "(class_id)*weight";
  TCut Bckg = "(!class_id)*weight";
  
  string file = "output_1l_10ji3bi_0.root";
  //string file = "higgs.root";
  

  //double ne=500;double xmin=0;double xmax=5800; string selection = "HT_all/1e3";
  double ne=100;double xmin=-0.2;double xmax=1.2; string selection = "score";
  
  
  TFile *f = new TFile(file.c_str());
  TTree *test_tr = (TTree*)f->Get("TestTree");
  TTree *train_tr = (TTree*)f->Get("TrainTree");
 
  TH1D *S_t = new TH1D("S_t","Signal",ne,xmin,xmax);
  TH1D *B_t = new TH1D("B_t","Background",ne,xmin,xmax);
  TH1D *S_tr = new TH1D("S_tr","Signal",ne,xmin,xmax);
  TH1D *B_tr = new TH1D("B_tr","Background",ne,xmin,xmax);
  test_tr->Project("S_t",selection.c_str(),Sig);
  test_tr->Project("B_t",selection.c_str(),Bckg);
  train_tr->Project("S_tr",selection.c_str(),Sig);
  train_tr->Project("B_tr",selection.c_str(),Bckg);
  double sfac = ne/(xmax-xmin);
  S_t->Scale(sfac/S_t->Integral());
  B_t->Scale(sfac/B_t->Integral());
  //Sval_t->Scale(3*sfac/Sval_t->Integral());
  //Bval_t->Scale(3*sfac/Bval_t->Integral());
  S_tr->Scale(sfac/S_tr->Integral());
  B_tr->Scale(sfac/B_tr->Integral());



  
  TH2F *hcan = new TH2F("hcan","",100,S_t->GetXaxis()->GetXmin(),S_t->GetXaxis()->GetXmax(),100,0,10.95);
  hcan->GetXaxis()->SetTitle("BDTG Response");
  //hcan->GetXaxis()->SetTitle("H_{T} (GeV)");
  hcan->GetYaxis()->SetTitle("1/N (dN/dx)");


  TCanvas *c = new TCanvas("c");
  c->SetTicky(); c->SetTickx(); 

  S_t->SetFillStyle(1001); S_t->SetFillColor(38); S_t->SetLineColor(4); S_t->SetLineWidth(3);
  B_t->SetFillStyle(3445); B_t->SetFillColor(2); B_t->SetLineColor(2); B_t->SetLineWidth(3);
  //Sval_t->SetLineColor(kYellow-7); Sval_t->SetLineWidth(3);
  //Bval_t->SetLineColor(kGreen-7); Bval_t->SetLineWidth(3);
  S_tr->SetMarkerStyle(20); S_tr->SetMarkerSize(1.7); S_tr->SetMarkerColor(4); S_tr->SetLineColor(4); S_tr->SetLineWidth(3);
  B_tr->SetMarkerStyle(20); B_tr->SetMarkerSize(1.7); B_tr->SetMarkerColor(2); B_tr->SetLineColor(2); B_tr->SetLineWidth(3);


  hcan->Draw();
  S_t->Draw("hist same");
  B_t->Draw("hist same");
  S_tr->Draw("e1 same");
  B_tr->Draw("e1 same");
  //Sval_t->Draw("hist same");
  //Bval_t->Draw("hist same");

  string namex="10ji3bi";
  TLegend *l_1 = new TLegend(0.0740375,0.776119,0.368213,0.951493);
  l_1->SetLineWidth(3);
  l_1->SetFillStyle(0);
  l_1->AddEntry(S_t,("Signal Test @"+namex).c_str(),"f");
  l_1->AddEntry(B_t,("Background Test @"+namex).c_str(),"f");
  l_1->Draw();
  TLegend *l_2 = new TLegend(0.368213,0.776119,0.682132,0.951493);
  l_2->SetLineWidth(3);
  l_2->SetFillStyle(0);
  l_2->AddEntry(S_tr,("Signal Train @"+namex).c_str(),"p");
  l_2->AddEntry(B_tr,("Background Train @"+namex).c_str(),"p");
  l_2->Draw();
  TLegend *l_3 = new TLegend(0.682132,0.776119,0.995064,0.951493);
  l_3->SetLineWidth(3);
  l_3->SetFillStyle(0);
  //l_3->AddEntry(Sval_t,("Signal Eval @"+namex).c_str(),"l");
  //l_3->AddEntry(Bval_t,("Background Eval @"+namex).c_str(),"l");
  l_3->Draw();

  //-----------------------------------------------------------

  MVAExtra MVAE("MVAE",S_t,B_t);
  MVAExtra MVAE_tr("MVAE_tr",S_tr,B_tr);
  //MVAExtra MVAE_val("MVAE_val",Sval_t,Bval_t);
  
  cout << "Seperation Eval: " << MVAE.GetSeparation()
       << " Seperation Trained: " << MVAE_tr.GetSeparation() << endl;

  //-----------------------------------------------------------

  TCanvas *cg = new TCanvas("cg");
  cg->SetTicky(); cg->SetTickx();

  TH1D *rc_1 = MVAE.GetROC(2);
  TH1D *rc_3 = MVAE_tr.GetROC(2);
  //TH1D *rc_2 = MVAE_val.GetROC(2);

  rc_1->SetLineColor(kGray+3); rc_1->SetLineWidth(5); rc_1->SetLineStyle(1);
  //rc_2->SetLineColor(kRed-4); rc_2->SetLineWidth(5); rc_2->SetLineStyle(9);
  rc_3->SetLineColor(kAzure-3); rc_3->SetLineWidth(5); rc_3->SetLineStyle(1);
  rc_3->SetMarkerColor(kGray+3); rc_3->SetMarkerStyle(20); rc_3->SetMarkerSize(1.6);

  rc_1->Draw("L hist");
  //rc_2->Draw("L hist same");
  rc_3->Draw("hist L same");
  
  TLegend *l_4 = new TLegend(0.11846,0.205224,0.82231,0.376866);
  l_4->SetLineWidth(3);
  l_4->SetFillStyle(0);
  l_4->SetHeader(("ROC Curves @"+namex).c_str());
  l_4->AddEntry(rc_1,"Test Sample","l");
  //l_4->AddEntry(rc_2,"Evaluation Sample","l");
  l_4->AddEntry(rc_3,"Training Sample","l");
  l_4->Draw();

}

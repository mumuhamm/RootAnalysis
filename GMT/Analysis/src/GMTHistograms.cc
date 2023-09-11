#include <iostream>
#include <cmath>

#include "GMTHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TRandom3.h"
#include "TGaxis.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooJohnson.h"
#include "TLatex.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TColor.h"
#include "RooHist.h"
#include "RooCrystalBall.h"
#include "RooBreitWigner.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
using namespace std;
using namespace RooFit ;


#include "utilsL1RpcStyle.h"

int nPtBins = 32;
const float GMTHistograms::ptBins[36]={0., 0.1,
  		 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8.,
		 10., 12., 14., 16., 18., 20., 22, 24, 26, 30., 35., 40., 45.,
  		 50., 60., 70., 80., 90., 100., 120., 140.,
		 160., 200};

const int GMTHistograms::color[6] = {kBlack, kBlue, kRed, kMagenta, kTeal, kGreen};
//Single mu
const int GMTHistograms::ptCutsOMTF[4] =   {0, 14, 19, 20};

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
GMTHistograms::GMTHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
GMTHistograms::GMTHistograms(TDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
GMTHistograms::GMTHistograms(TDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
GMTHistograms::~GMTHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string GMTHistograms::getTemplateName(const std::string& name){

  std::string templateName = "";

  if(name.find("DeltaEta")!=std::string::npos) templateName = "h1DDeltaEtaTemplate";
  if(name.find("DeltaPhi")!=std::string::npos) templateName = "h1DDeltaPhiTemplate";
  if(name.find("Pt")!=std::string::npos) templateName = "h2DPtTemplate";
  if(name.find("HighPt")!=std::string::npos) templateName = "h2DHighPtTemplate";
  if(name.find("PtTag")!=std::string::npos) templateName = "h1DPtTagTemplate";
  if(name.find("AbsEtaTag")!=std::string::npos) templateName = "h1DAbsEtaTagTemplate";
  if(name.find("PtProbe")!=std::string::npos) templateName = "h1DPtProbeTemplate";
  if(name.find("AbsEtaProbe")!=std::string::npos) templateName = "h1DAbsEtaProbeTemplate";
  if(name.find("DiMuonMassTagProbe")!=std::string::npos) templateName = "h1DDiMuonMassTagProbeTemplate";
   
  if(name.find("Eta")!=std::string::npos) templateName = "h2DEtaTemplate";
  if(name.find("Phi")!=std::string::npos) templateName = "h2DPhiTemplate";
  if(name.find("RecoMuonPtVsL1Pt")!=std::string::npos) templateName = "h2DRecoMuonPtVsL1PtTemplate";

  if(name.find("Quality")!=std::string::npos) templateName = "h2DQualityTemplate";
  if(name.find("OMTFRateTot")!=std::string::npos) templateName = "h2DOMTFRateTotTemplate";
  if(name.find("uGMTRateTot")!=std::string::npos) templateName = "h2DuGMTRateTotTemplate";
  if(name.find("RateVsEta")!=std::string::npos) templateName = "h2DRateVsEtaTemplate";
  if(name.find("RateVsPt")!=std::string::npos) templateName = "h2DRateVsPtTemplate";
  if(name.find("RateVsQuality")!=std::string::npos) templateName = "h2DRateVsQualityTemplate";
  
  if(name.find("LLH")!=std::string::npos) templateName = "h1DLLHTemplate";
  if(name.find("HitsPattern")!=std::string::npos) templateName = "h1DHitsPatternTemplate";
  
  return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void GMTHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){

 //Make template histos
 add1DHistogram("h1DDeltaEtaTemplate","",11,-0.83,0.83,file_);
 add1DHistogram("h1DDeltaPhiTemplate","",5*32,-M_PI,M_PI,file_);
 add1DHistogram("h1DPtTagTemplate", "", 100, 0, 100, file_);
 add1DHistogram("h1DAbsEtaTagTemplate", "", 60, -2.4, 2.4, file_);
 add1DHistogram("h1DPtProbeTemplate", "", 100, 0, 100, file_);
 add1DHistogram("h1DAbsEtaProbeTemplate", "", 60, 0.8, 1.4, file_);
 add1DHistogram("h1DDiMuonMassTagProbeTemplate", "", 100, 70, 110, file_); 



 ///Efficiency histos
 add2DHistogram("h2DPtTemplate","",150,0,150,2,-0.5,1.5,file_);
 add2DHistogram("h2DHighPtTemplate","",50,50,550,2,-0.5,1.5,file_);
 add2DHistogram("h2DRecoMuonPtVsL1PtTemplate", "", 100, 0, 120, 100, 0, 120, file_);


 add2DHistogram("h2DEtaTemplate","",60,-3,3,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiTemplate","",4*32,-3.2,3.2,2,-0.5,1.5,file_);
 add2DHistogram("h2DQualityTemplate","",201,-0.5,200.5,2,-0.5,1.5,file_);

 //Rate histos
 add2DHistogram("h2DOMTFRateTotTemplate","",404,1,202,404,1,202,file_);
 add2DHistogram("h2DuGMTRateTotTemplate","",404,1,202,404,1,202,file_);
 add2DHistogram("h2DRateVsEtaTemplate","",404,1,202,60,-3,3,file_);
 add2DHistogram("h2DRateVsPtTemplate","",404,1,202,100,0,50,file_);
 add2DHistogram("h2DRateVsQualityTemplate","",404,1,202,201,-0.5,200.5,file_);

 //Likelihood histos
 add1DHistogram("h1DLLHTemplate","",40,0,20,file_);
 add1DHistogram("h1DHitsPatternTemplate","",101,-0.5,100.5,file_);

 
 histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void GMTHistograms::finalizeHistograms(){

  AnalysisHistograms::finalizeHistograms();
  utilsL1RpcStyle()->cd();
  gErrorIgnoreLevel = kError;
  
  plotEffPanel("OMTF");// Many turn on curves 
  plotEffPanel("OMTF", true);// High pT turn on curves 
 
  //1D histograms 
  
  plotSingleHistogram("h1DPtTag");
  plotSingleHistogram("h1DAbsEtaTag");
  plotSingleHistogram("h1DPtProbe");
  plotSingleHistogram("h1DAbsEtaProbe");
  plotSingleHistogram("h1DDiMuonMassTagProbe"); 
  plotSingleHistogram("h2DOMTFRecoMuonPtVsL1Pt");
  plotSingleHistogram("h2DOMTFRateTot");
  plotSingleHistogram("h2DuGMTRateTot");
 
  //Efficiency as a function of eta and phis, Lines for selected points on the turn on curve shown
  
  plotEffVsEta("OMTF");
  plotEffVsVar("OMTF", "Eta");
  plotEffVsVar("OMTF", "Phi");
    
  //Turn on curves for many pT thresholds.
  ///Lines for reference - Phase2 uGMT, and other algorithm shown
  for(int iPtCode=1;iPtCode<=30;++iPtCode){
      plotGMTVsOther(iPtCode,"OMTF");
  }

  plotRate("Tot");

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1* GMTHistograms::Integrate(TH1 * histoD) {

  TH1* histoI = (TH1*)histoD->Clone("hIntegrated");

  Double_t * cont = new Double_t [histoD->GetNbinsX()+2];  //with under+overflow
  Double_t * errs = new Double_t [histoD->GetNbinsX()+2];  //with under+overflow
  histoI->Reset();

  // bin=0 underf
  // bin 1-GetNbinsX() -conten
  // bin GetNbinsX()+1 overflow

  for (Int_t i = 0; i <= histoD->GetNbinsX()+1; i++) {
    cont[i] = (Double_t)(histoD->GetBinContent(i));
    errs[i] = (Double_t)(histoD->GetBinError(i));
  }
  Double_t sum=0.;
  Double_t sume2=0.;
  for (Int_t i = histoD->GetNbinsX()+1; i > 0; i--) {
    sum+=cont[i];
    sume2+=errs[i]*errs[i];
    histoI->SetBinContent(i,sum);
    histoI->SetBinError(i,sqrt(sume2));
   }
  delete [] cont;
  delete [] errs;
  return histoI;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TEfficiency * GMTHistograms::DivideErr(TH1D * h1, TH1D * h2,const char * name,const char * optErr){


      if (!h1) std::cout <<"DivideErr called, but histogram (h1) pointer is:"<<h1<<std::endl;
      if (!h2) std::cout <<"DivideErr called, but histogram (h2) pointer is:"<<h1<<std::endl;
      if (!h1 || !h2) return 0;
      TEfficiency * hout = 0;
      if(TEfficiency::CheckConsistency(*h1,*h2)){
                         hout = new TEfficiency(*h1,*h2);
        }

       hout->SetName(name);
       return hout;

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void GMTHistograms::DrawLabels(TCanvas* c){//, const TString& eraLabel) {
    TLatex* cmsLabel = new TLatex();
    cmsLabel->SetTextFont(42);
    cmsLabel->SetTextSize(0.04);
    cmsLabel->SetTextAlign(11); // Left-align
    cmsLabel->DrawLatexNDC(0.17, 0.92, "#bf{CMS} #it{Preliminary}");

    TLatex* lumiLabel = new TLatex();
    lumiLabel->SetTextFont(42);
    lumiLabel->SetTextSize(0.04);
    lumiLabel->SetTextAlign(31); // Right-align
    TString lumiText =  "20.1fb^{-1}(13.6 TeV)";//"DrellYan";// now just era we have lumi info though eraLabel;
    lumiLabel->DrawLatexNDC(0.94444, 0.92, lumiText);

    c->Update();
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void GMTHistograms::plotEffPanel(const std::string & sysType, bool doHigh){

  TCanvas* c = new TCanvas(TString::Format("EffVsPt_%s",sysType.c_str()),
			   TString::Format("EffVsPt_%s",sysType.c_str()),
			   600,600);

  TLegend l(0.6513158,0.1673729,0.8903509,0.470339,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(1,1);

  TString hName("");
  const int *ptCuts = ptCutsOMTF;
 
  for (int icut=0; icut <=3;++icut){
    float ptCut = GMTHistograms::ptBins[ptCuts[icut]];
    hName = "h2D"+sysType+"Pt"+std::to_string((int)ptCut);
    if(doHigh) hName = "h2D"+sysType+"HighPt"+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName.Data());
    if(!h2D) return;
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);    
    hDenom->Add(hNum);
    TEfficiency* hEff =DivideErr(hNum,hDenom,"Pt_Int","B"); //TH1D
    hEff->SetMarkerStyle(21+icut);
    hEff->SetMarkerColor(color[icut]);
    hEff->SetTitle("; p_{T}^{reco} (GeV/c);L1 Muon efficiency");
    if (icut==0)hEff->Draw();
    else hEff->Draw("same");
    TString nameCut = TString::Format("%d", (int)GMTHistograms::ptBins[ptCuts[icut]])+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.Data());
    c->Update();
    auto graph = hEff->GetPaintedGraph();
    graph->GetXaxis()->SetRangeUser(0.0,60.0);
    graph->GetYaxis()->SetRangeUser(0.0,1.1);
    c->Update();
    DrawLabels(c); 
    c->Update();
  }
  l.DrawClone();
  
  if(!doHigh){ 
               c->Print(TString::Format("fig_png/PanelVsPt_%s.png",sysType.c_str()).Data());
               c->Print(TString::Format("fig_png/PanelVsPt_%s.pdf",sysType.c_str()).Data());
             }
  else{  
        c->Print(TString::Format("fig_png/PanelVsHighPt_%s.png",sysType.c_str()).Data());
	c->Print(TString::Format("fig_png/PanelVsHighPt_%s.pdf",sysType.c_str()).Data());
      }
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotVar(const std::string & sysType,
			    const std::string & varName){

  TCanvas* c = new TCanvas(TString::Format("Var%s_%s",varName.c_str(),sysType.c_str()),
			   TString::Format("Var%s_%s",varName.c_str(),sysType.c_str()),
			   460,500);


  TString hName = "h1D"+varName+sysType;
  TH1F* h1D = this->get1DHistogram(hName.Data());
  h1D->SetLineWidth(3);
  h1D->SetStats(kFALSE);

  c->SetLogy();

  h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
  h1D->SetXTitle("#Delta R(RECO - RECO)");
  h1D->Draw();

  std::cout<<sysType<<" integral total: "<<h1D->Integral(0,h1D->GetNbinsX()+1)<<std::endl;
  std::cout<<sysType<<" integral 0 - 1.0: "<<h1D->Integral(0,h1D->FindBin(1.0))<<std::endl;
  std::cout<<sysType<<" integral 0 - 0.5: "<<h1D->Integral(0,h1D->FindBin(0.5))<<std::endl;
  std::cout<<sysType<<" integral 0 - 0.05: "<<h1D->Integral(0,h1D->FindBin(0.05))<<std::endl;
  std::cout<<sysType<<" integral 0 - 0.02: "<<h1D->Integral(0,h1D->FindBin(0.02))<<std::endl;
 DrawLabels(c);
 c->Print(TString::Format("fig_eps/Var%s_%s.eps",varName.c_str(), sysType.c_str()).Data());
 c->Print(TString::Format("fig_png/Var%s_%s.png",varName.c_str(), sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotEffVsVar(const std::string & sysType,
		const std::string & varName){

  TCanvas* c = new TCanvas(TString::Format("EffVs%s_%s",varName.c_str(),sysType.c_str()),
			   TString::Format("EffVs%s_%s",varName.c_str(),sysType.c_str()),
			   460,500);

  TLegend l(0.6513158,0.1673729,0.8903509,0.470339,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(1,1);

  TString hName("");
  const int *ptCuts = ptCutsOMTF;
 
  for (int icut=0; icut<2;++icut){
    float ptCut = GMTHistograms::ptBins[ptCuts[icut]];
    hName = "h2D"+sysType+varName+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName.Data());
    if(!h2D) return;
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    hDenom->Add(hNum);
    TEfficiency* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
    hEff->SetMarkerStyle(21+icut);
    hEff->SetMarkerColor(color[icut]); // Please fix the title accriding to histograms 
    if(varName.find("Phi")!=std::string::npos)hEff->SetTitle(";|#phi^{#mu}| (rad); L1 Muon Efficiency");
    if(varName.find("Eta")!=std::string::npos)hEff->SetTitle(";|#eta^{#mu}| (a.u.); L1 Muon Efficiency");
    if (icut==0)hEff->Draw();
    else hEff->Draw("same");
    TString nameCut = TString::Format("%d", (int)GMTHistograms::ptBins[ptCuts[icut]])+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.Data());
    c->Update();
    auto graph = hEff->GetPaintedGraph();
    graph->GetXaxis()->SetRangeUser(0.0,100.0);
    graph->GetYaxis()->SetRangeUser(0.0,1.0);
    c->Update();

  }
  l.DrawClone();
  DrawLabels(c);
  c->Print(TString::Format("fig_eps/EffVs%s_%s.eps",varName.c_str(), sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVs%s_%s.png",varName.c_str(), sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVs%s_%s.pdf",varName.c_str(), sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotEffVsEta(const std::string & sysType){

  TCanvas* c = new TCanvas(TString::Format("EffVsEta_%s",sysType.c_str()),
			   TString::Format("EffVsEta_%s",sysType.c_str()),
			   600,600);
  
  TLegend l(0.6915995,0.5930233,0.7422325,0.8972868,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(1,1);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.35);

  int iCut = 18;
  std::string hName = "";
  for (int iType=0; iType<3;++iType){
    float ptCut = GMTHistograms::ptBins[iCut];
    hName = "h2D"+sysType+"Type" + std::to_string(iType) + "EtauGMT"+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName);
    if(!h2D) return;
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    if(iType==2) hNum->Scale(50.0);
    hDenom->Add(hNum);
    TEfficiency* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
    hEff->SetMarkerStyle(21+iCut);
    hEff->SetMarkerColor(color[iCut]);
    hEff->SetTitle(";#eta^{reco} (a.u.); L1 Muon-Efficiency");
    if (iType==0)hEff->Draw();
    else hEff->Draw("same");
    std::string nameCut = std::to_string((int)GMTHistograms::ptBins[iCut])+" GeV/c";
    if (iType==0) nameCut = "p_{T}^{#mu}>p_{T}^{cut} + 20 GeV/c";
    if (iType==1) nameCut = "p_{T}^{cut}<p_{T}^{#mu}<#dot p_{T}^{cut} + 5 GeV/c";
    if (iType==2) nameCut = "p_{T}^{#mu}<10 GeV/c (#epsilon #times 50)";
    l.AddEntry(hEff,nameCut.c_str());

    c->Update();
    auto graph = hEff->GetPaintedGraph();
    graph->GetXaxis()->SetRangeUser(0.0,100.0);
    graph->GetYaxis()->SetRangeUser(0.0,1.0);
    c->Update();
  }
  ///OMTF eta range used for generating patterns.
  TLine *aLine = new TLine(0,0,0,0);
  aLine->SetLineWidth(2);
  aLine->SetLineColor(2);
  aLine->DrawLine(0.83,0,0.83,1.0);
  aLine->DrawLine(-0.83,0,-0.83,1.0);
  aLine->DrawLine(1.24,0,1.24,1.0);
  aLine->DrawLine(-1.24,0,-1.24,1.0);

  l.SetHeader(TString::Format("p_{T} = %d  GeV/c",(int)GMTHistograms::ptBins[iCut]).Data());
  l.DrawClone();
  DrawLabels(c);
  c->Print(TString::Format("fig_png/EffVsEta_%s.png",sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVsEta_%s.pdf",sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotGMTVsOther(int iPtCut,
				     const std::string sysType){

  float ptCut = ptBins[iPtCut];

  TCanvas* c = new TCanvas(TString::Format("OMTFVsOther_%d",(int)ptCut).Data(),
			   TString::Format("OMTFVsOther_%d",(int)ptCut).Data(),
			   600,600);

  TLegend l(0.2,0.65,0.44,0.86,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetLogx(1);
  c->SetGrid(1,1);

  std::string hName = "h2D"+sysType+"Pt"+std::to_string((int)ptCut);
  TH2F* h2D = get2DHistogram(hName);
  if(!h2D) return;
  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);
  TEfficiency* hEffOther =DivideErr(hNum,hDenom,"hEffOther","B");
  hEffOther->SetMarkerStyle(23);
  hEffOther->SetMarkerColor(2);
  


  hName = "h2DuGMTPt"+std::to_string((int)ptCut);
  h2D = get2DHistogram(hName);
  hNum = h2D->ProjectionX("hNum",2,2);
  hDenom = h2D->ProjectionX("hDenom",1,1);    
  hDenom->Add(hNum);

  TEfficiency* hEffGMT =DivideErr(hNum,hDenom,"hEffGMTTmp","B");
  hEffGMT->SetMarkerStyle(8);
  hEffGMT->SetMarkerColor(1);
  hEffGMT->SetTitle(";p_{T}^{Reco} (GeV/c);L1 muon efficiency");
  hEffGMT->Draw(); 
  hEffOther->Draw("same");

  c->Update();
  auto graph = hEffGMT->GetPaintedGraph();
  graph->GetXaxis()->SetRangeUser(0.0,100.0);
  graph->GetYaxis()->SetRangeUser(0.0,1.0);
  c->Update();

  c->Update();
  auto graph2 = hEffOther->GetPaintedGraph();
  graph2->GetXaxis()->SetRangeUser(0.0,100.0);
  graph2->GetYaxis()->SetRangeUser(0.0,1.0);
  c->Update();

  std::string tmp = "p_{T} #geq ";
  if(int(ptCut*10)%10==5) tmp += "%1.1f GeV/c";
  else   tmp += "%1.0f GeV/c";
  l.AddEntry((TObject*)0, TString::Format(tmp.c_str(),ptCut).Data(), "");
  l.AddEntry((TObject*)0, "", "");
  l.AddEntry(hEffOther, sysType.c_str());
  l.AddEntry((TObject*)0, "", "");
  l.AddEntry(hEffGMT, "Phase1GMT");
  l.DrawClone();

  TLine aLine(0,0,0,0);
  aLine.SetLineColor(2);
  aLine.SetLineWidth(3);
  aLine.DrawLine(ptCut,0,ptCut,1.04);
  DrawLabels(c);
  c->Print(TString::Format("fig_eps/uGMTVs%s_%d.eps",sysType.c_str(),(int)ptCut).Data());
  c->Print(TString::Format("fig_png/uGMTVs%s_%d.png",sysType.c_str(),(int)ptCut).Data());

  c->SetLogy();
  c->Print(TString::Format("fig_eps/uGMTVs%s_%d_log.eps",sysType.c_str(),(int)ptCut).Data());
  c->Print(TString::Format("fig_png/uGMTVs%s_%d_log.png",sysType.c_str(),(int)ptCut).Data());
  c->Print(TString::Format("fig_png/uGMTVs%s_%d_log.pdf",sysType.c_str(),(int)ptCut).Data());

}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TH2F* GMTHistograms::makeRateWeights(TH2 *hOrig){

  TF1 *fIntVxMuRate = new TF1("fIntVxMuRate","0.1*TMath::Power(x,[0]*TMath::Log(x))*TMath::Power(x,[1])*TMath::Exp([2])",1,1000);
  fIntVxMuRate->SetParameters(-0.235801, -2.82346, 17.162);

  TH2F *hWeights = (TH2F*)hOrig->Clone("hWeights");
  hWeights->Reset();
  TH1D *hPtGen = hOrig->ProjectionX("hPtGen");

  int nEvInBin;
  float ptLow, ptHigh, weight;

  for (int iBin = 1; iBin <= hPtGen->GetNbinsX(); ++iBin){
    ptLow = hPtGen->GetXaxis()->GetBinLowEdge(iBin);
    ptHigh = hPtGen->GetXaxis()->GetBinUpEdge(iBin);
    nEvInBin = hPtGen->GetBinContent(iBin);
    if(nEvInBin<1) nEvInBin = 1;
    weight = (fIntVxMuRate->Eval(ptLow) - fIntVxMuRate->Eval(ptHigh))/nEvInBin;
    for (int iBinY = 0; iBinY<=hOrig->GetNbinsY()+1;++iBinY) hWeights->SetBinContent(iBin,iBinY,weight);
  }
  return hWeights;
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TH1* GMTHistograms::getRateHisto(std::string sysType,
				 std::string type){

  std::string hName = "h2D"+sysType+"Rate"+type;
  std::cout << " the name of the histogram , if exists : "<< hName << "\n";
  TH2F* h2D_original = (TH2F*)this->get2DHistogram(hName);
  if(! h2D_original) return 0;
  TH2F* h2D = (TH2F*)h2D_original->Clone("h2D");
  if(!h2D) return 0;

  if(selectionFlavours_.size() &&
     selectionFlavours_[0].find("NU_RATE")==std::string::npos){
    TH2F *hWeights = makeRateWeights(h2D);
    h2D->Multiply(hWeights);
  }
  
  TH1D *hRate = h2D->ProjectionY(("hRate"+sysType).c_str());
  if(sysType=="OMTF" || sysType=="uGMT") hRate = h2D->ProjectionX("hRate");
  std::cout<< " the name after the projection : "<< hRate->GetName()<< "\n";
  hRate->SetYTitle("Arbitrary units");
  hRate->SetLineWidth(3);
  if(type.find("Tot")!=std::string::npos) return(TH1*)hRate->Clone("hRateClone");
  if(type.find("VsEta")!=std::string::npos) return (TH1*)hRate->Clone("hRateClone");
  if(type.find("VsPt")!=std::string::npos) return (TH1*)hRate->Clone("hRateClone");
  if(type.find("VsQuality")!=std::string::npos) return (TH1*)hRate->Clone("hRateClone");

  return Integrate(hRate);
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotRate(std::string type){
  
  TH1 *hRateuGMT = getRateHisto("uGMT",type);
  TH1 *hRateOMTF = getRateHisto("OMTF",type);
  std::cout<< " the name of the histrogram : when plot rate :: printing OMTF : "<<  hRateOMTF->GetName() << "\n";
  if(!hRateOMTF || !hRateuGMT) return;
  std::cout<< " again printing the OMTF : "<< hRateOMTF->GetName()<< "\n";
  hRateOMTF->SetLineWidth(3); 
  hRateOMTF->SetLineColor(1);
  
  hRateuGMT->SetLineWidth(3);
  hRateuGMT->SetLineColor(2);
  hRateuGMT->SetLineStyle(2);
 
  TCanvas* c = new TCanvas("cRate","cRate",1.5*420,1.5*500);
  c->SetLogy(1);
  c->SetGrid(1,1);

  TLegend *leg = new TLegend(0.60,0.75,0.85,0.9,NULL,"brNDC");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  if(type.find("Tot")!=std::string::npos){
    hRateOMTF->SetAxisRange(2,50);
    hRateuGMT->SetAxisRange(2,50);
    hRateOMTF->SetMinimum(1E1);
    hRateOMTF->SetMaximum(2E5);
    hRateOMTF->SetXTitle("p_{T}^{cut} (GeV/c)");
    c->SetLogy();
    c->SetGrid(1,0);
    hRateOMTF->Draw();
    hRateuGMT->DrawCopy("same");
   
    std::cout<<"Rate OMTF @ 20 GeV: "<< hRateOMTF->GetBinContent(hRateOMTF->FindBin(20-0.01))<<std::endl;
  }  
  else if(type.find("VsEta")!=std::string::npos){
    c->SetLogy(0);
    hRateuGMT->SetXTitle("muon #eta");
    double max = hRateuGMT->GetMaximum();
    hRateuGMT->SetMaximum(1.5*max);
    hRateuGMT->Draw();
  }
  if(type=="VsPt"){
    c->SetLogy(1);
    hRateuGMT->SetXTitle("p_{T}^{Reco} (GeV/c)");
    hRateuGMT->SetAxisRange(2,100);
    double max = hRateuGMT->GetMaximum();
    hRateuGMT->SetMaximum(10*max);
    hRateuGMT->Draw();
  }

  leg->AddEntry(hRateuGMT,"OMTF");
  leg->Draw();
  DrawLabels(c);
  c->Print(("fig_eps/Rate"+type+".eps").c_str());
  c->Print(("fig_png/Rate"+type+".png").c_str());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
/*float  GMTHistograms::getEfficiency(TH2F *h2D, float ptCut){

  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);
  TEfficiency* hEffTmp =DivideErr(hNum,hDenom,"hEffTmp","B");
  
  float range = hEffTmp->GetBinLowEdge(binHigh+1) - hEffTmp->GetBinLowEdge(binLow);
  float eff = hEffTmp->Integral(binLow,binHigh,"width")/range;
  return eff;
}*/
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TEfficiency * GMTHistograms::getEfficiencyHisto(const std::string & hName){

  TH2F* h2D = this->get2DHistogram(hName);
  if(!h2D) return 0;
  
  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);  
  hDenom->Add(hNum);
  TEfficiency* hEffTmp =DivideErr(hNum,hDenom,"hEffTmp","B");
  return hEffTmp;
}
////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
void GMTHistograms::plotSingleHistogram(std::string hName){

  TH2F* h2D = get2DHistogram(hName);
  TH1F* h1D = get1DHistogram(hName);
  if(!h2D && !h1D) return;
	
  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
			   600,600);

  TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  if(h2D) {
    h2D->SetDirectory(myDirCopy);
    h2D->SetLineWidth(3);
    h2D->Scale(1.0/h2D->Integral());
    h2D->SetXTitle("p_{T}^{Reco} (GeV/c)");
    h2D->SetYTitle("p_{T}^{L1Matched} (GeV/c)");
    h2D->GetYaxis()->SetTitleOffset(1.4);
    h2D->SetStats(kFALSE);
    gStyle->SetPalette(kRainBow);
    h2D->Draw("box colz");
    DrawLabels(c);
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
    c->Print(TString::Format("fig_png/%s.pdf",hName.c_str()).Data());
  }
  if(h1D) {
    h1D->SetDirectory(myDirCopy);
    h1D->SetLineWidth(2);
    TH1D *diMuonClone = (TH1D*)(h1D->Clone("diMuonClone"));
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));    
    h1D->GetXaxis()->SetRange(1,h1D->GetNbinsX()+1);
    const std::string  histName = hName.c_str();
    std::cout<< " the name of the hist string : "<< histName << "\n";
    if(histName.find("PtTag")!=std::string::npos)h1D->SetXTitle("Tag p_{T}^{#mu RECO} (GeV/c)");
    if(histName.find("AbsEtaTag")!=std::string::npos)h1D->SetXTitle("Tag |#eta^{#mu RECO}| (a.u.)");
    if(histName.find("PtProbe")!=std::string::npos)h1D->SetXTitle("Probe p_{T}^{#mu RECO} (GeV/c)");
    if(histName.find("AbsEtaProbe")!=std::string::npos)h1D->SetXTitle("Probe |#eta^{#mu RECO}| (a.u.)");
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    gStyle->SetOptStat(0);
    h1D->Draw("");
    DrawLabels(c);
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
    c->Print(TString::Format("fig_png/%s.pdf",hName.c_str()).Data());
    
    if(histName.find("DiMuonMassTagProbe")!=std::string::npos){
    RooRealVar mass("mass", "m_{Z}(#mu^{+}#mu^{-}) (GeV/c^{2})", 80, 100);
    RooDataHist dh("dh", "dh", mass, Import(*diMuonClone));
    RooRealVar mean("mean", "mean", 91.18, 90, 93);
    RooRealVar width("width", "width", 5.5, 0, 20);
    RooBreitWigner brwg("brwg", "brwg", mass, mean, width);
    
    RooRealVar a0("a0", "a0", 0.1, -0.2, 0.3);
    RooRealVar a1("a1", "a1", 0.2, 0.0, 1.);
    RooChebychev pol("pol", "pol", mass, RooArgList(a0));
    

    RooRealVar mu("mu", "mu", 91.18, 90, 93);
    RooRealVar lambda ("lambda", "lambda",  0.5, 0.0, 10);
    RooRealVar gamma ("gamma", "gamma",  0., 0.0, 10);
    RooRealVar delta ("delta", "delta", 1., 0.0, 20);
    RooJohnson john ("john", "john", mass, mu, lambda, gamma, delta);
  
    RooRealVar conts ("conts","conts",  -3.0, 0.0);
    RooExponential expoBg ("expoBg","expoBG",mass,conts);

    RooRealVar nSig("nSig", "nSig",500,0,(int)h1D->GetEntries());
    RooRealVar nBkg("nBkg", "nBkg",400,0,(int)h1D->GetEntries());
    RooRealVar nlsb("nlsb", "nlsb",200,0,(int)h1D->GetEntries());
    RooAddPdf zpdf("zpdf","zpdf",RooArgList(john,expoBg), RooArgList(nSig,nBkg));
    //RooAddPdf zpdf("zpdf","zpdf",RooArgList(brwg,expoBg,pol), RooArgList(nSig,nBkg,nlsb));
    RooFitResult* fitRes = zpdf.fitTo(dh,Save(true),NumCPU(8), RooFit::Minimizer("Minuit2", "Migrad"));
    fitRes->Print("v");
    
    TCanvas* cDimuon = new TCanvas("cDimuon","the dimuon distribition", 800,800);
    RooPlot* zmassf = mass.frame(Title("m_{Z}(#mu^{+}#mu^{-}) (GeV/c^{2})"));
    dh.plotOn(zmassf,DataError(RooAbsData::SumW2), MarkerSize(0.7), LineColor(12));
    zpdf.plotOn(zmassf,RooFit::LineColor(kBlue),RooFit::Components("zpdf"), RooFit::Name("Total"), LineWidth(2), LineStyle(9) ) ;

    RooPlot* pullframezmass = mass.frame(RooFit::Title("Zmass"));
    RooHist* hpullmass = zmassf->pullHist();
    pullframezmass->addPlotable(hpullmass,"P0");
    pullframezmass->SetMinimum(-3) ;
    pullframezmass->SetMaximum(+3) ;
    pullframezmass->GetYaxis()->SetTitle("pull");
    pullframezmass->GetYaxis()->SetTitleSize(0.2);
    pullframezmass->GetYaxis()->SetTitleOffset(0.3);
    pullframezmass->GetXaxis()->SetTitle("m_{Z}(#mu^{+}#mu^{-}) (GeV/c^{2})");
    pullframezmass->GetXaxis()->SetTitleSize(0.2);
    pullframezmass->GetXaxis()->SetTitleOffset(0.85);
    pullframezmass->SetMarkerStyle(20);
    pullframezmass->SetNdivisions(10);
    double chisquare_mass = zmassf->chiSquare();
    std::cout<<"Chi square of mass fit is :   "<< chisquare_mass<<"\n";

    zpdf.plotOn(zmassf, RooFit::LineColor(kGreen),RooFit::Components("john"), RooFit::Name("signal"), LineWidth(2), LineStyle(4));
    zpdf.plotOn(zmassf,RooFit::LineColor(kRed),RooFit::Components("expoBg"), RooFit::Name("combinatorial"), LineWidth(2), LineStyle(6));
    //zpdf.plotOn(zmassf,RooFit::LineColor(kGreen),RooFit::Components("pol"), RooFit::Name("lowerSB"), LineWidth(2), LineStyle(8)); 
  
    TLegend *leg = new TLegend(0.15,0.55,0.45,0.85, NULL, "brNDC");
    leg->SetTextSize(0.05);
    leg->SetFillStyle(4000);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->AddEntry(zmassf->findObject("Total"),"Total PDF","l");
    leg->AddEntry(zmassf->findObject("signal"),"Z^{0}#rightarrow #mu^{+}#mu^{-}","l");
    leg->AddEntry(zmassf->findObject("combinatorial"),"Combinatorial","l");
    //leg->AddEntry(zmassf->findObject("lowerSB"),"lowerSB","l");
    leg->SetBorderSize(0);
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,0.97);
    TPad *pad2 = new TPad("pad2","pad2",0.0,0.1,1,0.27);
    pad1->SetBottomMargin(0.05);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.08);  
    pad2->SetBottomMargin(0.45);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    gStyle->SetOptTitle(0);
    cDimuon->SetFillColor(0);
    cDimuon->SetBorderSize(2);
    cDimuon->SetLeftMargin(0.14);
    cDimuon->SetRightMargin(0.044);
    cDimuon->SetBottomMargin(0.001);
    zmassf->SetStats(0);
    gStyle->SetTitle(0);
    zmassf->SetTitle("");
    zmassf->Draw();
    leg->Draw("same");
    pad2->cd();
    pullframezmass->SetStats(0);
    pullframezmass->SetLabelSize(0.2, "X");
    pullframezmass->SetLabelSize(0.1, "Y");
    pullframezmass->Draw();
    cDimuon->cd();
    DrawLabels(cDimuon);   
    cDimuon->Print("fig_png/diMuonMassTagandProbe.png", "png");
    cDimuon->Print("fig_png/diMuonMassTagandProbe.pdf", "pdf");
    }
  }
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>

bool makePlots_png  = true;
bool makePlots_pdf  = false;
bool makePlots_root = false;

TFile* openFile(const std::string& inputFilePath, const std::string& inputFileName)
{
  TString inputFileName_full = inputFilePath.data();
  if ( !inputFileName_full.EndsWith("/") ) inputFileName_full.Append("/");
  inputFileName_full.Append(inputFileName.data());
  TFile* inputFile = new TFile(inputFileName_full.Data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = " << inputFileName_full.Data() << " !!" << std::endl;
    assert(0);
  }
  return inputFile;
}

TTree* loadTree(TFile* inputFile, const std::string& directory, const std::string& treeName)
{  
  TString treeName_full = directory.data();
  if ( !treeName_full.EndsWith("/") ) treeName_full.Append("/");
  treeName_full.Append(treeName.data());
  TTree* tree = (TTree*)inputFile->Get(treeName_full.Data());
  if ( !tree ) {
    std::cerr << "Failed to load tree = " << treeName_full.Data() << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  std::cout << "Successfully loaded tree = '" << treeName_full.Data() << "' from file = " << inputFile->GetName() << std::endl;
  std::cout << " #entries = " << tree->GetEntries() << std::endl;
  return tree;
}

struct histogramEntryType
{
  histogramEntryType()
  {
    histogram_memLR_    = new TH1D("memLR",    "memLR",    50,   0.,    1.);
    histogram_memLRerr_ = new TH1D("memLRerr", "memLRerr", 50,   0.,    1.);
    histogram_memProbS_ = new TH1D("memProbS", "memProbS", 40, -70.,  +10.);
    histogram_memProbB_ = new TH1D("memProbB", "memProbB", 40, -70.,  +10.);
    histogram_drbb_     = new TH1D("drbb",     "drbb",     50,   0.,    5.);
    histogram_mbb_      = new TH1D("mbb",      "mbb",      40,   0.,  200.);
    histogram_drll_     = new TH1D("drll",     "drll",     50,   0.,    5.);
    histogram_dphill_   = new TH1D("dphill",   "dphill",   36,   0., TMath::Pi());
    histogram_mll_      = new TH1D("mll",      "mll",      40,   0.,  200.); 
  }
  ~histogramEntryType()
  {
    delete histogram_memLR_;
    delete histogram_memLRerr_;
    delete histogram_memProbS_;
    delete histogram_memProbB_;
    delete histogram_drbb_;
    delete histogram_mbb_;
    delete histogram_drll_;
    delete histogram_dphill_;
    delete histogram_mll_;
  }
  TH1* histogram_memLR_;
  TH1* histogram_memLRerr_;
  TH1* histogram_memProbS_;
  TH1* histogram_memProbB_;
  TH1* histogram_drbb_;
  TH1* histogram_mbb_;
  TH1* histogram_drll_;
  TH1* histogram_dphill_;
  TH1* histogram_mll_;
};

TH1* getHistogram(histogramEntryType* histograms, const std::string& histogramName)
{
  if      ( histogramName == "memLR"    ) return histograms->histogram_memLR_;
  else if ( histogramName == "memLRerr" ) return histograms->histogram_memLRerr_;
  else if ( histogramName == "memProbS" ) return histograms->histogram_memProbS_;
  else if ( histogramName == "memProbB" ) return histograms->histogram_memProbB_;
  else if ( histogramName == "drbb"     ) return histograms->histogram_drbb_;
  else if ( histogramName == "mbb"      ) return histograms->histogram_mbb_;
  else if ( histogramName == "drll"     ) return histograms->histogram_drll_;
  else if ( histogramName == "dphill"   ) return histograms->histogram_dphill_;
  else if ( histogramName == "mll"      ) return histograms->histogram_mll_;
  else {
    std::cerr << "Invalid histogram name = " << histogramName << " !!" << std::endl;
    assert(0);
  }
  return nullptr;
}

void fillWithOverFlow(TH1* histogram, double x, double evtWeight)
{
  if ( !histogram ) return;
  const TAxis* const xAxis = histogram->GetXaxis();
  int idxBin = xAxis->FindBin(x);
  if ( idxBin < 1                 ) idxBin = 1;
  if ( idxBin > xAxis->GetNbins() ) idxBin = xAxis->GetNbins();
  double binCenter = histogram->GetBinCenter(idxBin);
  histogram->Fill(binCenter, evtWeight);
}

void fillWithOverFlow_logx(TH1* histogram, double x, double evtWeight)
{    
  const double nonzero = 1.e-30;
  fillWithOverFlow(histogram, TMath::Log(TMath::Max(nonzero, x)), evtWeight);
}

histogramEntryType*
fillHistograms(TTree* tree, const std::string& selection)
{
  int numEntries = tree->GetEntries();

  TTreeFormula treeFormula("treeFormula", selection.data(), tree);
  int numEntries_selected = 0;
  for ( int idxEntry = 0; idxEntry < numEntries; ++idxEntry ) {
    tree->GetEntry(idxEntry);
    if ( !treeFormula.EvalInstance() ) continue;
    ++numEntries_selected;
  }

  std::cout << "Applying selection= '" << selection << "'" << std::endl;
  std::cout << " " << numEntries_selected << " out of " << numEntries << " entries selected." << std::endl;

  Double_t memLR;
  tree->SetBranchAddress("memLR", &memLR);
  Double_t memLRerr;
  tree->SetBranchAddress("memLRerr", &memLRerr);
  Double_t memProbS;
  tree->SetBranchAddress("memProbS", &memProbS);
  Double_t memProbB;
  tree->SetBranchAddress("memProbB", &memProbB);
  Float_t drbb;
  tree->SetBranchAddress("drbb", &drbb);
  Float_t mbb;
  tree->SetBranchAddress("mbb", &mbb);
  Float_t drll;
  tree->SetBranchAddress("drll", &drll);
  Float_t dphill;
  tree->SetBranchAddress("dphill", &dphill);
  Float_t mll;
  tree->SetBranchAddress("mll", &mll);

  histogramEntryType* histograms = new histogramEntryType();

  for ( int idxEntry = 0; idxEntry < numEntries; ++idxEntry ) {
    tree->GetEntry(idxEntry);
    if ( !treeFormula.EvalInstance() ) continue;
    
    const double evtWeight = 1.;

    if ( (idxEntry % 100) == 0 ) { 
      std::cout << "entry #" << idxEntry << ":" 
                << " memLR = " << memLR << " (memProbS = " << memProbS << ", memProbB = " << memProbB << ")" << std::endl;
    }

    fillWithOverFlow(histograms->histogram_memLR_,         memLR,              evtWeight);
    fillWithOverFlow(histograms->histogram_memLRerr_,      memLRerr,           evtWeight);
    fillWithOverFlow_logx(histograms->histogram_memProbS_, memProbS,           evtWeight);
    fillWithOverFlow_logx(histograms->histogram_memProbB_, memProbB,           evtWeight);
    fillWithOverFlow(histograms->histogram_drbb_,          drbb,               evtWeight);
    fillWithOverFlow(histograms->histogram_mbb_,           mbb,                evtWeight);
    fillWithOverFlow(histograms->histogram_drll_,          drll,               evtWeight);
    fillWithOverFlow(histograms->histogram_dphill_,        TMath::Abs(dphill), evtWeight);
    fillWithOverFlow(histograms->histogram_mll_,           mll,                evtWeight);
  }

  return histograms;
}

void showHistograms(double canvasSizeX, double canvasSizeY,
		    TH1* histogram1, const std::string& legendEntry1,
		    TH1* histogram2, const std::string& legendEntry2,
		    TH1* histogram3, const std::string& legendEntry3,
		    TH1* histogram4, const std::string& legendEntry4,
		    int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
		    double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions, 
		    const std::string& labelText, double labelTextSize,
		    double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		    const std::string& xAxisTitle, double xAxisOffset,
		    bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		    const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = ( outputFileName.find("_prob") != std::string::npos ) ? 0.162 : 0.150;
  const double margin_bottom   = ( outputFileName.find("_prob") != std::string::npos ) ? 0.140 : 0.125;
  const double margin_right    = 0.015;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(margin_top);
  canvas->SetLeftMargin(margin_left);
  canvas->SetBottomMargin(margin_bottom);
  canvas->SetRightMargin(margin_right); 
  canvas->SetLogx(false);
  canvas->SetLogy(useLogScale);
  canvas->Draw();
  canvas->cd();
  
  assert(histogram1);
  histogram1->SetFillColor(0);
  histogram1->SetFillStyle(0);
  histogram1->SetLineColor(colors[0]);
  histogram1->SetLineStyle(lineStyles[0]);
  histogram1->SetLineWidth(lineWidths[0]);
  histogram1->SetMarkerColor(colors[0]);
  histogram1->SetMarkerStyle(markerStyles[0]);
  histogram1->SetMarkerSize(markerSizes[0]);

  assert(histogram2);
  histogram2->SetFillColor(0);
  histogram2->SetFillStyle(0);
  histogram2->SetLineColor(colors[1]);
  histogram2->SetLineStyle(lineStyles[1]);
  histogram2->SetLineWidth(lineWidths[1]);
  histogram2->SetMarkerColor(colors[1]);
  histogram2->SetMarkerStyle(markerStyles[1]);
  histogram2->SetMarkerSize(markerSizes[1]);

  if ( histogram3 ) {
    histogram3->SetFillColor(0);
    histogram3->SetFillStyle(0);
    histogram3->SetLineColor(colors[2]);
    histogram3->SetLineStyle(lineStyles[2]);
    histogram3->SetLineWidth(lineWidths[2]);
    histogram3->SetMarkerColor(colors[2]);
    histogram3->SetMarkerStyle(markerStyles[2]);
    histogram3->SetMarkerSize(markerSizes[2]);
  }

  if ( histogram4 ) {
    histogram4->SetFillColor(0);
    histogram4->SetFillStyle(0);
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineStyle(lineStyles[3]);
    histogram4->SetLineWidth(lineWidths[3]);
    histogram4->SetMarkerColor(colors[3]);
    histogram4->SetMarkerStyle(markerStyles[3]);
    histogram4->SetMarkerSize(markerSizes[3]);
  }
  
  TAxis* xAxis = histogram1->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(60);
  xAxis->SetTitleFont(43);
  //xAxis->SetLabelOffset(-0.01);
  xAxis->SetLabelSize(0.050);
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = histogram1->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(60);
  yAxis->SetTitleFont(43);
  if ( yMax > yMin ) {
    yAxis->SetRangeUser(yMin, yMax);
  }
  //yAxis->SetLabelOffset(0.010);
  yAxis->SetLabelSize(0.055);
  yAxis->SetLabelFont(42);
  yAxis->SetTickLength(0.040);  
  yAxis->SetNdivisions(505);

  histogram1->SetTitle("");
  histogram1->SetStats(false);

  histogram1->Draw(Form("%ssame", drawOptions[0].data()));
  histogram2->Draw(Form("%ssame", drawOptions[1].data()));
  if ( histogram3 ) histogram3->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogram4 ) histogram4->Draw(Form("%ssame", drawOptions[3].data()));
  histogram1->Draw("axissame");

  TPaveText* label = nullptr;
  if ( labelText != "" ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->AddText(labelText.data());
    label->SetTextFont(42);
    label->SetTextSize(labelTextSize);
    label->SetTextColor(1);
    label->SetTextAlign(13);
    label->Draw();
  }

  TLegend* legend = nullptr;
  if ( legendEntry1 != "" && legendEntry2 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(histogram1, legendEntry1.data(), legendOptions[0].data());
    legend->AddEntry(histogram2, legendEntry2.data(), legendOptions[1].data());
    if ( histogram3 ) legend->AddEntry(histogram3, legendEntry3.data(), legendOptions[2].data());
    if ( histogram4 ) legend->AddEntry(histogram4, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  if ( makePlots_png  ) canvas->Print(std::string(outputFileName_plot).append(".png").data());
  if ( makePlots_pdf  ) canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  if ( makePlots_root ) canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete label;
  delete legend;
  delete canvas;
}

void debug_bbww_dilepton()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);
 
  std::string inputFilePath = "/home/veelken/CMSSW_11_1_2/CMSSW_11_1_2/src/hhAnalysis/bbwwMEMPerformanceStudies/test/DEBUG/";

  std::string inputFileName_signal = "analyze_hh_bbwwMEM_dilepton_signal_ggf_nonresonant_node_sm_hh_2b2v_jetSmearingDisabled_metSmearingDisabled_all.root";
  TFile* inputFile_signal = openFile(inputFilePath, inputFileName_signal);

  std::string inputFileName_background = "analyze_hh_bbwwMEM_dilepton_TTJets_DiLept_ext1_jetSmearingDisabled_metSmearingDisabled_all.root";
  TFile* inputFile_background = openFile(inputFilePath, inputFileName_background);

  std::string directory = "ntuples";
  std::string treeName  = "mem";

  TTree* tree_signal     = loadTree(inputFile_signal,     directory, treeName);
  TTree* tree_background = loadTree(inputFile_background, directory, treeName);

  std::string selection_low_memLR  = "memLR < 0.1";
  std::string selection_high_memLR = "memLR > 0.9";

  histogramEntryType* histograms_signal_low_memLR      = fillHistograms(tree_signal,     selection_low_memLR);
  histogramEntryType* histograms_signal_high_memLR     = fillHistograms(tree_signal,     selection_high_memLR);
  histogramEntryType* histograms_background_low_memLR  = fillHistograms(tree_background, selection_low_memLR);
  histogramEntryType* histograms_background_high_memLR = fillHistograms(tree_background, selection_high_memLR);
  
  std::string labelText_signal = "HH #rightarrow b#bar{b} WW^{*} #rightarrow b#bar{b} l^{+}#nu l^{-}#bar{#nu}";
  std::string labelText_background = "t#bar{t} #rightarrow bW #bar{b}W #rightarrow b l^{+}#nu #bar{b} l^{-}#bar{#nu}";

  int showHistograms_canvasSizeX = 1050;
  int showHistograms_canvasSizeY =  950;
  double showHistograms_xAxisOffset = 0.96;
  double showHistograms_yAxisOffset = 1.21;
  int showHistograms_colors[4]       = { kGreen - 6, kBlack, kBlue - 7, 28 };
  int showHistograms_markerStyles[4] = { 20, 24, 21, 25 };
  int showHistograms_markerSizes[4]  = { 2, 2, 2, 2 };
  int showHistograms_lineStyles[4]   = { 1, 1, 1, 7 };
  int showHistograms_lineWidths[4]   = { 3, 2, 2, 3 };
  std::vector<std::string> showHistograms_drawOptions = { "hist", "ep", "ep", "hist" };
  std::vector<std::string> showHistograms_legendOptions = { "l", "p", "p", "l" };

  std::vector<std::string> histogramNames;
  histogramNames.push_back("memLR");
  histogramNames.push_back("memLRerr");
  histogramNames.push_back("memProbS");
  histogramNames.push_back("memProbB");
  histogramNames.push_back("drbb");
  histogramNames.push_back("mbb");
  histogramNames.push_back("drll");
  histogramNames.push_back("dphill");
  histogramNames.push_back("mll");

  std::map<std::string, std::string> xAxisTitle; // key = histogramName
  std::map<std::string, std::string> yAxisTitle; // key = histogramName
  std::map<std::string, double>      yMin;       // key = histogramName
  std::map<std::string, double>      yMax;       // key = histogramName

  xAxisTitle["memLR"]    = "P";
  yAxisTitle["memLR"]    = "dN/dP";
  yMin["memLR"]          = 1.1e-5;
  yMax["memLR"]          = 1.9e0;

  xAxisTitle["memLRerr"] = "#sigma P";
  yAxisTitle["memLRerr"] = "dN/d#sigma P";
  yMin["memLRerr"]       = 1.1e-5;
  yMax["memLRerr"]       = 1.9e0;

  xAxisTitle["memProbS"] = "log w_{0}";
  yAxisTitle["memProbS"] = "dN/dlog w_{0}";
  yMin["memProbS"]       = 1.1e-5;
  yMax["memProbS"]       = 1.9e0;

  xAxisTitle["memProbB"] = "log w_{1}";
  yAxisTitle["memProbB"] = "dN/dlog w_{1}";
  yMin["memProbB"]       = 1.1e-5;
  yMax["memProbB"]       = 1.9e0;

  xAxisTitle["drbb"]     = "#Delta r_{bb}";
  yAxisTitle["drbb"]     = "dN/d#Delta r_{bb}";
  yMin["drbb"]           = 1.1e-5;
  yMax["drbb"]           = 1.9e0;

  xAxisTitle["mbb"]      = "m_{bb} [GeV]";
  yAxisTitle["mbb"]      = "dN/dm_{bb} [1/GeV]";
  yMin["mbb"]            = 1.1e-5;
  yMax["mbb"]            = 1.9e0;

  xAxisTitle["drll"]     = "#Delta r_{ll}";
  yAxisTitle["drll"]     = "dN/d#Delta r_{ll}";
  yMin["drll"]           = 1.1e-5;
  yMax["drll"]           = 1.9e0;
  
  xAxisTitle["dphill"]   = "#Delta #phi_{ll} [Rad]";
  yAxisTitle["dphill"]   = "dN/d#Delta #phi_{ll} [1/Rad]";
  yMin["dphill"]         = 1.1e-5;
  yMax["dphill"]         = 1.9e0;

  xAxisTitle["mll"]      = "m_{ll} [GeV]";
  yAxisTitle["mll"]      = "dN/dm_{ll} [1/GeV]";
  yMin["mll"]            = 1.1e-5;
  yMax["mll"]            = 1.9e0;

  for ( std::vector<std::string>::const_iterator histogramName = histogramNames.begin();
        histogramName != histogramNames.end(); ++histogramName ) {
    TH1* histogram_signal_low_memLR      = getHistogram(histograms_signal_low_memLR,      *histogramName);
    histogram_signal_low_memLR->Scale(1./histogram_signal_low_memLR->Integral());
    TH1* histogram_signal_high_memLR     = getHistogram(histograms_signal_high_memLR,     *histogramName);
    histogram_signal_high_memLR->Scale(1./histogram_signal_high_memLR->Integral());
    
    TH1* histogram_background_low_memLR  = getHistogram(histograms_background_low_memLR,  *histogramName);
    histogram_background_low_memLR->Scale(1./histogram_background_low_memLR->Integral());
    TH1* histogram_background_high_memLR = getHistogram(histograms_background_high_memLR, *histogramName);
    histogram_background_high_memLR->Scale(1./histogram_background_high_memLR->Integral());

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_signal_low_memLR,  "P < 0.10",
      histogram_signal_high_memLR, "P > 0.90",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_legendOptions,
      "", 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      xAxisTitle[*histogramName], showHistograms_xAxisOffset,
      true, yMin[*histogramName], yMax[*histogramName], yAxisTitle[*histogramName], showHistograms_yAxisOffset, 
      Form("debug_bbww_dilepton_signal_%s_unsmeared.pdf", histogramName->data()));

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_background_low_memLR,  "P < 0.10",
      histogram_background_high_memLR, "P > 0.90",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_legendOptions,
      "", 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      xAxisTitle[*histogramName], showHistograms_xAxisOffset,
      true, yMin[*histogramName], yMax[*histogramName], yAxisTitle[*histogramName], showHistograms_yAxisOffset, 
      Form("debug_bbww_dilepton_background_%s_unsmeared.pdf", histogramName->data()));
  }

  delete inputFile_signal;
  delete inputFile_background;
}

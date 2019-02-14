
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>

enum { kUndefined, kSignal, kBackground };

enum { kProbSignal, kProbBackground, kLR };

TH1* loadHistogram(TFile* inputFile, const std::string& directory_part1, const std::string& directory_part2, int signal_or_background, const std::string& histogramName)
{  
  TString histogramName_full = directory_part1.data();
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  histogramName_full.Append(directory_part2.data());
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  if      ( signal_or_background == kSignal     ) histogramName_full.Append("signal_ggf_nonresonant_node_sm_hh_bbvv");
  else if ( signal_or_background == kBackground ) histogramName_full.Append("TT");
  else assert(0);
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  histogramName_full.Append(histogramName.data());
  TH1* histogram = (TH1*)inputFile->Get(histogramName_full.Data());
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName_full.Data() << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  double integral = histogram->Integral(1, histogram->GetNbinsX()); // CV: exclude underflow and overflow bins
  if ( integral > 0. ) histogram->Scale(1./integral);
  return histogram;
}

TH1* addHistograms(const std::string& histogramSumName, const TH1* histogram1, const TH1* histogram2, const TH1* histogram3)
{
  TH1* histogramSum = (TH1*)histogram1->Clone(histogramSumName.data());
  histogramSum->Reset();
  if ( !histogramSum->GetSumw2N() ) histogramSum->Sumw2();
  histogramSum->Add(histogram1);
  histogramSum->Add(histogram2);
  if ( histogram3 ) histogramSum->Add(histogram2);
  double integral = histogramSum->Integral(1, histogramSum->GetNbinsX()); // CV: exclude underflow and overflow bins
  if ( integral > 0. ) histogramSum->Scale(1./integral);
  return histogramSum;
}

TH1* rebinHistogram(const TH1* histogram, int numBinsX)
{
  TH1* histogram_rebinned = (TH1*)histogram->Clone(Form("%s_rebinned", histogram->GetName()));
  if ( !histogram_rebinned->GetSumw2N() ) histogram_rebinned->Sumw2();
  if ( numBinsX < histogram->GetNbinsX() && (histogram->GetNbinsX() % numBinsX) == 0 ) {
    histogram_rebinned->Rebin(histogram->GetNbinsX() / numBinsX);
  }
  return histogram_rebinned;
}

void showHistograms(double canvasSizeX, double canvasSizeY,
		    TH1* histogram1, const std::string& legendEntry1,
		    TH1* histogram2, const std::string& legendEntry2,
		    TH1* histogram3, const std::string& legendEntry3,
		    TH1* histogram4, const std::string& legendEntry4,
		    int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[],
		    double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
		    const std::string& labelText, double labelTextSize,
		    double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		    int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		    bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		    const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.065);
  canvas->SetLeftMargin(0.17);
  canvas->SetBottomMargin(0.155);
  canvas->SetRightMargin(0.015); 
  canvas->SetLogx(false);
  canvas->SetLogy(useLogScale);
  canvas->Draw();
  canvas->cd();

  assert(histogram1);
  TH1* histogram1_rebinned = rebinHistogram(histogram1, numBinsX);
  histogram1_rebinned->SetFillColor(0);
  histogram1_rebinned->SetFillStyle(0);
  histogram1_rebinned->SetLineColor(colors[0]);
  histogram1_rebinned->SetLineStyle(lineStyles[0]);
  histogram1_rebinned->SetLineWidth(lineWidths[0]);
  histogram1_rebinned->SetMarkerColor(colors[0]);
  //histogram1_rebinned->SetMarkerStyle(markerStyles[0]);
  //histogram1_rebinned->SetMarkerSize(markerSizes[0]);

  assert(histogram2);
  TH1* histogram2_rebinned = rebinHistogram(histogram2, numBinsX);
  histogram2_rebinned->SetFillColor(0);
  histogram2_rebinned->SetFillStyle(0);
  histogram2_rebinned->SetLineColor(colors[1]);
  histogram2_rebinned->SetLineStyle(lineStyles[1]);
  histogram2_rebinned->SetLineWidth(lineWidths[1]);
  histogram2_rebinned->SetMarkerColor(colors[1]);
  histogram2_rebinned->SetMarkerStyle(markerStyles[1]);
  histogram2_rebinned->SetMarkerSize(markerSizes[1]);

  TH1* histogram3_rebinned = nullptr;     
  if ( histogram3 ) {
    histogram3_rebinned = rebinHistogram(histogram3, numBinsX);
    histogram3_rebinned->SetFillColor(0);
    histogram3_rebinned->SetFillStyle(0);
    histogram3_rebinned->SetLineColor(colors[2]);
    histogram3_rebinned->SetLineStyle(lineStyles[2]);
    histogram3_rebinned->SetLineWidth(lineWidths[2]);
    histogram3_rebinned->SetMarkerColor(colors[2]);
    histogram3_rebinned->SetMarkerStyle(markerStyles[2]);
    histogram3_rebinned->SetMarkerSize(markerSizes[2]);
  }

  TH1* histogram4_rebinned = nullptr; 
  if ( histogram4 ) {
    histogram4_rebinned = rebinHistogram(histogram4, numBinsX);
    histogram4_rebinned->SetFillColor(0);
    histogram4_rebinned->SetFillStyle(0);
    histogram4_rebinned->SetLineColor(colors[3]);
    histogram4_rebinned->SetLineStyle(lineStyles[3]);
    histogram4_rebinned->SetLineWidth(lineWidths[3]);
    histogram4_rebinned->SetMarkerColor(colors[3]);
    histogram4_rebinned->SetMarkerStyle(markerStyles[3]);
    histogram4_rebinned->SetMarkerSize(markerSizes[3]);
  }
  
  TAxis* xAxis = histogram1_rebinned->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(60);
  xAxis->SetTitleFont(43);
  //xAxis->SetLabelOffset(-0.01);
  xAxis->SetLabelSize(0.050);
  if ( xMax > xMin ) {
    xAxis->SetRangeUser(xMin, xMax);
  }
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = histogram1_rebinned->GetYaxis();
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

  histogram1_rebinned->SetTitle("");
  histogram1_rebinned->SetStats(false);

  histogram1_rebinned->Draw("hist");
  histogram2_rebinned->Draw("histsame");
  if ( histogram3_rebinned ) histogram3_rebinned->Draw("epsame");
  if ( histogram4_rebinned ) histogram4_rebinned->Draw("epsame");
  histogram1_rebinned->Draw("axissame");

  TPaveText* label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->AddText(labelText.data());
  label->SetTextFont(42);
  label->SetTextSize(labelTextSize);
  label->SetTextColor(1);
  label->SetTextAlign(13);
  label->Draw();

  TLegend* legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
  legend->SetFillColor(10);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(legendTextSize);
  legend->SetTextColor(1);
  legend->SetMargin(0.20);
  legend->AddEntry(histogram1_rebinned, legendEntry1.data(), "l");
  legend->AddEntry(histogram2_rebinned, legendEntry2.data(), "l");
  if ( histogram3 ) legend->AddEntry(histogram3_rebinned, legendEntry3.data(), "p");
  if ( histogram4 ) legend->AddEntry(histogram4_rebinned, legendEntry4.data(), "p");
  legend->Draw();

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete label;
  delete legend;
  delete histogram1_rebinned;
  delete histogram2_rebinned;
  delete histogram3_rebinned;
  delete histogram4_rebinned;
  delete canvas;
}

TGraph* compGraphEfficiency(const std::string& graphName, const TH1* histogram)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int numPoints = xAxis->GetNbins();
  TGraph* graphEfficiency = new TGraph(numPoints + 2);
  graphEfficiency->SetName(graphName.data());
  graphEfficiency->SetPoint(0, xAxis->GetXmin(), 1.0);
  double integral = histogram->Integral(1, histogram->GetNbinsX());
  double sum = 0.;
  for ( int idxPoint = 1; idxPoint <= numPoints; ++idxPoint ) {
    int idxBin = idxPoint;
    double binCenter = xAxis->GetBinCenter(idxBin);
    double binContent = histogram->GetBinContent(idxBin);
    sum += binContent;
    graphEfficiency->SetPoint(idxPoint, binCenter, 1.0 - (sum/integral));
  }
  graphEfficiency->SetPoint(numPoints + 1, xAxis->GetXmax(), 0.0);
  return graphEfficiency;
} 

TGraph* compGraphROC(const std::string& graphName, const TGraph* graphEfficiency_signal, const TGraph* graphEfficiency_background, bool useLogScale)
{
  std::cout << "<compGraphROC>:" << std::endl;
  std::cout << " graphName = " << graphName << std::endl;
  assert(graphEfficiency_signal->GetN() == graphEfficiency_background->GetN());
  int numPoints = graphEfficiency_signal->GetN();
  TGraph* graphROC = new TGraph(numPoints);
  graphROC->SetName(graphName.data());
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, efficiency_signal;
    graphEfficiency_signal->GetPoint(idxPoint, x, efficiency_signal);
    double efficiency_background = graphEfficiency_background->Eval(x);
    double xROC = efficiency_signal;
    double yROC;
    if ( useLogScale ) yROC = efficiency_background;
    else yROC = 1.0 - efficiency_background;
    std::cout << "point #" << idxPoint << ": x = " << xROC << ", y = " << yROC << std::endl;
    graphROC->SetPoint(idxPoint, xROC, yROC);
  }
  return graphROC;
}

TGraph* compGraphROC(const std::string& graphName, const TH1* histogram_signal, const TH1* histogram_background, bool useLogScale)
{
  std::cout << "<compGraphROC>:" << std::endl;
  std::cout << " graphName = " << graphName << std::endl;
  assert(histogram_signal->GetNbinsX() == histogram_background->GetNbinsX());
  TGraph* graphEfficiency_signal = compGraphEfficiency(Form("%s_signal", graphName.data()), histogram_signal);
  TGraph* graphEfficiency_background = compGraphEfficiency(Form("%s_background", graphName.data()), histogram_background);
  TGraph* graphROC = compGraphROC(graphName, graphEfficiency_signal, graphEfficiency_background, useLogScale);
  delete graphEfficiency_signal;
  delete graphEfficiency_background;
  return graphROC;
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[],
		double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
		const std::string& labelText, double labelTextSize,
		double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.065);
  canvas->SetLeftMargin(0.17);
  canvas->SetBottomMargin(0.155);
  canvas->SetRightMargin(0.015); 
  canvas->SetLogx(false);
  canvas->SetLogy(useLogScale);
  canvas->Draw();
  canvas->cd();

  assert(graph1);
  graph1->SetLineColor(colors[0]);
  graph1->SetLineStyle(lineStyles[0]);
  graph1->SetLineWidth(lineWidths[0]);
  graph1->SetMarkerColor(colors[0]);
  //graph1->SetMarkerStyle(markerStyles[0]);
  //graph1->SetMarkerSize(markerSizes[0]);

  assert(graph2);
  graph2->SetLineColor(colors[1]);
  graph2->SetLineStyle(lineStyles[1]);
  graph2->SetLineWidth(lineWidths[1]);
  graph2->SetMarkerColor(colors[1]);
  graph2->SetMarkerStyle(markerStyles[1]);
  graph2->SetMarkerSize(markerSizes[1]);
  
  if ( graph3 ) {
    graph3->SetLineColor(colors[2]);
    graph3->SetLineStyle(lineStyles[2]);
    graph3->SetLineWidth(lineWidths[2]);
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetMarkerSize(markerSizes[2]);
  }

  if ( graph4 ) {
    graph4->SetLineColor(colors[3]);
    graph4->SetLineStyle(lineStyles[3]);
    graph4->SetLineWidth(lineWidths[3]);
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetMarkerSize(markerSizes[3]);
  }

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", numBinsX, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  assert(yMax > yMin);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(60);
  xAxis->SetTitleFont(43);
  //xAxis->SetLabelOffset(-0.01);
  xAxis->SetLabelSize(0.050);
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(60);
  yAxis->SetTitleFont(43);
  //yAxis->SetLabelOffset(0.010);
  yAxis->SetLabelSize(0.055);
  yAxis->SetLabelFont(42);
  yAxis->SetTickLength(0.040);  
  yAxis->SetNdivisions(505);
  
  dummyHistogram->Draw("axis");
  graph1->Draw("L");
  graph2->Draw("L");
  if ( graph3 ) graph3->Draw("p");
  if ( graph4 ) graph4->Draw("p");
  dummyHistogram->Draw("axissame");

  TPaveText* label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->AddText(labelText.data());
  label->SetTextFont(42);
  label->SetTextSize(labelTextSize);
  label->SetTextColor(1);
  label->SetTextAlign(13);
  label->Draw();

  TLegend* legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
  legend->SetFillColor(10);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(legendTextSize);
  legend->SetTextColor(1);
  legend->SetMargin(0.20);
  legend->AddEntry(graph1, legendEntry1.data(), "l");
  legend->AddEntry(graph2, legendEntry2.data(), "l");
  if ( graph3 ) legend->AddEntry(graph3, legendEntry3.data(), "p");
  if ( graph4 ) legend->AddEntry(graph4, legendEntry4.data(), "p");
  legend->Draw();

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete label;
  delete legend;
  delete canvas;
}

void makeMEMPerformancePlots_bbww_dilepton()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);

  bool makePlots_effectOfFakes = true;
  bool makePlots_effectOfSmearing = true;

  std::string inputFilePath = "/hdfs/local/veelken/hhAnalysis/2017/2019Feb06/histograms/hh_bbwwMEM_dilepton/";
  std::string inputFileName = "histograms_harvested_stage2_hh_bbwwMEM_dilepton.root";
  TString inputFileName_full = inputFilePath.data();
  if ( !inputFileName_full.EndsWith("/") ) inputFileName_full.Append("/");
  inputFileName_full.Append(inputFileName.data());
  TFile* inputFile = new TFile(inputFileName_full.Data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = " << inputFileName_full.Data() << " !!" << std::endl;
    assert(0);
  }

  std::map<bool, std::map<bool, std::string>> directories_part1; // key = apply_jetSmearing, apply_metSmearing
  directories_part1[false][false]      = "hh_bbwwMEM_dilepton_jetSmearingDisabled_metSmearingDisabled";
  directories_part1[false][true]       = "hh_bbwwMEM_dilepton_jetSmearingDisabled_metSmearingEnabled";
  directories_part1[true][false]       = "hh_bbwwMEM_dilepton_jetSmearingEnabled_metSmearingDisabled";
  directories_part1[true][true]        = "hh_bbwwMEM_dilepton_jetSmearingEnabled_metSmearingEnabled";
  
  std::map<bool, std::map<bool, std::string>> directories_part2; // key = selGenBJet_lead_isFake, selGenBJet_sublead_isFake
  directories_part2[false][false]      = "sel/mem_genuineLeadingBJet_genuineSubleadingBJet";
  directories_part2[false][true]       = "sel/mem_genuineLeadingBJet_fakeSubleadingBJet";
  directories_part2[true][false]       = "sel/mem_fakeLeadingBJet_genuineSubleadingBJet";
  directories_part2[true][true]        = "sel/mem_fakeLeadingBJet_fakeSubleadingBJet";

  std::map<bool, std::string> directories_part2_missingBJet; // key = selGenBJet_isFake
  directories_part2_missingBJet[false] = "sel/mem_missingBJet_genuineBJet";
  directories_part2_missingBJet[true]  = "sel/mem_missingBJet_fakeBJet";

  std::map<int, std::string> histogramNames; // key = { kProbSignal, kProbBackground, kLR }
  histogramNames[kProbSignal]        = "log_memProb_signal";
  histogramNames[kProbBackground]    = "log_memProb_background";
  histogramNames[kLR]                = "memLR";

  int colors[4]       = { kGreen - 6, 28, kBlue - 7, kBlack };
  int markerStyles[4] = { 20, 25, 21, 24 };
  int markerSizes[4]  = { 2, 2, 2, 2 };
  int lineStyles[4]   = { 1, 7, 1, 1 };
  int lineWidths[4]   = { 3, 3, 2, 2 };

  std::string labelText_signal = "SM HH #rightarrow b#bar{b} W^{+}W^{-} #rightarrow b#bar{b} l^{+}#nu l^{-}#bar{#nu}";
  std::string labelText_background = "SM t#bar{t} #rightarrow bW^{+} #bar{b}W^{-} #rightarrow b l^{+}#nu #bar{b} l^{-}#bar{#nu}";
  std::string labelText_signal_vs_background = Form("%s vs %s", labelText_signal.data(), labelText_background.data());

  int showHistograms_canvasSizeX = 1050;
  int showHistograms_canvasSizeY =  800;
  double showHistograms_xAxisOffset = 0.96;
  double showHistograms_yAxisOffset = 1.12;

  int showGraphs_canvasSizeX = 1050;
  int showGraphs_canvasSizeY =  800;
  double showGraphs_xAxisOffset = 0.96;
  double showGraphs_yAxisOffset = 1.12;

  int numBinsX_memLR = 60;
  double xMin_memLR = 0.;
  double xMax_memLR = 1.;
  double yMin_memLR = 1.e-5;
  double yMax_memLR = 1.e0;

  int numBinsX_probS = 75;
  double xMin_probS = -70.;
  double xMax_probS =  +5.;
  double yMin_probS = 1.e-5;
  double yMax_probS = 1.e0;

  int numBinsX_probB = 75;
  double xMin_probB = -70.;
  double xMax_probB =  +5.;
  double yMin_probB = 1.e-5;
  double yMax_probB = 1.e0;

  if ( makePlots_effectOfFakes ) {
    //-------------------------------------------------------------------------------------------------
    TH1* histogram_probS_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_fakeLeadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_fakeSubleadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_2fakeBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kSignal, histogramNames[kProbSignal]);
  
    TH1* histogram_probS_noSmearing_1fakeBJet_signal = addHistograms("histogram_probS_noSmearing_1fakeBJet_signal",
      histogram_probS_noSmearing_fakeLeadBJet_signal, 
      histogram_probS_noSmearing_fakeSubleadBJet_signal, 
      histogram_probS_noSmearing_2fakeBJets_signal);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probS_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
		   histogram_probS_noSmearing_1fakeBJet_signal, "#geq 1 fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.73, 0.33, 0.15,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset, 
		   "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_probS_signal.pdf");
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probS_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
		   histogram_probS_noSmearing_fakeLeadBJet_signal, "fake lead. b-jet",
		   histogram_probS_noSmearing_fakeSubleadBJet_signal, "fake sublead. b-jet",
		   histogram_probS_noSmearing_2fakeBJets_signal, "2 fake b-jets",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset, 
		   "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_probS_signal.pdf");

    TH1* histogram_probS_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_noSmearing_fakeBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kSignal, histogramNames[kProbSignal]);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probS_missingBJet_noSmearing_genuineBJet_signal, "genuine b-jet",
		   histogram_probS_missingBJet_noSmearing_fakeBJet_signal, "fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.73, 0.33, 0.15,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_probS_missingBJet_signal.pdf");
    
    TH1* histogram_probB_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_fakeLeadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_fakeSubleadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_2fakeBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kSignal, histogramNames[kProbBackground]);
  
    TH1* histogram_probB_noSmearing_1fakeBJet_signal = addHistograms("histogram_probB_noSmearing_1fakeBJet_signal",
      histogram_probB_noSmearing_fakeLeadBJet_signal, 
      histogram_probB_noSmearing_fakeSubleadBJet_signal, 
      histogram_probB_noSmearing_2fakeBJets_signal);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
		   histogram_probB_noSmearing_1fakeBJet_signal, "#geq 1 fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.73, 0.33, 0.15,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}", showHistograms_xAxisOffset,
		   true, yMin_probB, yMax_probB, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_probB_signal.pdf");
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
		   histogram_probB_noSmearing_fakeLeadBJet_signal, "fake lead. b-jet",
		   histogram_probB_noSmearing_fakeSubleadBJet_signal, "fake sublead. b-jet",
		   histogram_probB_noSmearing_2fakeBJets_signal, "2 fake b-jets",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}",showHistograms_xAxisOffset,
		   true, 1.e-4, 1.e0, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_probB_signal.pdf");

    TH1* histogram_probB_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_noSmearing_fakeBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kSignal, histogramNames[kProbBackground]);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_missingBJet_noSmearing_genuineBJet_signal, "genuine b-jet",
		   histogram_probB_missingBJet_noSmearing_fakeBJet_signal, "fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.73, 0.33, 0.15,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}",showHistograms_xAxisOffset,
		   true, yMin_probB, yMax_probB, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_probB_missingBJet_signal.pdf");
    //-------------------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------------------------
    TH1* histogram_ProbS_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_ProbS_noSmearing_fakeLeadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_ProbS_noSmearing_fakeSubleadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_ProbS_noSmearing_2fakeBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kBackground, histogramNames[kProbSignal]);
  
    TH1* histogram_ProbS_noSmearing_1fakeBJet_background = addHistograms("histogram_ProbS_noSmearing_1fakeBJet_background",
      histogram_ProbS_noSmearing_fakeLeadBJet_background, 
      histogram_ProbS_noSmearing_fakeSubleadBJet_background, 
      histogram_ProbS_noSmearing_2fakeBJets_background);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_ProbS_noSmearing_2genuineBJets_background, "2 genuine b-jets",
		   histogram_ProbS_noSmearing_1fakeBJet_background, "#geq 1 fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.73, 0.33, 0.15,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_probS_background.pdf");
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_ProbS_noSmearing_2genuineBJets_background, "2 genuine b-jets",
		   histogram_ProbS_noSmearing_fakeLeadBJet_background, "fake lead. b-jet",
		   histogram_ProbS_noSmearing_fakeSubleadBJet_background, "fake sublead. b-jet",
		   histogram_ProbS_noSmearing_2fakeBJets_background, "2 fake b-jets",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_probS_background.pdf");

    TH1* histogram_probS_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_noSmearing_fakeBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kBackground, histogramNames[kProbSignal]);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probS_missingBJet_noSmearing_genuineBJet_background, "genuine b-jet",
		   histogram_probS_missingBJet_noSmearing_fakeBJet_background, "fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.73, 0.33, 0.15,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_probS_missingBJet_background.pdf");

    TH1* histogram_probB_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_fakeLeadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_fakeSubleadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_2fakeBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kBackground, histogramNames[kProbBackground]);
  
    TH1* histogram_probB_noSmearing_1fakeBJet_background = addHistograms("histogram_probB_noSmearing_1fakeBJet_background",
      histogram_probB_noSmearing_fakeLeadBJet_background, 
      histogram_probB_noSmearing_fakeSubleadBJet_background, 
      histogram_probB_noSmearing_2fakeBJets_background);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_noSmearing_2genuineBJets_background, "2 genuine b-jets",
		   histogram_probB_noSmearing_1fakeBJet_background, "#geq 1 fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.73, 0.33, 0.15,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}",showHistograms_xAxisOffset,
		   true, yMin_probB, yMax_probB, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_probB_background.pdf");
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_noSmearing_2genuineBJets_background, "2 genuine b-jets",
		   histogram_probB_noSmearing_fakeLeadBJet_background, "fake lead. b-jet",
		   histogram_probB_noSmearing_fakeSubleadBJet_background, "fake sublead. b-jet",
		   histogram_probB_noSmearing_2fakeBJets_background, "2 fake b-jets",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}",showHistograms_xAxisOffset,
		   true, yMin_probB, yMax_probB, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_probB_background.pdf");

    TH1* histogram_probB_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_noSmearing_fakeBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kBackground, histogramNames[kProbBackground]);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_missingBJet_noSmearing_genuineBJet_background, "genuine b-jet",
		   histogram_probB_missingBJet_noSmearing_fakeBJet_background, "fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.73, 0.33, 0.15,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}",showHistograms_xAxisOffset,
		   true, yMin_probB, yMax_probB, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_probB_missingBJet_background.pdf");
    //-------------------------------------------------------------------------------------------------
  
    //-------------------------------------------------------------------------------------------------
    TH1* histogram_memLR_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_fakeLeadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_fakeSubleadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_2fakeBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kSignal, histogramNames[kLR]);
  
    TH1* histogram_memLR_noSmearing_1fakeBJet_signal = addHistograms("histogram_memLR_noSmearing_1fakeBJet_signal",
      histogram_memLR_noSmearing_fakeLeadBJet_signal, 
      histogram_memLR_noSmearing_fakeSubleadBJet_signal, 
      histogram_memLR_noSmearing_2fakeBJets_signal);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
		   histogram_memLR_noSmearing_1fakeBJet_signal, "#geq 1 fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.73, 0.33, 0.15,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_memLR_signal.pdf");
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
		   histogram_memLR_noSmearing_fakeLeadBJet_signal, "fake lead. b-jet",
		   histogram_memLR_noSmearing_fakeSubleadBJet_signal, "fake sublead. b-jet",
		   histogram_memLR_noSmearing_2fakeBJets_signal, "2 fake b-jets",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.62, 0.33, 0.30,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_memLR_signal.pdf");

    TH1* histogram_memLR_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_fakeLeadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_fakeSubleadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_2fakeBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kBackground, histogramNames[kLR]);
  
    TH1* histogram_memLR_noSmearing_1fakeBJet_background = addHistograms("histogram_memLR_noSmearing_1fakeBJet_background",
      histogram_memLR_noSmearing_fakeLeadBJet_background, 
      histogram_memLR_noSmearing_fakeSubleadBJet_background, 
      histogram_memLR_noSmearing_2fakeBJets_background);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_noSmearing_2genuineBJets_background, "2 genuine b-jets",
		   histogram_memLR_noSmearing_1fakeBJet_background, "#geq 1 fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.73, 0.33, 0.15,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_memLR_background.pdf");
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_noSmearing_2genuineBJets_background, "2 genuine b-jets",
		   histogram_memLR_noSmearing_fakeLeadBJet_background, "fake lead. b-jet",
		   histogram_memLR_noSmearing_fakeSubleadBJet_background, "fake sublead. b-jet",
		   histogram_memLR_noSmearing_2fakeBJets_background, "2 fake b-jets",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.62, 0.33, 0.30,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_memLR_background.pdf");

    TGraph* graph_ROC_noSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_noSmearing_2genuineBJets",
      histogram_memLR_noSmearing_2genuineBJets_signal, histogram_memLR_noSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_noSmearing_fakeLeadBJet_logScale = compGraphROC("graph_ROC_noSmearing_fakeLeadBJet",
      histogram_memLR_noSmearing_fakeLeadBJet_signal, histogram_memLR_noSmearing_fakeLeadBJet_background, true);
    TGraph* graph_ROC_noSmearing_fakeSubleadBJet_logScale = compGraphROC("graph_ROC_noSmearing_fakeSubleadBJet",
      histogram_memLR_noSmearing_fakeSubleadBJet_signal, histogram_memLR_noSmearing_fakeSubleadBJet_background, true);
    TGraph* graph_ROC_noSmearing_2fakeBJets_logScale = compGraphROC("graph_ROC_noSmearing_2fakeBJets",
      histogram_memLR_noSmearing_2fakeBJets_signal, histogram_memLR_noSmearing_2fakeBJets_background, true);

    TGraph* graph_ROC_noSmearing_1fakeBJet_logScale = compGraphROC("graph_ROC_noSmearing_",
      histogram_memLR_noSmearing_1fakeBJet_signal, histogram_memLR_noSmearing_1fakeBJet_background, true);

    showGraphs(showGraphs_canvasSizeX, showGraphs_canvasSizeY,
	       graph_ROC_noSmearing_2genuineBJets_logScale, "2 genuine b-jets",
	       graph_ROC_noSmearing_1fakeBJet_logScale, "#geq 1 fake b-jet",
	       nullptr, "",
	       nullptr, "",
	       colors, markerStyles, markerSizes, lineStyles, lineWidths,
	       0.055, 0.41, 0.73, 0.33, 0.15,
	       labelText_signal_vs_background, 0.055,
	       0.1800, 0.9525, 0.2900, 0.0900,
	       10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
	       true, 4.e-4, 1.5e0, "Background Rate", showGraphs_yAxisOffset, 
	       "hh_bbwwMEM_dilepton_effectOfFakes_2graphs_ROC.pdf");
    showGraphs(showGraphs_canvasSizeX, showGraphs_canvasSizeY,
	       graph_ROC_noSmearing_2genuineBJets_logScale, "2 genuine b-jets",
	       graph_ROC_noSmearing_fakeLeadBJet_logScale, "fake lead. b-jet",
	       graph_ROC_noSmearing_fakeSubleadBJet_logScale, "fake sublead. b-jet",
	       graph_ROC_noSmearing_2fakeBJets_logScale, "2 fake b-jets",
	       colors, markerStyles, markerSizes, lineStyles, lineWidths,
	       0.055, 0.41, 0.62, 0.33, 0.30,
	       labelText_signal_vs_background, 0.055,
	       0.1800, 0.9525, 0.2900, 0.0900,
	       10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
	       true, 4.e-4, 1.5e0, "Background Rate", showGraphs_yAxisOffset, 
	       "hh_bbwwMEM_dilepton_effectOfFakes_4graphs_ROC.pdf");

    TH1* histogram_memLR_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_noSmearing_fakeBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kSignal, histogramNames[kLR]);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_missingBJet_noSmearing_genuineBJet_signal, "genuine b-jet",
		   histogram_memLR_missingBJet_noSmearing_fakeBJet_signal, "fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.73, 0.33, 0.15,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_memLR_missingBJet_signal.pdf");

    TH1* histogram_memLR_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_noSmearing_fakeBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kBackground, histogramNames[kLR]);

    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_missingBJet_noSmearing_genuineBJet_background, "genuine b-jet",
		   histogram_memLR_missingBJet_noSmearing_fakeBJet_background, "fake b-jet",
		   nullptr, "",
		   nullptr, "",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.73, 0.33, 0.15,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfFakes_memLR_missingBJet_background.pdf");

    TGraph* graph_ROC_missingBJet_noSmearing_genuineBJet_logScale = compGraphROC("graph_ROC_missingBJet_noSmearing_genuineBJet",
      histogram_memLR_missingBJet_noSmearing_genuineBJet_signal, histogram_memLR_missingBJet_noSmearing_genuineBJet_background, true);
    TGraph* graph_ROC_missingBJet_noSmearing_fakeBJet_logScale = compGraphROC("graph_ROC_missingBJet_noSmearing_fakeBJet",
      histogram_memLR_missingBJet_noSmearing_fakeBJet_signal, histogram_memLR_missingBJet_noSmearing_fakeBJet_background, true);

    showGraphs(showGraphs_canvasSizeX, showGraphs_canvasSizeY,
	       graph_ROC_missingBJet_noSmearing_genuineBJet_logScale, "genuine b-jet",
	       graph_ROC_missingBJet_noSmearing_fakeBJet_logScale, "fake b-jet",
	       nullptr, "",
	       nullptr, "",
	       colors, markerStyles, markerSizes, lineStyles, lineWidths,
	       0.055, 0.41, 0.73, 0.33, 0.15,
	       labelText_signal_vs_background, 0.055,
	       0.1800, 0.9525, 0.2900, 0.0900,
	       10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
	       true, 4.e-4, 1.5e0, "Background Rate", showGraphs_yAxisOffset, 
	       "hh_bbwwMEM_dilepton_effectOfFakes_ROC_missingBJet.pdf");
    //-------------------------------------------------------------------------------------------------
  }

  if ( makePlots_effectOfSmearing ) {
    //-------------------------------------------------------------------------------------------------
    TH1* histogram_probS_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_jetSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_jet_and_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kSignal, histogramNames[kProbSignal]);
  
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probS_noSmearing_2genuineBJets_signal, "MC truth",
		   histogram_probS_jetSmearing_2genuineBJets_signal, "Jet smearing",
		   histogram_probS_metSmearing_2genuineBJets_signal, "MET smearing",
		   histogram_probS_jet_and_metSmearing_2genuineBJets_signal, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_probS_signal.pdf");

    TH1* histogram_probS_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_jetSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_metSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kSignal, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_jet_and_metSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kSignal, histogramNames[kProbSignal]);
    
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probS_missingBJet_noSmearing_genuineBJet_signal, "MC truth",
		   histogram_probS_missingBJet_jetSmearing_genuineBJet_signal, "Jet smearing",
		   histogram_probS_missingBJet_metSmearing_genuineBJet_signal, "MET smearing",
		   histogram_probS_missingBJet_jet_and_metSmearing_genuineBJet_signal, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_probS_missingBJet_signal.pdf");

    TH1* histogram_probB_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_jetSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_jet_and_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kSignal, histogramNames[kProbBackground]);
  
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_noSmearing_2genuineBJets_signal, "MC truth",
		   histogram_probB_jetSmearing_2genuineBJets_signal, "Jet smearing",
		   histogram_probB_metSmearing_2genuineBJets_signal, "MET smearing",
		   histogram_probB_jet_and_metSmearing_2genuineBJets_signal, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}",showHistograms_xAxisOffset,
		   true, yMin_probB, yMax_probB, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_probB_signal.pdf");

    TH1* histogram_probB_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_jetSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_metSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kSignal, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_jet_and_metSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kSignal, histogramNames[kProbBackground]);
    
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_missingBJet_noSmearing_genuineBJet_signal, "MC truth",
		   histogram_probB_missingBJet_jetSmearing_genuineBJet_signal, "Jet smearing",
		   histogram_probB_missingBJet_metSmearing_genuineBJet_signal, "MET smearing",
		   histogram_probB_missingBJet_jet_and_metSmearing_genuineBJet_signal, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}",showHistograms_xAxisOffset,
		   true, yMin_probB, yMax_probB, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_probB_missingBJet_signal.pdf");
    //-------------------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------------------------
    TH1* histogram_probS_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_probS_jetSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_probS_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_probS_jet_and_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kBackground, histogramNames[kProbSignal]);
  
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probS_noSmearing_2genuineBJets_background, "MC truth",
		   histogram_probS_jetSmearing_2genuineBJets_background, "Jet smearing",
		   histogram_probS_metSmearing_2genuineBJets_background, "MET smearing",
		   histogram_probS_jet_and_metSmearing_2genuineBJets_background, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_probS_background.pdf");

    TH1* histogram_probS_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_jetSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_metSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kBackground, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_jet_and_metSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kBackground, histogramNames[kProbSignal]);
    
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probS_missingBJet_noSmearing_genuineBJet_background, "MC truth",
		   histogram_probS_missingBJet_jetSmearing_genuineBJet_background, "Jet smearing",
		   histogram_probS_missingBJet_metSmearing_genuineBJet_background, "MET smearing",
		   histogram_probS_missingBJet_jet_and_metSmearing_genuineBJet_background, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probS, xMin_probS, xMax_probS, "w_{0}", showHistograms_xAxisOffset,
		   true, yMin_probS, yMax_probS, "dN/dw_{0}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_probS_missingBJet_background.pdf");

    TH1* histogram_probB_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_jetSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_jet_and_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kBackground, histogramNames[kProbBackground]);
  
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_noSmearing_2genuineBJets_background, "MC truth",
		   histogram_probB_jetSmearing_2genuineBJets_background, "Jet smearing",
		   histogram_probB_metSmearing_2genuineBJets_background, "MET smearing",
		   histogram_probB_jet_and_metSmearing_2genuineBJets_background, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}",showHistograms_xAxisOffset,
		   true, yMin_probB, yMax_probB, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_probB_background.pdf");

    TH1* histogram_probB_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_jetSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_metSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kBackground, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_jet_and_metSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kBackground, histogramNames[kProbBackground]);
    
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_probB_missingBJet_noSmearing_genuineBJet_background, "MC truth",
		   histogram_probB_missingBJet_jetSmearing_genuineBJet_background, "Jet smearing",
		   histogram_probB_missingBJet_metSmearing_genuineBJet_background, "MET smearing",
		   histogram_probB_missingBJet_jet_and_metSmearing_genuineBJet_background, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.31, 0.62, 0.33, 0.30,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_probB, xMin_probB, xMax_probB, "w_{1}", showHistograms_xAxisOffset,
		   true, yMin_probB, yMax_probB, "dN/dw_{1}", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_probB_missingBJet_background.pdf");
    //-------------------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------------------------
    TH1* histogram_memLR_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_jetSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_jet_and_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kSignal, histogramNames[kLR]);
  
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_noSmearing_2genuineBJets_signal, "MC truth",
		   histogram_memLR_jetSmearing_2genuineBJets_signal, "Jet smearing",
		   histogram_memLR_metSmearing_2genuineBJets_signal, "MET smearing",
		   histogram_memLR_jet_and_metSmearing_2genuineBJets_signal, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.62, 0.33, 0.30,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_memLR_signal.pdf");

    TH1* histogram_memLR_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_jetSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_jet_and_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kBackground, histogramNames[kLR]);
  
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_noSmearing_2genuineBJets_background, "MC truth",
		   histogram_memLR_jetSmearing_2genuineBJets_background, "Jet smearing",
		   histogram_memLR_metSmearing_2genuineBJets_background, "MET smearing",
		   histogram_memLR_jet_and_metSmearing_2genuineBJets_background, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.62, 0.33, 0.30,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_memLR_background.pdf");

    TGraph* graph_ROC_noSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_noSmearing_2genuineBJets",
      histogram_memLR_noSmearing_2genuineBJets_signal, histogram_memLR_noSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_jetSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_jetSmearing_2genuineBJets",
      histogram_memLR_jetSmearing_2genuineBJets_signal, histogram_memLR_jetSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_metSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_metSmearing_2genuineBJets",
      histogram_memLR_metSmearing_2genuineBJets_signal, histogram_memLR_metSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_jet_and_metSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_jet_and_metSmearing_2genuineBJets",
      histogram_memLR_jet_and_metSmearing_2genuineBJets_signal, histogram_memLR_jet_and_metSmearing_2genuineBJets_background, true);

    showGraphs(showGraphs_canvasSizeX, showGraphs_canvasSizeY,
	       graph_ROC_noSmearing_2genuineBJets_logScale, "MC truth",
	       graph_ROC_jetSmearing_2genuineBJets_logScale, "Jet smearing",
	       graph_ROC_metSmearing_2genuineBJets_logScale, "MET smearing",
	       graph_ROC_jet_and_metSmearing_2genuineBJets_logScale, "Jet+MET smearing",
	       colors, markerStyles, markerSizes, lineStyles, lineWidths,
	       0.055, 0.41, 0.62, 0.33, 0.30,
	       labelText_signal_vs_background, 0.055,
	       0.1800, 0.9525, 0.2900, 0.0900,
	       10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
	       true, 4.e-4, 1.5e0, "Background Rate", showGraphs_yAxisOffset, 
	       "hh_bbwwMEM_dilepton_effectOfSmearing_ROC.pdf");

    TH1* histogram_memLR_missingBJet_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_jetSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kSignal, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kSignal, histogramNames[kLR]);
  
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_missingBJet_noSmearing_2genuineBJets_signal, "MC truth",
		   histogram_memLR_missingBJet_jetSmearing_2genuineBJets_signal, "Jet smearing",
		   histogram_memLR_missingBJet_metSmearing_2genuineBJets_signal, "MET smearing",
		   histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_signal, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.62, 0.33, 0.30,
		   labelText_signal, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_memLR_missingBJet_signal.pdf");

    TH1* histogram_memLR_missingBJet_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_jetSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kBackground, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kBackground, histogramNames[kLR]);
  
    showHistograms(showHistograms_canvasSizeX, showHistograms_canvasSizeY,
		   histogram_memLR_missingBJet_noSmearing_2genuineBJets_background, "MC truth",
		   histogram_memLR_missingBJet_jetSmearing_2genuineBJets_background, "Jet smearing",
		   histogram_memLR_missingBJet_metSmearing_2genuineBJets_background, "MET smearing",
		   histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_background, "Jet+MET smearing",
		   colors, markerStyles, markerSizes, lineStyles, lineWidths,
		   0.055, 0.41, 0.62, 0.33, 0.30,
		   labelText_background, 0.055,
		   0.1800, 0.9525, 0.2900, 0.0900,
		   numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
		   true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
		   "hh_bbwwMEM_dilepton_effectOfSmearing_memLR_missingBJet_background.pdf");

    TGraph* graph_ROC_missingBJet_noSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_missingBJet_noSmearing_2genuineBJets",
      histogram_memLR_missingBJet_noSmearing_2genuineBJets_signal, histogram_memLR_missingBJet_noSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_missingBJet_jetSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_missingBJet_jetSmearing_2genuineBJets",
      histogram_memLR_missingBJet_jetSmearing_2genuineBJets_signal, histogram_memLR_missingBJet_jetSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_missingBJet_metSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_missingBJet_metSmearing_2genuineBJets",
      histogram_memLR_missingBJet_metSmearing_2genuineBJets_signal, histogram_memLR_missingBJet_metSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_missingBJet_jet_and_metSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_missingBJet_jet_and_metSmearing_2genuineBJets",
      histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_signal, histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_background, true);

    showGraphs(showGraphs_canvasSizeX, showGraphs_canvasSizeY,
	       graph_ROC_missingBJet_noSmearing_2genuineBJets_logScale, "MC truth",
	       graph_ROC_missingBJet_jetSmearing_2genuineBJets_logScale, "Jet smearing",
	       graph_ROC_missingBJet_metSmearing_2genuineBJets_logScale, "MET smearing",
	       graph_ROC_missingBJet_jet_and_metSmearing_2genuineBJets_logScale, "Jet+MET smearing",
	       colors, markerStyles, markerSizes, lineStyles, lineWidths,
	       0.055, 0.41, 0.62, 0.33, 0.30,
	       labelText_signal_vs_background, 0.055,
	       0.1800, 0.9525, 0.2900, 0.0900,
	       10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
	       true, 4.e-4, 1.5e0, "Background Rate", showGraphs_yAxisOffset, 
	       "hh_bbwwMEM_dilepton_effectOfSmearing_ROC_missingBJet.pdf");
    //-------------------------------------------------------------------------------------------------
  }

  delete inputFile;
}

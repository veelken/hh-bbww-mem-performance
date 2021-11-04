
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
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

enum { kUndefined, kSignal_lo, kSignal_nlo, kBackground_lo, kBackground_nlo };

enum { kMbb, kMll, kMll_vs_Mbb };

bool makePlots_png  = true;
bool makePlots_pdf  = true;
bool makePlots_root = true;

std::string getHistogramKey(int idxHistogram)
{
  if      ( idxHistogram == kMbb        ) return "mbb";
  else if ( idxHistogram == kMll        ) return "mll";
  else if ( idxHistogram == kMll_vs_Mbb ) return "mbb_vs_mll";
  else assert(0);
  return "";
}

TH1* loadHistogram(TFile* inputFile, const std::string& directory_part1, const std::string& directory_part2, int signal_or_background, const std::string& histogramName)
{  
  TString histogramName_full = directory_part1.data();
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  histogramName_full.Append(directory_part2.data());
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  if      ( signal_or_background == kSignal_lo      ) histogramName_full.Append("signal_lo");
  else if ( signal_or_background == kSignal_nlo     ) histogramName_full.Append("signal_nlo");
  else if ( signal_or_background == kBackground_lo  ) histogramName_full.Append("background_lo");
  else if ( signal_or_background == kBackground_nlo ) histogramName_full.Append("background_nlo");
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

TH1* addHistograms(const std::string& histogramSumName, const TH1* histogram1, const TH1* histogram2, const TH1* histogram3 = nullptr)
{
  TH1* histogramSum = (TH1*)histogram1->Clone(histogramSumName.data());
  histogramSum->Reset();
  if ( !histogramSum->GetSumw2N() ) histogramSum->Sumw2();
  histogramSum->Add(histogram1);
  histogramSum->Add(histogram2);
  if ( histogram3 ) histogramSum->Add(histogram3);
  double integral = 0;
  if ( dynamic_cast<TH2*>(histogramSum) )
  {
    TH2* histogramSum2d = dynamic_cast<TH2*>(histogramSum);
    integral = histogramSum2d->Integral(1, histogramSum2d->GetNbinsX(), 1, histogramSum2d->GetNbinsY()); // CV: exclude underflow and overflow bins
  }
  else
  {
    integral = histogramSum->Integral(1, histogramSum->GetNbinsX()); // CV: exclude underflow and overflow bins
  }
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
		    int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
		    double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions, 
		    const std::string& labelText, double labelTextSize,
		    double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		    int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
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
  double numBinsX_forRebinning;
  if ( xMax > xMin ) {
    const TAxis* xAxis = histogram1->GetXaxis();
    assert(xAxis->GetXmax() > xAxis->GetXmin());
    numBinsX_forRebinning = TMath::Nint(numBinsX*(xAxis->GetXmax() - xAxis->GetXmin())/(xMax - xMin));
    //std::cout << "numBinsX = " << numBinsX << ": numBinsX_forRebinning = " << numBinsX_forRebinning << std::endl;
  } else {
    numBinsX_forRebinning = numBinsX;
  }
  
  assert(histogram1);
  TH1* histogram1_rebinned = rebinHistogram(histogram1, numBinsX_forRebinning);
  histogram1_rebinned->SetFillColor(0);
  histogram1_rebinned->SetFillStyle(0);
  histogram1_rebinned->SetLineColor(colors[0]);
  histogram1_rebinned->SetLineStyle(lineStyles[0]);
  histogram1_rebinned->SetLineWidth(lineWidths[0]);
  histogram1_rebinned->SetMarkerColor(colors[0]);
  histogram1_rebinned->SetMarkerStyle(markerStyles[0]);
  histogram1_rebinned->SetMarkerSize(markerSizes[0]);

  assert(histogram2);
  TH1* histogram2_rebinned = rebinHistogram(histogram2, numBinsX_forRebinning);
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
    histogram3_rebinned = rebinHistogram(histogram3, numBinsX_forRebinning);
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
    histogram4_rebinned = rebinHistogram(histogram4, numBinsX_forRebinning);
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

  histogram1_rebinned->Draw(Form("%ssame", drawOptions[0].data()));
  histogram2_rebinned->Draw(Form("%ssame", drawOptions[1].data()));
  if ( histogram3_rebinned ) histogram3_rebinned->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogram4_rebinned ) histogram4_rebinned->Draw(Form("%ssame", drawOptions[3].data()));
  histogram1_rebinned->Draw("axissame");

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
    legend->AddEntry(histogram1_rebinned, legendEntry1.data(), legendOptions[0].data());
    legend->AddEntry(histogram2_rebinned, legendEntry2.data(), legendOptions[1].data());
    if ( histogram3 ) legend->AddEntry(histogram3_rebinned, legendEntry3.data(), legendOptions[2].data());
    if ( histogram4 ) legend->AddEntry(histogram4_rebinned, legendEntry4.data(), legendOptions[3].data());
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
  delete histogram1_rebinned;
  delete histogram2_rebinned;
  delete histogram3_rebinned;
  delete histogram4_rebinned;
  delete canvas;
}

double square(double x)
{
  return x*x;
}

TH1* compRatioHistogram(const TH1* histogram_numerator, const TH1* histogram_denominator)
{
  assert(histogram_numerator->GetNbinsX() == histogram_denominator->GetNbinsX());

  std::string histogramRatioName = Form("%s_div_%s", histogram_numerator->GetName(), histogram_denominator->GetName());
  TH1* histogramRatio = (TH1*)histogram_numerator->Clone(histogramRatioName.data());
  histogramRatio->Reset();
  if ( !histogramRatio->GetSumw2N() ) histogramRatio->Sumw2();
  histogramRatio->SetTitle("");
  histogramRatio->SetStats(false);

  int numBinsX = histogram_numerator->GetNbinsX();
  for ( int idxBin = 1; idxBin <= numBinsX; ++idxBin ) {
    double binContent_numerator   = histogram_numerator->GetBinContent(idxBin);
    double binError_numerator     = histogram_numerator->GetBinError(idxBin);
    double binContent_denominator = histogram_denominator->GetBinContent(idxBin);
    double binError_denominator   = histogram_denominator->GetBinError(idxBin);

    if ( binContent_denominator > 0. ) {
      //double binContent_ratio = binContent_numerator/binContent_denominator - 1.;
      double binContent_ratio = binContent_numerator/binContent_denominator;
      double binErr2_ratio = 0.;
      binErr2_ratio += square(binError_numerator/binContent_denominator);
      binErr2_ratio += square(binError_denominator*binContent_numerator/square(binContent_denominator));
      double binErr_ratio = TMath::Sqrt(binErr2_ratio);
      histogramRatio->SetBinContent(idxBin, binContent_ratio);
      histogramRatio->SetBinError(idxBin, binErr_ratio);
    }
  }
  
  return histogramRatio;
}

void copyHistogramStyle(const TH1* histogram_source, TH1* histogram_target)
{
  histogram_target->SetMarkerColor(histogram_source->GetMarkerColor());
  histogram_target->SetMarkerSize(histogram_source->GetMarkerSize());
  histogram_target->SetMarkerStyle(histogram_source->GetMarkerStyle());
  histogram_target->SetLineColor(histogram_source->GetLineColor());
  histogram_target->SetLineWidth(histogram_source->GetLineWidth());
  histogram_target->SetLineStyle(histogram_source->GetLineStyle());
  histogram_target->SetFillColor(histogram_source->GetFillColor());
  histogram_target->SetFillStyle(histogram_source->GetFillStyle());
}

void showHistograms_wRatio(double canvasSizeX, double canvasSizeY,
			   TH1* histogramRef, const std::string& legendEntryRef,
			   TH1* histogram2, const std::string& legendEntry2,
			   TH1* histogram3, const std::string& legendEntry3,
			   TH1* histogram4, const std::string& legendEntry4,
			   int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
			   double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions, 
			   const std::string& labelText, double labelTextSize,
			   double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
			   int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
			   bool useLogScale, double yMin, double yMax, double yMin_ratio, double yMax_ratio, const std::string& yAxisTitle, double yAxisOffset,
			   const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = 0.170;
  const double margin_bottom   = 0.155;
  const double margin_right    = 0.015;
  
  const double topPad_sizeY    = 0.69;
  const double clearance_sizeY = 0.02;
  const double bottomPad_sizeY = 0.29;
  const double sum_sizeY       = topPad_sizeY + clearance_sizeY + bottomPad_sizeY;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.);
  canvas->SetLeftMargin(0.);
  canvas->SetBottomMargin(0.);
  canvas->SetRightMargin(0.); 
  canvas->Draw();
  canvas->cd();

  const double topPad_x0 = 0.;
  const double topPad_x1 = 1.;
  const double topPad_y1 = 1.;
  const double topPad_y0 = topPad_y1 - topPad_sizeY/sum_sizeY;
  assert(topPad_x1 > topPad_x0 && topPad_y1 > topPad_y0);

  TPad* topPad = new TPad("topPad", "topPad", topPad_x0, topPad_y0, topPad_x1, topPad_y1);
  topPad->SetFillColor(10);
  topPad->SetBorderSize(0);
  topPad->SetTopMargin(margin_top);
  topPad->SetLeftMargin(margin_left);
  topPad->SetBottomMargin(0.000);
  topPad->SetRightMargin(margin_right); 
  topPad->SetLogx(false);
  topPad->SetLogy(useLogScale);
  topPad->Draw();
  topPad->cd();

  assert(histogramRef);
  double numBinsX_forRebinning;
  if ( xMax > xMin ) {
    const TAxis* xAxis = histogramRef->GetXaxis();
    assert(xAxis->GetXmax() > xAxis->GetXmin());
    numBinsX_forRebinning = TMath::Nint(numBinsX*(xAxis->GetXmax() - xAxis->GetXmin())/(xMax - xMin));
    //std::cout << "numBinsX = " << numBinsX << ": numBinsX_forRebinning = " << numBinsX_forRebinning << std::endl;
  } else {
    numBinsX_forRebinning = numBinsX;
  }
  
  assert(histogramRef);
  TH1* histogramRef_rebinned = rebinHistogram(histogramRef, numBinsX_forRebinning);
  histogramRef_rebinned->SetFillColor(0);
  histogramRef_rebinned->SetFillStyle(0);
  histogramRef_rebinned->SetLineColor(colors[0]);
  histogramRef_rebinned->SetLineStyle(lineStyles[0]);
  histogramRef_rebinned->SetLineWidth(lineWidths[0]);
  histogramRef_rebinned->SetMarkerColor(colors[0]);
  histogramRef_rebinned->SetMarkerStyle(markerStyles[0]);
  histogramRef_rebinned->SetMarkerSize(markerSizes[0]);

  assert(histogram2);
  TH1* histogram2_rebinned = rebinHistogram(histogram2, numBinsX_forRebinning);
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
    histogram3_rebinned = rebinHistogram(histogram3, numBinsX_forRebinning);
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
    histogram4_rebinned = rebinHistogram(histogram4, numBinsX_forRebinning);
    histogram4_rebinned->SetFillColor(0);
    histogram4_rebinned->SetFillStyle(0);
    histogram4_rebinned->SetLineColor(colors[3]);
    histogram4_rebinned->SetLineStyle(lineStyles[3]);
    histogram4_rebinned->SetLineWidth(lineWidths[3]);
    histogram4_rebinned->SetMarkerColor(colors[3]);
    histogram4_rebinned->SetMarkerStyle(markerStyles[3]);
    histogram4_rebinned->SetMarkerSize(markerSizes[3]);
  }
  
  TAxis* xAxis_top = histogramRef_rebinned->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitleSize(65);
  xAxis_top->SetTitleFont(43);
  //xAxis_top->SetLabelOffset(-0.01);
  xAxis_top->SetLabelSize(0.050);
  if ( xMax > xMin ) {
    xAxis_top->SetRangeUser(xMin, xMax);
  }
  xAxis_top->SetLabelFont(42);
  xAxis_top->SetTickLength(0.040);
  xAxis_top->SetNdivisions(505);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = histogramRef_rebinned->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);
  yAxis_top->SetTitleSize(65);
  yAxis_top->SetTitleFont(43);
  if ( yMax > yMin ) {
    yAxis_top->SetRangeUser(yMin, yMax);
  }
  //yAxis_top->SetLabelOffset(0.010);
  yAxis_top->SetLabelSize(0.055);
  yAxis_top->SetLabelFont(42);
  yAxis_top->SetTickLength(0.040);  
  yAxis_top->SetNdivisions(505);

  histogramRef_rebinned->SetTitle("");
  histogramRef_rebinned->SetStats(false);

  histogramRef_rebinned->Draw(drawOptions[0].data());
  histogram2_rebinned->Draw(Form("%ssame", drawOptions[1].data()));
  if ( histogram3_rebinned ) histogram3_rebinned->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogram4_rebinned ) histogram4_rebinned->Draw(Form("%ssame", drawOptions[3].data()));
  histogramRef_rebinned->Draw("axissame");

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
  if ( legendEntryRef != "" && legendEntry2 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(histogramRef_rebinned, legendEntryRef.data(), legendOptions[0].data());
    legend->AddEntry(histogram2_rebinned, legendEntry2.data(), legendOptions[1].data());
    if ( histogram3 ) legend->AddEntry(histogram3_rebinned, legendEntry3.data(), legendOptions[2].data());
    if ( histogram4 ) legend->AddEntry(histogram4_rebinned, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->cd();

  const double bottomPad_x0 = 0.;
  const double bottomPad_x1 = 1.;
  const double bottomPad_y0 = 0.;
  const double bottomPad_y1 = bottomPad_y0 + bottomPad_sizeY/sum_sizeY;
  assert(bottomPad_x1 > bottomPad_x0 && bottomPad_y1 > bottomPad_y0);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", bottomPad_x0, bottomPad_y0, bottomPad_x1, bottomPad_y1);
  bottomPad->SetFillColor(10);
  bottomPad->SetBorderSize(0);
  bottomPad->SetTopMargin(0.000);
  bottomPad->SetLeftMargin(margin_left);
  bottomPad->SetBottomMargin(margin_bottom*topPad_sizeY/bottomPad_sizeY);
  bottomPad->SetRightMargin(margin_right);
  bottomPad->SetLogx(false);
  bottomPad->SetLogy(false);
  bottomPad->Draw();
  bottomPad->cd();

  TH1* histogramRatio2_rebinned = compRatioHistogram(histogram2_rebinned, histogramRef_rebinned);
  copyHistogramStyle(histogram2_rebinned, histogramRatio2_rebinned);

  TH1* histogramRatio3_rebinned = nullptr;
  if ( histogram3_rebinned ) {
    histogramRatio3_rebinned = compRatioHistogram(histogram3_rebinned, histogramRef_rebinned);
    copyHistogramStyle(histogram3_rebinned, histogramRatio3_rebinned);
  }

  TH1* histogramRatio4_rebinned = nullptr;
  if ( histogram4_rebinned ) {
    histogramRatio4_rebinned = compRatioHistogram(histogram4_rebinned, histogramRef_rebinned);
    copyHistogramStyle(histogram4_rebinned, histogramRatio4_rebinned);
  }

  TAxis* xAxis_bottom = histogramRatio2_rebinned->GetXaxis();
  xAxis_bottom->SetTitle(xAxisTitle.data());
  xAxis_bottom->SetTitleOffset(xAxisOffset);
  xAxis_bottom->SetTitleSize(65);
  xAxis_bottom->SetTitleFont(43);
  //xAxis_bottom->SetLabelOffset(-0.01);
  xAxis_bottom->SetLabelSize(0.132);
  if ( xMax > xMin ) {
    xAxis_bottom->SetRangeUser(xMin, xMax);
  }
  xAxis_bottom->SetLabelFont(42);
  xAxis_bottom->SetTickLength(0.105);
  xAxis_bottom->SetNdivisions(505);

  TAxis* yAxis_bottom = histogramRatio2_rebinned->GetYaxis();
  yAxis_bottom->SetTitle("Ratio");
  yAxis_bottom->SetTitleOffset(yAxisOffset);
  yAxis_bottom->SetTitleSize(65);
  yAxis_bottom->SetTitleFont(43);
  //yAxis_bottom->SetLabelOffset(0.010);
  yAxis_bottom->SetLabelSize(0.132);
  yAxis_bottom->SetLabelFont(42);
  yAxis_bottom->SetTickLength(0.060);
  yAxis_bottom->SetNdivisions(505);

  histogramRatio2_rebinned->SetMinimum(yMin_ratio);
  histogramRatio2_rebinned->SetMaximum(yMax_ratio);

  histogramRatio2_rebinned->Draw(drawOptions[1].data());
  
  TGraph* graph_line = new TGraph(2);
  graph_line->SetPoint(0, xAxis_bottom->GetXmin(), 1.);
  graph_line->SetPoint(1, xAxis_bottom->GetXmax(), 1.);
  graph_line->SetLineColor(colors[0]);
  graph_line->SetLineStyle(lineStyles[0]);
  graph_line->SetLineWidth(lineWidths[0]);
  graph_line->Draw("L");

  histogramRatio2_rebinned->Draw("same");

  if ( histogramRatio3_rebinned ) histogramRatio3_rebinned->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogramRatio4_rebinned ) histogramRatio4_rebinned->Draw(Form("%ssame", drawOptions[3].data()));
  histogramRatio2_rebinned->Draw("axissame");

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
  delete histogramRef_rebinned;
  delete histogram2_rebinned;
  delete histogramRatio2_rebinned;
  delete histogram3_rebinned;
  delete histogramRatio3_rebinned;
  delete histogram4_rebinned;
  delete histogramRatio4_rebinned;
  delete graph_line;
  delete topPad;
  delete bottomPad;
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
  //std::cout << "<compGraphROC>:" << std::endl;
  //std::cout << " graphName = " << graphName << std::endl;
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
    //std::cout << "point #" << idxPoint << ": x = " << xROC << ", y = " << yROC << std::endl;
    graphROC->SetPoint(idxPoint, xROC, yROC);
  }
  return graphROC;
}

struct S_and_B
{  
  S_and_B(double S, double B)
    : S_(S)
    , B_(B)
  {}
  ~S_and_B() {}
  double S_;
  double B_;
  double S_over_B_;
};

bool isHigherS_over_B(const S_and_B& bin1, const S_and_B& bin2)
{
  return bin1.S_over_B_ > bin2.S_over_B_;
}

TGraph* compGraphROC(const std::string& graphName, const TH1* histogram_signal, const TH1* histogram_background, bool useLogScale)
{
  assert(histogram_signal && histogram_background);
  std::vector<S_and_B> bins;
  double integralS = 0.;
  double integralB = 0.;
  if ( dynamic_cast<const TH2*>(histogram_signal) && dynamic_cast<const TH2*>(histogram_background) )
  {
    const TH2* histogram2d_signal = dynamic_cast<const TH2*>(histogram_signal);
    const TH2* histogram2d_background = dynamic_cast<const TH2*>(histogram_background);
    assert(histogram2d_signal->GetNbinsX() == histogram2d_background->GetNbinsX());
    assert(histogram2d_signal->GetNbinsY() == histogram2d_background->GetNbinsY());
    int numBinsX = histogram2d_signal->GetNbinsX();
    int numBinsY = histogram2d_signal->GetNbinsY();
    for ( int idxBinX = 0; idxBinX <= (numBinsX + 1); ++idxBinX ) { // CV: include underflow and overflow bins
      for ( int idxBinY = 0; idxBinY <= (numBinsY + 1); ++idxBinY ) { // CV: include underflow and overflow bins
        double binContent_signal = histogram2d_signal->GetBinContent(idxBinX, idxBinY);      
        double binContent_background = histogram2d_background->GetBinContent(idxBinX, idxBinY);
        bins.push_back(S_and_B(binContent_signal, binContent_background));
        integralS += binContent_signal;
        integralB += binContent_background;
      }
    }
  }
  else
  {
    assert(histogram_signal->GetNbinsX() == histogram_background->GetNbinsX());
    int numBinsX = histogram_signal->GetNbinsX();
    for ( int idxBinX = 0; idxBinX <= (numBinsX + 1); ++idxBinX ) { // CV: include underflow and overflow bins
      double binContent_signal = histogram_signal->GetBinContent(idxBinX);      
      double binContent_background = histogram_background->GetBinContent(idxBinX);
      bins.push_back(S_and_B(binContent_signal, binContent_background));
      integralS += binContent_signal;
      integralB += binContent_background;
    }
  }
  std::sort(bins.begin(), bins.end(), isHigherS_over_B);
  int numPoints = bins.size();
  std::string tmpGraphName_signal = Form("%s_signal", graphName.data());
  TGraph* tmpGraph_signal = new TGraph(numPoints + 2);
  tmpGraph_signal->SetName(tmpGraphName_signal.data());
  tmpGraph_signal->SetPoint(0, 0, 1.0);
  std::string tmpGraphName_background = Form("%s_background", graphName.data());
  TGraph* tmpGraph_background = new TGraph(numPoints + 2);
  tmpGraph_background->SetName(tmpGraphName_background.data());
  tmpGraph_background->SetPoint(0, 0, 1.0);
  double sumS = 0.;
  double sumB = 0.;
  int idxPoint = 1;
  for ( std::vector<S_and_B>::const_iterator bin = bins.begin();
        bin != bins.end(); ++bin ) {
    sumS += bin->S_;
    tmpGraph_signal->SetPoint(idxPoint, idxPoint, 1.0 - (sumS/integralS));
    sumB += bin->B_;
    tmpGraph_background->SetPoint(idxPoint, idxPoint, 1.0 - (sumB/integralB));
    ++idxPoint;
  }
  tmpGraph_signal->SetPoint(numPoints + 1, numPoints + 1, 0.0);
  tmpGraph_background->SetPoint(numPoints + 1, numPoints + 1, 0.0);
  TGraph* graphROC = compGraphROC(graphName, tmpGraph_signal, tmpGraph_background, useLogScale);
  delete tmpGraph_signal;
  delete tmpGraph_background;
  return graphROC;
}

struct graphPoint
{
  graphPoint(double x, double y)
    : x_(x)
    , y_(y)
  {}
  ~graphPoint() {}
  double x_;
  double y_;
};

TGraph* sparsifyGraph(TGraph* graph, double minDeltaX = 0.025, double maxDeltaY = 0.100)
{
  std::vector<graphPoint> graphPoints_sparsified;
  double x_last = -1.e+3;
  double y_last = -1.e+3;
  int numPoints = graph->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y;
    graph->GetPoint(idxPoint, x, y);
    if ( x < 0.01 ) continue; // CV: prevent point @ zero signal efficiency and zero background rate from being drawn
    //if ( x > 0.99 ) continue; // CV: prevent point @ 100% signal efficiency and 100% background rate from being drawn
    if ( idxPoint == 0 || TMath::Abs(x - x_last) > minDeltaX || TMath::Abs(y - y_last) > maxDeltaY || idxPoint == (numPoints - 1) ) {
      graphPoints_sparsified.push_back(graphPoint(x, y));
      x_last = x;
      y_last = y;
    }
  }
  int numPoints_sparsified = graphPoints_sparsified.size();
  TGraph* graph_sparsified = new TGraph(numPoints_sparsified);
  graph_sparsified->SetName(Form("%s_sparsified", graph->GetName()));
  for ( int idxPoint = 0; idxPoint < numPoints_sparsified; ++idxPoint ) {
    const graphPoint& graphPoint_sparsified = graphPoints_sparsified[idxPoint];
    graph_sparsified->SetPoint(idxPoint, graphPoint_sparsified.x_, graphPoint_sparsified.y_);
  }
  return graph_sparsified;
}

void setLineStyle(TLine* line)
{
  line->SetLineColor(14);
  line->SetLineStyle(7);
  line->SetLineWidth(1);
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
		double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions,
		const std::string& labelText, double labelTextSize,
		double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = 0.170;
  const double margin_bottom   = 0.145;
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

  assert(graph1);
  TGraph* graph1_sparsified = sparsifyGraph(graph1);
  graph1_sparsified->SetLineColor(colors[0]);
  graph1_sparsified->SetLineStyle(lineStyles[0]);
  graph1_sparsified->SetLineWidth(lineWidths[0]);
  graph1_sparsified->SetMarkerColor(colors[0]);
  graph1_sparsified->SetMarkerStyle(markerStyles[0]);
  graph1_sparsified->SetMarkerSize(markerSizes[0]);

  TGraph* graph2_sparsified = nullptr;
  if ( graph2 ) {    
    graph2_sparsified = sparsifyGraph(graph2);
    graph2_sparsified->SetLineColor(colors[1]);
    graph2_sparsified->SetLineStyle(lineStyles[1]);
    graph2_sparsified->SetLineWidth(lineWidths[1]);
    graph2_sparsified->SetMarkerColor(colors[1]);
    graph2_sparsified->SetMarkerStyle(markerStyles[1]);
    graph2_sparsified->SetMarkerSize(markerSizes[1]);
  }
  
  TGraph* graph3_sparsified = nullptr;
  if ( graph3 ) {    
    graph3_sparsified = sparsifyGraph(graph3);
    graph3_sparsified->SetLineColor(colors[2]);
    graph3_sparsified->SetLineStyle(lineStyles[2]);
    graph3_sparsified->SetLineWidth(lineWidths[2]);
    graph3_sparsified->SetMarkerColor(colors[2]);
    graph3_sparsified->SetMarkerStyle(markerStyles[2]);
    graph3_sparsified->SetMarkerSize(markerSizes[2]);
  }
  
  TGraph* graph4_sparsified = nullptr;
  if ( graph4 ) {
    graph4_sparsified = sparsifyGraph(graph4);
    graph4_sparsified->SetLineColor(colors[3]);
    graph4_sparsified->SetLineStyle(lineStyles[3]);
    graph4_sparsified->SetLineWidth(lineWidths[3]);
    graph4_sparsified->SetMarkerColor(colors[3]);
    graph4_sparsified->SetMarkerStyle(markerStyles[3]);
    graph4_sparsified->SetMarkerSize(markerSizes[3]);
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
  xAxis->SetTitleSize(55);
  xAxis->SetTitleFont(43);
  //xAxis->SetLabelOffset(-0.01);
  xAxis->SetLabelSize(0.050);
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(55);
  yAxis->SetTitleFont(43);
  //yAxis->SetLabelOffset(0.010);
  yAxis->SetLabelSize(0.050);
  yAxis->SetLabelFont(42);
  yAxis->SetTickLength(0.040);  
  yAxis->SetNdivisions(505);
  
  dummyHistogram->Draw("axis");
  graph1_sparsified->Draw(drawOptions[0].data());
  if ( graph2_sparsified ) graph2_sparsified->Draw(drawOptions[1].data());
  if ( graph3_sparsified ) graph3_sparsified->Draw(drawOptions[2].data());
  if ( graph4_sparsified ) graph4_sparsified->Draw(drawOptions[3].data());

  std::vector<TLine*> lines;
  if ( outputFileName == "hh_bbwwMEM_dilepton_effectOfFakes_2graphs_ROC.pdf" ) {
    double x = 0.35;
    assert(graph1);
    double y1 = graph1->Eval(x);
    TLine* line1_horizontal = new TLine(0., y1, x, y1);
    setLineStyle(line1_horizontal);
    line1_horizontal->Draw();
    lines.push_back(line1_horizontal);
    assert(graph2);
    double y2 = graph2->Eval(x);
    TLine* line2_horizontal = new TLine(0., y2, x, y2);
    setLineStyle(line2_horizontal);
    line2_horizontal->Draw();
    lines.push_back(line2_horizontal);
    TLine* line_vertical = new TLine(x, yMin, x, TMath::Max(y1, y2));
    setLineStyle(line_vertical);
    line_vertical->Draw();
    lines.push_back(line_vertical);
  } else if ( outputFileName == "hh_bbwwMEM_dilepton_effectOfFakes_ROC_missingBJet.pdf" ) {
    double x = 0.35;
    assert(graph1);
    double y = graph1->Eval(x);
    TLine* line_horizontal = new TLine(0., y, x, y);
    setLineStyle(line_horizontal);
    line_horizontal->Draw();
    lines.push_back(line_horizontal);
    TLine* line_vertical = new TLine(x, yMin, x, y);
    setLineStyle(line_vertical);
    line_vertical->Draw();
    lines.push_back(line_vertical);
  }

  dummyHistogram->Draw("axissame");

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
  if ( graph2_sparsified || graph3_sparsified || graph4_sparsified ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(graph1_sparsified, legendEntry1.data(), legendOptions[0].data());
    if ( graph2_sparsified ) legend->AddEntry(graph2_sparsified, legendEntry2.data(), legendOptions[1].data());
    if ( graph3_sparsified ) legend->AddEntry(graph3_sparsified, legendEntry3.data(), legendOptions[2].data());
    if ( graph4_sparsified ) legend->AddEntry(graph4_sparsified, legendEntry4.data(), legendOptions[3].data());
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
  delete graph1_sparsified;
  delete graph2_sparsified;
  delete graph3_sparsified;
  delete graph4_sparsified;
  for ( std::vector<TLine*>::iterator line = lines.begin();
        line != lines.end(); ++line ) {
    delete (*line);
  }
  delete canvas;
}

void makeControlPlots_bbww_dilepton()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);

  bool makePlots_signal_vs_background = true;
  bool makePlots_effectOfFakes = true;
  bool makePlots_effectOfSmearing = true;
  bool makePlots_effectOfHigherOrders = true;

  std::string inputFilePath = "/hdfs/local/veelken/hhAnalysis/2016/2021Aug31v2/histograms/hh_bbwwMEM_dilepton/";
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
  
  std::map<int, std::string> directories_part2; // key = number of genuine b-jets
  directories_part2[2]                 = "sel/evt_2genuineBJets_";
  directories_part2[1]                 = "sel/evt_1genuineBJet";
  directories_part2[0]                 = "sel/evt_0genuineBJets_";
  
  std::map<int, std::string> histogramNames; // key = { kMbb, kMll, kMll_vs_Mbb }
  histogramNames[kMbb]                 = "mbb";
  histogramNames[kMll]                 = "mll";
  histogramNames[kMll_vs_Mbb]          = "mll_vs_mbb";

  std::string labelText_signal = "HH #rightarrow b#bar{b} WW^{*} #rightarrow b#bar{b} l^{+}#nu l^{-}#bar{#nu}";
  std::string labelText_background = "t#bar{t} #rightarrow bW #bar{b}W #rightarrow b l^{+}#nu #bar{b} l^{-}#bar{#nu}";
  //std::string labelText_signal_vs_background = Form("%s vs %s", labelText_signal.data(), labelText_background.data());
  std::string labelText_signal_vs_background = "";

  int showHistograms_canvasSizeX = 1050;
  int showHistograms_canvasSizeY =  950;
  int showHistograms_canvasSizeY_wRatio = 1150;
  double showHistograms_xAxisOffset = 0.96;
  double showHistograms_xAxisOffset_wRatio = 2.90;
  double showHistograms_yAxisOffset = 1.21;
  double showHistograms_yAxisOffset_wRatio = 1.45;
  int showHistograms_colors[4]       = { kGreen - 6, kBlack, kBlue - 7, 28 };
  int showHistograms_markerStyles[4] = { 20, 24, 21, 25 };
  int showHistograms_markerSizes[4]  = { 2, 2, 2, 2 };
  int showHistograms_lineStyles[4]   = { 1, 1, 1, 7 };
  int showHistograms_lineWidths[4]   = { 3, 2, 2, 3 };
  std::vector<std::string> showHistograms_drawOptions = { "hist", "ep", "ep", "hist" };
  std::vector<std::string> showHistograms_legendOptions = { "l", "p", "p", "l" };

  int showHistograms_signal_vs_background_colors[2]       = { kBlack, kRed };
  int showHistograms_signal_vs_background_markerStyles[2] = { 20, 24 };
  int showHistograms_signal_vs_background_markerSizes[2]  = { 2, 2 };
  int showHistograms_signal_vs_background_lineStyles[2]   = { 1, 1 };
  int showHistograms_signal_vs_background_lineWidths[2]   = { 2, 2 };
  std::vector<std::string> showHistograms_signal_vs_background_drawOptions = { "ep", "ep" };
  std::vector<std::string> showHistograms_signal_vs_background_legendOptions = { "p", "p" };

  std::map<std::string, std::string> xAxisTitle;               // key = getHistogramKey(idxHistogram)
  std::map<std::string, int>         numBinsX;                 // key = getHistogramKey(idxHistogram)
  std::map<std::string, double>      xMin;                     // key = getHistogramKey(idxHistogram)
  std::map<std::string, double>      xMax;
  std::map<std::string, std::string> yAxisTitle;               // key = getHistogramKey(idxHistogram)
  std::map<std::string, double>      yMin;                     // key = getHistogramKey(idxHistogram)
  std::map<std::string, double>      yMin_wRatio;              // key = getHistogramKey(idxHistogram)
  std::map<std::string, double>      yMax;                     // key = getHistogramKey(idxHistogram)

  xAxisTitle["mbb"]                 = "m_{bb} [GeV]";
  numBinsX["mbb"]                   =  40;
  xMin["mbb"]                       =   0.;
  xMax["mbb"]                       = 200.;
  yAxisTitle["mbb"]                 = "dN/dm_{bb} [1/GeV]";
  yMin["mbb"]                       = 1.1e-5;
  yMin_wRatio["mbb"]                = 1.1e-5;
  yMax["mbb"]                       = 1.9e0;

  xAxisTitle["mll"]                 = "m_{ll} [GeV]";
  numBinsX["mll"]                   =  40;
  xMin["mll"]                       =   0.;
  xMax["mll"]                       = 200.;
  yAxisTitle["mll"]                 = "dN/dm_{ll} [1/GeV]";
  yMin["mll"]                       = 1.1e-5;
  yMin_wRatio["mll"]                = 1.1e-5;
  yMax["mll"]                       = 1.9e0;
  
  int showGraphs_canvasSizeX = 1050;
  int showGraphs_canvasSizeY =  950;
  double showGraphs_xAxisOffset = 1.18;
  double showGraphs_yAxisOffset = 1.44;
  int showGraphs_colors[4]       = { kGreen - 6, kBlack, kBlue - 7, 28 };
  int showGraphs_markerStyles[4] = { 20, 24, 21, 25 };
  int showGraphs_markerSizes[4]  = { 2, 2, 2, 2 };
  int showGraphs_lineStyles[4]   = { 1, 1, 1, 1 };
  int showGraphs_lineWidths[4]   = { 2, 2, 2, 2 };
  std::vector<std::string> showGraphs_drawOptions = { "Lp", "Lp", "Lp", "Lp" };
  std::vector<std::string> showGraphs_legendOptions = { "lp", "lp", "lp", "lp" };

  if ( makePlots_signal_vs_background ) {
    for ( int idxHistogram = kMbb; idxHistogram <= kMll; ++idxHistogram ) {
      const std::string& histogramName = histogramNames[idxHistogram];
      std::string histogramKey = getHistogramKey(idxHistogram);

      TH1* histogram_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kSignal_lo, histogramName);
      TH1* histogram_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kBackground_lo, histogramName);

      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_noSmearing_2genuineBJets_signal, "Signal",
        histogram_noSmearing_2genuineBJets_background, "Background",
        nullptr, "",
        nullptr, "",
        showHistograms_signal_vs_background_colors, showHistograms_signal_vs_background_markerStyles, showHistograms_signal_vs_background_markerSizes, 
        showHistograms_signal_vs_background_lineStyles, showHistograms_signal_vs_background_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_signal_vs_background_legendOptions,
        "", 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        numBinsX[histogramKey], xMin[histogramKey], xMax[histogramKey], xAxisTitle[histogramKey], showHistograms_xAxisOffset,
        true, yMin[histogramKey], yMax[histogramKey], yAxisTitle[histogramKey], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_signal_vs_background_%s_unsmeared.pdf", histogramKey.data()));

      if ( idxHistogram == kMbb || idxHistogram == kMll || idxHistogram == kMll_vs_Mbb ) {
        TGraph* graph_ROC_noSmearing_2genuineBJets_logScale = compGraphROC(
          "graph_ROC_noSmearing_2genuineBJets",
          histogram_noSmearing_2genuineBJets_signal, 
          histogram_noSmearing_2genuineBJets_background, true);

        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          graph_ROC_noSmearing_2genuineBJets_logScale, "",
          nullptr, "",
          nullptr, "",
          nullptr, "",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.055, 0.23, 0.86, 0.33, 0.08, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 1.e-3, 1.e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_ROC_unsmeared.pdf");
      }
    }
  }

  if ( makePlots_effectOfFakes ) {
    for ( int idxHistogram = kMbb; idxHistogram <= kMbb; ++idxHistogram ) {
      const std::string& histogramName = histogramNames[idxHistogram];
      std::string histogramKey = getHistogramKey(idxHistogram);

      TH1* histogram_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kSignal_lo, histogramName);
      TH1* histogram_noSmearing_1genuineBJet_signal = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[1], kSignal_lo, histogramName);
      TH1* histogram_noSmearing_0genuineBJets_signal = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[0], kSignal_lo, histogramName);
  
      TH1* histogram_noSmearing_geq1fakeBJet_signal = addHistograms(
        Form("histogram_%s_noSmearing_geq1fakeBJet_signal", histogramKey.data()),
        histogram_noSmearing_1genuineBJet_signal, 
        histogram_noSmearing_0genuineBJets_signal);

      TH1* histogram_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kBackground_lo, histogramName);
      TH1* histogram_noSmearing_1genuineBJet_background = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[1], kBackground_lo, histogramName);
      TH1* histogram_noSmearing_0genuineBJets_background = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[0], kBackground_lo, histogramName);
  
      TH1* histogram_noSmearing_geq1fakeBJet_background = addHistograms(
        Form("histogram_%s_noSmearing_geq1fakeBJet_background", histogramKey.data()),
        histogram_noSmearing_1genuineBJet_background, 
        histogram_noSmearing_0genuineBJets_background);

      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
        histogram_noSmearing_geq1fakeBJet_signal, "#geq 1 fake b-jet",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        numBinsX[histogramKey], xMin[histogramKey], xMax[histogramKey], xAxisTitle[histogramKey], showHistograms_xAxisOffset,
        true, yMin[histogramKey], yMax[histogramKey], yAxisTitle[histogramKey], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_effectOfFakes_2histograms_%s_signal.pdf", histogramKey.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
        histogram_noSmearing_1genuineBJet_signal, "1 genuine + 1 fake b-jet",
        histogram_noSmearing_0genuineBJets_signal, "2 fake b-jets",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.045, 0.23, 0.66, 0.33, 0.28, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        numBinsX[histogramKey], xMin[histogramKey], xMax[histogramKey], xAxisTitle[histogramKey], showHistograms_xAxisOffset,
        true, yMin[histogramKey], yMax[histogramKey], yAxisTitle[histogramKey], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_effectOfFakes_3histograms_%s_signal.pdf", histogramKey.data()));showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_noSmearing_2genuineBJets_background, "2 genuine b-jets",
        histogram_noSmearing_geq1fakeBJet_background, "#geq 1 fake b-jet",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        numBinsX[histogramKey], xMin[histogramKey], xMax[histogramKey], xAxisTitle[histogramKey], showHistograms_xAxisOffset,
        true, yMin[histogramKey], yMax[histogramKey], yAxisTitle[histogramKey], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_effectOfFakes_2histograms_%s_background.pdf", histogramKey.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_noSmearing_2genuineBJets_background, "2 genuine b-jets",
        histogram_noSmearing_1genuineBJet_background, "1 genuine + 1 fake b-jet",
        histogram_noSmearing_0genuineBJets_background, "2 fake b-jets",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.045, 0.23, 0.66, 0.33, 0.28, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        numBinsX[histogramKey], xMin[histogramKey], xMax[histogramKey], xAxisTitle[histogramKey], showHistograms_xAxisOffset,
        true, yMin[histogramKey], yMax[histogramKey], yAxisTitle[histogramKey], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_effectOfFakes_3histograms_%s_background.pdf", histogramKey.data()));
    }
  }

  if ( makePlots_effectOfSmearing ) {
    for ( int idxHistogram = kMbb; idxHistogram <= kMbb; ++idxHistogram ) {
      const std::string& histogramName = histogramNames[idxHistogram];
      std::string histogramKey = getHistogramKey(idxHistogram);

      TH1* histogram_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kSignal_lo, histogramName);
      TH1* histogram_jetSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
        directories_part1[true][false], directories_part2[2], kSignal_lo, histogramName);

      TH1* histogram_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kBackground_lo, histogramName);
      TH1* histogram_jetSmearing_2genuineBJets_background = loadHistogram(inputFile, 
        directories_part1[true][false], directories_part2[2], kBackground_lo, histogramName);
  
      showHistograms_wRatio(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
        histogram_noSmearing_2genuineBJets_signal, "MC truth",
        histogram_jetSmearing_2genuineBJets_signal, "E_{b} smearing",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.065, 0.23, 0.75, 0.52, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        numBinsX[histogramKey], xMin[histogramKey], xMax[histogramKey], xAxisTitle[histogramKey], showHistograms_xAxisOffset_wRatio,
        true, yMin_wRatio[histogramKey], yMax[histogramKey], 1. - 0.59, 1. + 0.59, yAxisTitle[histogramKey], showHistograms_yAxisOffset_wRatio,
        Form("hh_bbwwMEM_dilepton_effectOfSmearing_%s_signal.pdf", histogramKey.data()));
      showHistograms_wRatio(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
        histogram_noSmearing_2genuineBJets_background, "MC truth",
        histogram_jetSmearing_2genuineBJets_background, "E_{b} smearing",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.065, 0.23, 0.75, 0.52, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        numBinsX[histogramKey], xMin[histogramKey], xMax[histogramKey], xAxisTitle[histogramKey], showHistograms_xAxisOffset_wRatio,
        true, yMin_wRatio[histogramKey], yMax[histogramKey], 1. - 0.59, 1. + 0.59, yAxisTitle[histogramKey], showHistograms_yAxisOffset_wRatio,
        Form("hh_bbwwMEM_dilepton_effectOfSmearing_%s_background.pdf", histogramKey.data()));
    }
  }

  if ( makePlots_effectOfHigherOrders ) {
    for ( int idxHistogram = kMbb; idxHistogram <= kMll; ++idxHistogram ) {
      const std::string& histogramName = histogramNames[idxHistogram];
      std::string histogramKey = getHistogramKey(idxHistogram);

      TH1* histogram_lo_signal = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kSignal_lo, histogramName);
      TH1* histogram_nlo_signal = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kSignal_nlo, histogramName);

      TH1* histogram_lo_background = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kBackground_lo, histogramName);
      TH1* histogram_nlo_background = loadHistogram(inputFile, 
        directories_part1[false][false], directories_part2[2], kBackground_nlo, histogramName);

      showHistograms_wRatio(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
        histogram_lo_signal, "LO",
        histogram_nlo_signal, "NLO",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.060, 0.26, 0.72, 0.27, 0.18, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        numBinsX[histogramKey], xMin[histogramKey], xMax[histogramKey], xAxisTitle[histogramKey], showHistograms_xAxisOffset_wRatio,
        true, yMin_wRatio[histogramKey], yMax[histogramKey], 1. - 0.29, 1. + 0.29, yAxisTitle[histogramKey], showHistograms_yAxisOffset_wRatio,
        Form("hh_bbwwMEM_dilepton_lo_vs_nlo_%s_signal.pdf", histogramKey.data()));
      showHistograms_wRatio(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
        histogram_lo_background, "LO",
        histogram_nlo_background, "NLO",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.060, 0.26, 0.72, 0.27, 0.18, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        numBinsX[histogramKey], xMin[histogramKey], xMax[histogramKey], xAxisTitle[histogramKey], showHistograms_xAxisOffset_wRatio,
        true, yMin_wRatio[histogramKey], yMax[histogramKey], 1. - 0.29, 1. + 0.29, yAxisTitle[histogramKey], showHistograms_yAxisOffset_wRatio,
        Form("hh_bbwwMEM_dilepton_lo_vs_nlo_%s_background.pdf", histogramKey.data()));
    }
  }

  delete inputFile;
}

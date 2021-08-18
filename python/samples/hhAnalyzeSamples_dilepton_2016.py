from hhAnalysis.bbww.samples.hhAnalyzeSamples_2016 import samples_2016
#
# Define MC samples to produce MEM performance plots for dilepton channel of the HH->bbWW analysis
#
for sample_name, sample_info in samples_2016.items():
  if sample_name == 'sum_events':
    continue
  if sample_info["process_name_specific"] in [
        "signal_ggf_nonresonant_node_sm_hh_2b2v", # HH signal sample @ LO 
      ]:
    sample_info["sample_category"] = "signal_lo"
    sample_info["use_it"] = True
  elif sample_info["process_name_specific"] in [
        "signal_ggf_nonresonant_cHHH1_hh_2b2v",   # HH signal sample @ NLO
      ]:
    sample_info["sample_category"] = "signal_nlo"
    sample_info["use_it"] = True
  elif sample_info["process_name_specific"] in [
        "TTJets_DiLept",                          # ttbar background sample @ LO
        "TTJets_DiLept_ext1",                     # ttbar background sample @ LO ("extension" sample with additional event statistics)
      ]:
    sample_info["sample_category"] = "background_lo"
    sample_info["use_it"] = True
  elif sample_info["process_name_specific"] in [
        "TTTo2L2Nu"                               # ttbar background sample @ NLO
      ]:
    sample_info["sample_category"] = "background_nlo"
    sample_info["use_it"] = True
  else:
    sample_info["use_it"] = False

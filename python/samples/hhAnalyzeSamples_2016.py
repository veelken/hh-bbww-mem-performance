from hhAnalysis.bbww.samples.hhAnalyzeSamples_2016 import samples_2016

for sample_name, sample_info in samples_2016.items():
  if sample_name == 'sum_events':
    continue
  if sample_info["process_name_specific"] in [
        "signal_ggf_nonresonant_node_sm_hh_2b2v", # HH signal sample @ LO 
        "signal_ggf_nonresonant_cHHH1_hh_2b2v",   # HH signal sample @ NLO
        "TTJets_DiLept",                          # ttbar background sample @ LO
        "TTJets_DiLept_ext1",                     # ttbar background sample @ LO ("extension" sample with additional event statistics)
        "TTTo2L2Nu"                               # ttbar background sample @ NLO
      ]:
    sample_info["use_it"] = True
  else:
    sample_info["use_it"] = False

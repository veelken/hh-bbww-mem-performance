from hhAnalysis.bbww.samples.hhAnalyzeSamples_2017 import samples_2017

for sample_name, sample_info in samples_2017.items():
  if sample_name == 'sum_events':
    continue
  if sample_info["process_name_specific"] in [
        "signal_ggf_nonresonant_node_sm_hh_2b2v", # HH signal sample @ LO 
        "signal_ggf_nonresonant_cHHH1_hh_2b2v",   # HH signal sample @ NLO
        "TTJets_DiLept",                          # ttbar background sample @ LO
        "TTTo2L2Nu",                              # ttbar background sample @ NLO
        "TTTo2L2Nu_PSweights"                     # ttbar background sample @ NLO ("extension" sample with additional event statistics) 
      ]:
    sample_info["use_it"] = True
  else:
    sample_info["use_it"] = False

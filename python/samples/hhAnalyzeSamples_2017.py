from hhAnalysis.bbww.samples.hhAnalyzeSamples_2017 import samples_2017

for sample_name, sample_info in samples_2017.items():
  if sample_name == 'sum_events':
    continue
  if sample_info["process_name_specific"] in [
        "signal_ggf_nonresonant_node_sm_hh_2b2v",
        "TTTo2L2Nu",
        "TTTo2L2Nu_PSweights"
      ]:
    sample_info["use_it"] = True
  else:
    sample_info["use_it"] = False

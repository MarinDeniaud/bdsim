execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# data stored using HepMC2 library but loaded using HepMC3 in BDSIM
if (USE_HEPMC3)
  simple_testing(bunch-eventgeneratorfile-hepmc2                 "--file=eventgeneratorfile-hepmc2.gmad"                    "")
  simple_testing(bunch-eventgeneratorfile-hepmc2-ngenerate       "--file=eventgeneratorfile-hepmc2-ngenerate.gmad --ngenerate=2"  "")	
  simple_testing(bunch-eventgeneratorfile-hepmc2-filtered        "--file=eventgeneratorfile-hepmc2-filtered.gmad"           "")
  simple_testing(bunch-eventgeneratorfile-hepmc2-filtered2       "--file=eventgeneratorfile-hepmc2-filtered2.gmad"          "")
  simple_testing(bunch-eventgeneratorfile-hepmc2-loop-ntimes     "--file=eventgeneratorfile-hepmc2-loop-ntimes.gmad"        "")
  simple_testing(bunch-eventgeneratorfile-hepmc2-lowphysics      "--file=eventgeneratorfile-hepmc2-lowphysics.gmad"         "")
  simple_testing(bunch-eventgeneratorfile-hepmc2-offset          "--file=eventgeneratorfile-hepmc2-offset.gmad"             "")
  simple_testing(bunch-eventgeneratorfile-hepmc2-recreate        "--file=eventgeneratorfile-hepmc2-recreate.gmad"           "")
  simple_testing(bunch-eventgeneratorfile-hepmc2-permit-unstable "--file=eventgeneratorfile-hepmc2-permit-unstable.gmad"    "")
  simple_testing_w_string(bunch-eventgeneratorfile-hepmc2-no-part-event   "--file=eventgeneratorfile-hepmc2-no-particles-event.gmad" "WARNING")
  simple_testing_w_string(bunch-eventgeneratorfile-hepmc2-no-part-file    "--file=eventgeneratorfile-hepmc2-no-particles-file.gmad"  "WARNING")
  simple_testing_w_string(bunch-eventgeneratorfile-hepmc2-loop            "--file=eventgeneratorfile-hepmc2-loop.gmad"               "BDSPrimaryGeneratorFileHEPMC::OpenFile")
  simple_testing_w_string(bunch-eventgeneratorfile-hepmc2-nevents-skip    "--file=eventgeneratorfile-hepmc2-nEventsSkip.gmad"        "BDSPrimaryGeneratorFileHEPMC::SkipEvents")
  simple_testing_w_string(bunch-eventgeneratorfile-hepmc2-end-run-early   "--file=eventgeneratorfile-hepmc2-end-run-early.gmad"      "2 events were skipped")
  simple_testing_w_string(bunch-eventgeneratorfile-hepmc2-loop-bad        "--file=eventgeneratorfile-hepmc2-loop-bad.gmad"           "ending")
endif()

int loadXephyr(){

  gROOT->ProcessLine(".L ../Xephyr/src/XeVersion.cxx+g");
  gROOT->ProcessLine(".L ../Xephyr/src/XeUtils.cxx+g");
  gROOT->ProcessLine(".L ../Xephyr/src/XeStat.cxx+g");
  gROOT->ProcessLine(".L ../Xephyr/src/dataHandler.cxx+g");
  gROOT->ProcessLine(".L ../Xephyr/src/XePdfObjects.cxx+g");
  gROOT->ProcessLine(".L ../Xephyr/src/XeLikelihoods.cxx+g");
  gROOT->ProcessLine(".L ../Xephyr/src/AsymptoticExclusion.cxx+g");
  gROOT->ProcessLine(".L ../Xephyr/src/ToyGenerator.cxx+g");
  gROOT->ProcessLine(".L ../Xephyr/src/plotHelpers.cxx+g");
  gROOT->ProcessLine(".L ../Xephyr/src/ToyFitterExclusion.cxx+g");
    gInterpreter->AddIncludePath("../Xephyr/src"); // in this case is just XEPHYR src from next dir.
    gROOT->ProcessLine(".L ../SR1/StatisticalAnalyses/xephyr_sr1_likelihood/src/likelihoodDef.cxx");
  return 0;
}

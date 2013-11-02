#!/bin/csh -f

rm testNSVfitTrackLikelihoods4_ZplusJets_madgraph_tautau_cfg.py
rm testNSVfitTrackLikelihoods4_Higgs_125_tautau_cfg.py
rm testNSVfitTrackLikelihoods4_ZplusJets_madgraph_mutau_cfg.py
rm testNSVfitTrackLikelihoods4_Higgs_125_mutau_cfg.py

sed 's/#__//g;s/#sample_type#/Z/g;s/#channel#/tautau/g;s/#massPoint#/None/g' testNSVfitTrackLikelihoods4_cfg.py > testNSVfitTrackLikelihoods4_ZplusJets_madgraph_tautau_cfg.py
sed 's/#__//g;s/#sample_type#/Higgs/g;s/#channel#/tautau/g;s/#massPoint#/125/g' testNSVfitTrackLikelihoods4_cfg.py > testNSVfitTrackLikelihoods4_Higgs_125_tautau_cfg.py
sed 's/#__//g;s/#sample_type#/Z/g;s/#channel#/mutau/g;s/#massPoint#/None/g' testNSVfitTrackLikelihoods4_cfg.py > testNSVfitTrackLikelihoods4_ZplusJets_madgraph_mutau_cfg.py
sed 's/#__//g;s/#sample_type#/Higgs/g;s/#channel#/mutau/g;s/#massPoint#/125/g' testNSVfitTrackLikelihoods4_cfg.py > testNSVfitTrackLikelihoods4_Higgs_125_mutau_cfg.py

rm testNSVfitTrackLikelihoods4_ZplusJets_madgraph_tautau_2013Aug08.out
rm testNSVfitTrackLikelihoods4_Higgs_125_tautau_2013Aug08.out
rm testNSVfitTrackLikelihoods4_ZplusJets_madgraph_mutau_2013Aug08.out
rm testNSVfitTrackLikelihoods4_Higgs_125_mutau_2013Aug08.out

nice cmsRun testNSVfitTrackLikelihoods4_ZplusJets_madgraph_tautau_cfg.py >&! testNSVfitTrackLikelihoods4_ZplusJets_madgraph_tautau_2013Aug08.out &
sleep 1800
nice cmsRun testNSVfitTrackLikelihoods4_Higgs_125_tautau_cfg.py >&! testNSVfitTrackLikelihoods4_Higgs_125_tautau_2013Aug08.out &
sleep 1800
nice cmsRun testNSVfitTrackLikelihoods4_ZplusJets_madgraph_mutau_cfg.py >&! testNSVfitTrackLikelihoods4_ZplusJets_madgraph_mutau_2013Aug08.out &
sleep 1800
nice cmsRun testNSVfitTrackLikelihoods4_Higgs_125_mutau_cfg.py >&! testNSVfitTrackLikelihoods4_Higgs_125_mutau_2013Aug08.out &


/*
Created:         4 September 2016
Last Updated:    4 September 2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109

-----

Make TEfficiencies for measuring efficiencies
https://root.cern.ch/doc/master/classTEfficiency.html

*/
#include "Analysis/CyMiniAna/interface/efficiency.h"


efficiency::efficiency(configuration &cmaConfig) : 
  efficiencyBase::efficiencyBase(cmaConfig){
    m_isTtbar = m_config->isTtbar();
    m_isOneLeptonAnalysis = m_config->isOneLeptonAnalysis();
}

efficiency::~efficiency() {}


void efficiency::initialize( TFile& outputFile ){
    /* Book efficiencies -- modify/inherit this function for analysis-specific efficiencies */
    outputFile.cd();

    efficiencyBase::init_eff("dilution_incl",    1, 0, 1);     // delta|y| reconstruction efficiency (inclusive)
    efficiencyBase::init_eff("dilution_mtt",    50, 0, 5000);  // delta|y| reconstruction efficiency (mttbar)
    efficiencyBase::init_eff("dilution_ytt",    10, 0, 5);     // delta|y| reconstruction efficiency (yttbar)
    efficiencyBase::init_eff("dilution_pttt",   10, 0, 1000);  // delta|y| reconstruction efficiency (pT-ttbar)
    efficiencyBase::init_eff("dilution_betatt", 10, 0, 1);     // delta|y| reconstruction efficiency (beta-ttbar)

    return;
}

void efficiency::fill( Event &event, const std::vector<unsigned int>& evtsel_decisions ){
    /* Fill efficiencies -- just use information from the event 
       This is the function to modify / inherit for analysis-specific purposes
       Example
       Fill an efficiency for jet trigger vs leading jet pT 
    */
    Ttbar1L ttbarSL = event.ttbar1L();

    double dy(0.0);
    TLorentzVector top_lep;
    TLorentzVector top_had;
    TLorentzVector ttbar;

    if (m_isOneLeptonAnalysis){
        dy      = ttbarSL.dy;
        top_lep = ttbarSL.leptop;
        top_had = ttbarSL.ljet.p4;
        ttbar   = top_had+top_lep;
    }

    float mtt    = ttbar.M();
    float pttt   = ttbar.Pt();
    float ytt    = std::abs(ttbar.Rapidity());
    float betatt = std::abs(top_had.Pz() + top_lep.Pz()) / (top_had.E() + top_lep.E());

    float true_dy(0.0);
    if (m_isTtbar){
        std::vector<TruthTop> truth_tops = event.truth();
        std::vector<Parton> truth_partons = event.truth_partons();
        if (truth_tops.size()==2){
            TruthTop top0 = truth_tops.at(0);
            TruthTop top1 = truth_tops.at(1);
            Parton ptop0  = truth_partons.at( top0.Top );
            Parton ptop1  = truth_partons.at( top1.Top );

            true_dy = (top0.isTop && top1.isAntiTop) ? std::abs(ptop0.p4.Rapidity()) - std::abs(ptop1.p4.Rapidity()) :
                                                       std::abs(ptop1.p4.Rapidity()) - std::abs(ptop0.p4.Rapidity());
            bool correct_dy_sign = true_dy*dy>0;

            efficiencyBase::fill("dilution_incl",       dy, correct_dy_sign );
            efficiencyBase::fill("dilution_mtt",       mtt, correct_dy_sign );
            efficiencyBase::fill("dilution_ytt",       ytt, correct_dy_sign );
            efficiencyBase::fill("dilution_pttt",     pttt, correct_dy_sign );
            efficiencyBase::fill("dilution_betatt", betatt, correct_dy_sign );
        }
    }

    return;
}

// THE END

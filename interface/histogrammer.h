#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <string>
#include <map>
#include <vector>

#include "Analysis/CyMiniAna/interface/histogrammerBase.h"
#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/Event.h"

class histogrammer : public histogrammerBase {
  public:

    histogrammer( configuration& cmaConfig, std::string name="" );

    // Default - so we can clean up;
    virtual ~histogrammer();

    /* initialize histograms (1D, 2D, & 3D) */
    virtual void initialize( TFile& outputFile, bool doSystWeights=false );
    virtual void bookHists( std::string name );

    void init_hists_4vec(const std::string& name, const std::string& prefix);
    void init_hists_leptons(const std::string& name, const std::string& prefix);
    void init_hists_jets(const std::string& name, const std::string& prefix);
    void init_hists_ljets(const std::string& name, const std::string& prefix);
    void init_hists_ttbarAC(const std::string& name);

    /* fill histograms */
    virtual void fill( Event& event, const std::vector<unsigned int>& evtsel_decisions=std::vector<unsigned int>() );
    virtual void fill( const std::string& name, Event& event, double event_weight );

  protected:

    bool m_isMC;
    bool m_doSystWeights;

    std::vector<std::string> m_lepton_charges = {"Qpos","Qneg"};
    std::vector<std::string> m_containments;
    std::map<int,std::string> m_mapContainmentRev;

    bool m_useJets;
    bool m_useLjets;
    bool m_useLeptons;
    bool m_useNeutrinos;
    bool m_isTtbar;
    bool m_isOneLeptonAnalysis;
};

#endif

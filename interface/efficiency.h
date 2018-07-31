#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TEfficiency.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <string>
#include <map>
#include <vector>

#include "Analysis/CyMiniAna/interface/efficiencyBase.h"
#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/Event.h"
#include "Analysis/CyMiniAna/interface/configuration.h"

class efficiency : public efficiencyBase {
  public:

    // Default
    efficiency(configuration &cmaConfig);

    // Default - so we can clean up;
    virtual ~efficiency();

    /* Book efficiencies */
    virtual void initialize( TFile& outputFile );

    /* fill efficiencies */
    virtual void fill( Event &event, const std::vector<unsigned int>& evtsel_decisions=std::vector<unsigned int>());

  protected:

    bool m_isTtbar;
    bool m_isOneLeptonAnalysis;
};

#endif

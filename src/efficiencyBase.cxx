/*
Created:        30 July 2018
Last Updated:   30 July 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Make TEfficiencies for measuring efficiencies
https://root.cern.ch/doc/master/classTEfficiency.html
*/
#include "Analysis/CyMiniAna/interface/efficiencyBase.h"


efficiencyBase::efficiencyBase(configuration &cmaConfig) : 
  m_config(&cmaConfig){
   m_map_efficiencies.clear();
  }

efficiencyBase::~efficiencyBase() {}


// -- Existing efficiencies
void efficiencyBase::init_eff(const TH1 &passed, const TH1 &total){
    /* Initialize efficiency -- existing efficiencies */
    std::string name(total.GetName());
    name+="_clone";
    m_map_efficiencies[name] = new TEfficiency(passed,total);

    return;
} 

// -- 1D efficiencies
void efficiencyBase::init_eff( const std::string &name, const unsigned int nBins, const double x_min, const double x_max ){
    /* Initialize efficiency -- equal bins */
    m_map_efficiencies[name] = new TEfficiency((name).c_str(), (name).c_str(),nBins,x_min,x_max);
//    m_map_efficiencies[name]->SetUseWeightedEvents();

    return;
}

void efficiencyBase::init_eff( const std::string &name, const unsigned int nBins, const double *xbins ){
    /* Initialize efficiency -- variable bins */
    m_map_efficiencies[name] = new TEfficiency((name).c_str(), (name).c_str(),nBins,xbins);
//    m_map_efficiencies[name]->SetUseWeightedEvents();

    return;
}

// -- 2D efficiencies
void efficiencyBase::init_eff( const std::string &name, const unsigned int nBinsX, const double x_min, const double x_max,
                              const unsigned int nBinsY, const double y_min, const double y_max ){
    /* Initialize efficiency -- equal bins */
    m_map_efficiencies[name] = new TEfficiency((name).c_str(), (name).c_str(),
                                               nBinsX,x_min,x_max,nBinsY,y_min,y_max);
//    m_map_efficiencies[name]->SetUseWeightedEvents();

    return;
}

void efficiencyBase::init_eff( const std::string &name, const unsigned int nBinsX, const double *xbins,
                              const unsigned int nBinsY, const double *ybins ){
    /* Initialize efficiency -- variable bins */
    m_map_efficiencies[name] = new TEfficiency((name).c_str(), (name).c_str(),
                                               nBinsX,xbins,nBinsY,ybins);
//    m_map_efficiencies[name]->SetUseWeightedEvents();

    return;
}




void efficiencyBase::fill( const std::string &name, const double &value, const bool &decision ){
    /* Fill efficiencies with values! */
    m_map_efficiencies.at(name)->Fill(decision, value);

    return;
}

void efficiencyBase::fill( const std::string &name, 
                         const double &xvalue, const double &yvalue, const bool &decision ){
    /* Fill efficiencies with values! */
    m_map_efficiencies.at(name)->Fill(decision, xvalue, yvalue);

    return;
}

// THE END

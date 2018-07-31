#ifndef HISTOGRAMMERBASE_H
#define HISTOGRAMMERBASE_H

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

#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/Event.h"

class histogrammerBase {
  public:

    // Default - so root can load based on a name;
    histogrammerBase( configuration& cmaConfig, std::string name="" );

    // Default - so we can clean up;
    virtual ~histogrammerBase();

    /* initialize histograms (1D, 2D, & 3D) */
    virtual void init_hist( const std::string& name, 
                                  const unsigned int nBins, const double x_min, const double x_max );
    virtual void init_hist( const std::string& name, 
                                  const unsigned int nBins, const double *xbins );

    virtual void init_hist( const std::string& name, 
                                  const unsigned int nBinsX, const double x_min, const double x_max,
                                  const unsigned int nBinsY, const double y_min, const double y_max );
    virtual void init_hist( const std::string& name, 
                                  const unsigned int nBinsX, const double *xbins,
                                  const unsigned int nBinsY, const double *ybins );
    virtual void init_hist( const std::string& name, const unsigned int nBinsX, const double x_min, const double x_max,
                              const unsigned int nBinsY, const double y_min, const double y_max,
                              const unsigned int nBinsZ, const double z_min, const double z_max );
    virtual void init_hist( const std::string& name, const unsigned int nBinsX, const double *xbins,
                              const unsigned int nBinsY, const double *ybins,  
                              const unsigned int nBinsZ, const double *zbins );

    /* fill histograms */
    virtual void fill( const std::string& name, const double& value, const double& weight );
    virtual void fill( const std::string& name, const double& xvalue, const double& yvalue, const double& weight );
    virtual void fill( const std::string& name, const double& xvalue, const double& yvalue, const double& zvalue, const double& weight );

    /* Put over/underflow in last/first bins.  Called from outside macro */
    void overUnderFlow();
    virtual void overFlow();
    virtual void underFlow();

  protected:

    configuration *m_config;
    std::string m_name;

    std::map<std::string, TH1D*> m_map_histograms1D;
    std::map<std::string, TH2D*> m_map_histograms2D;
    std::map<std::string, TH3D*> m_map_histograms3D;

    std::vector<std::string> m_names;
    bool m_putOverflowInLastBin;
    bool m_putUnderflowInFirstBin;
};

#endif

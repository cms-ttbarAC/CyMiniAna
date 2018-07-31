/*
Created:        --
Last Updated:   26 July 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Make histograms for systematic uncertainties (& nominal) 
to go into plots || TRexFitter

*/
#include "Analysis/CyMiniAna/interface/histogrammer.h"


histogrammer::histogrammer( configuration& cmaConfig, std::string name ) :
  histogrammerBase::histogrammerBase(cmaConfig,name){
    m_isMC  = m_config->isMC();

    m_useJets      = m_config->useJets();
    m_useLjets     = m_config->useLargeRJets();
    m_useLeptons   = m_config->useLeptons();
    m_useNeutrinos = m_config->useNeutrinos();

    m_isTtbar = m_config->isTtbar();
    m_isOneLeptonAnalysis = m_config->isOneLeptonAnalysis();

    m_mapContainmentRev = m_config->mapOfPartonContainmentRev();
    m_containments.clear();
    for (const auto& x : m_config->mapOfPartonContainment())
        m_containments.push_back(x.first);
  }

histogrammer::~histogrammer() {}


/**** INITIALIZE HISTOGRAMS ****/
void histogrammer::initialize( TFile& outputFile, bool doSystWeights ){
    /* Setup some values and book histograms */
    m_doSystWeights = doSystWeights;
    outputFile.cd();


    // loop over selections (typically only one treename)
    for (const auto& sel : m_config->selections() ){
        bookHists( m_name+sel );

        // weight systematics
        if (m_isMC && m_doSystWeights){
            for (const auto& syst : m_config->listOfWeightSystematics()){
                bookHists( m_name+sel+"_"+syst );
            } // end weight systematics

            // vector weight systematics
            for (const auto& syst : m_config->mapOfWeightVectorSystematics()){
                for (unsigned int el=0;el<syst.second;++el){
                    std::string weightIndex = std::to_string(el);
                    bookHists( m_name+sel+"_"+weightIndex+"_"+syst.first );
                } // end components of vector
            } // end vector weight systematics
        } // end if MC and save weight systematics
    } // end loop over selections

    return;
}


void histogrammer::init_hists_4vec(const std::string& name, const std::string& prefix){
    /* Histograms for basic four-vector values */
    init_hist(prefix+"_pt_"+name,  2000, 0.0,2000);
    init_hist(prefix+"_eta_"+name,   50,-2.5, 2.5);
    init_hist(prefix+"_phi_"+name,   64,-3.2, 3.2);

    init_hist(prefix+"_deltaPhi_MET_"+name, 64,-3.2,3.2);

    return;
}

void histogrammer::init_hists_leptons(const std::string& name, const std::string& prefix){
    /* Histograms for lepton distributions */
    init_hists_4vec(name,prefix);

    init_hist(prefix+"_charge_"+name,240,-1.2, 1.2);
    init_hist(prefix+"_ptrel_"+name, 500, 0.0, 500);
    init_hist(prefix+"_drmin_"+name,  50, 0.0,   5);

    return;
}


void histogrammer::init_hists_jets(const std::string& name, const std::string& prefix){
    /* Histograms for AK4 distributions */
    init_hists_4vec(name,prefix);

    init_hist(prefix+"_bdisc_"+name, 100,  0.0,    1.0);

    return;
}


void histogrammer::init_hists_ljets(const std::string& name, const std::string& prefix){
    /* Histograms for AK8 distributions -- put in function to avoid duplicating init_hist calls*/
    init_hists_4vec(name,prefix);

    init_hist(prefix+"_SDmass_"+name,  500,  0.0,  500.0);
    init_hist(prefix+"_charge_"+name, 1000, -5.0,    5.0);

    init_hist(prefix+"_BEST_t_"+name,  100,  0.0,    1.0);
    init_hist(prefix+"_BEST_w_"+name,  100,  0.0,    1.0);
    init_hist(prefix+"_BEST_z_"+name,  100,  0.0,    1.0);
    init_hist(prefix+"_BEST_h_"+name,  100,  0.0,    1.0);
    init_hist(prefix+"_BEST_j_"+name,  100,  0.0,    1.0);
    init_hist(prefix+"_BEST_t_j_"+name,100,  0.0,    1.0);

    init_hist(prefix+"_tau1_"+name,    200,  0.0,    2.0);
    init_hist(prefix+"_tau2_"+name,    200,  0.0,    2.0);
    init_hist(prefix+"_tau3_"+name,    200,  0.0,    2.0);
    init_hist(prefix+"_tau21_"+name,   100,  0.0,    1.0);
    init_hist(prefix+"_tau32_"+name,   100,  0.0,    1.0);

    init_hist(prefix+"_subjet0_bdisc_"+name, 100, 0.0, 1.0);
    init_hist(prefix+"_subjet0_charge_"+name,1000, -5.0, 5.0);
    init_hist(prefix+"_subjet0_tau21_"+name,100, 0.0, 1.0);
    init_hist(prefix+"_subjet0_tau32_"+name,100, 0.0, 1.0);
    init_hist(prefix+"_subjet1_bdisc_"+name, 100, 0.0, 1.0);
    init_hist(prefix+"_subjet1_charge_"+name,1000, -5.0, 5.0);
    init_hist(prefix+"_subjet1_tau21_"+name,100, 0.0, 1.0);
    init_hist(prefix+"_subjet1_tau32_"+name,100, 0.0, 1.0);

    init_hist(prefix+"_subjet-b_subjet-w_tau21_"+name,100, 0.0,1.0, 100,0.0,1.0);  // subjet b (=higher bdisc) tau21 vs subjet w (=lower bdisc) tau21
    init_hist(prefix+"_subjet0_subjet1_tau21_"+name,100, 0.0,1.0, 100,0.0,1.0);  // subjet0 tau21 vs subjet1 tau21
    init_hist(prefix+"_subjet0_subjet1_tau32_"+name,100, 0.0,1.0, 100,0.0,1.0);  // subjet0 tau32 vs subjet1 tau32
    init_hist(prefix+"_subjet0_charge_bdisc_"+name,1000,-5.0,5.0, 100,0.0,1.0);  // charge vs bdisc (charge=x-axis)
    init_hist(prefix+"_subjet1_charge_bdisc_"+name,1000,-5.0,5.0, 100,0.0,1.0);  // charge vs bdisc (charge=x-axis)

    init_hist(prefix+"_pt_eta_"+name,       200, 0.0,2000,  50,-2.5, 2.5);    // pt vs eta (pt=x-axis)
    init_hist(prefix+"_pt_SDmass_"+name,    200, 0.0,2000,  50, 0.0, 500);    // pt vs SDmass (pt=x-axis)
    init_hist(prefix+"_pt_BEST_t_"+name,    200, 0.0,2000, 100, 0.0, 1.0);    // pt vs BEST(t) (pt=x-axis)
    init_hist(prefix+"_SDmass_tau32_"+name, 500, 0.0, 500, 100, 0.0,   1);    // SDmass  vs tau32  (SDmass=x-axis)
    init_hist(prefix+"_BEST_t_SDmass_"+name,100, 0.0, 1.0, 500, 0.0, 500);    // BEST(t) vs SDmass (BEST=x-axis)
    init_hist(prefix+"_BEST_t_tau32_"+name, 100, 0.0, 1.0, 100, 0.0,   1);    // BEST(t) vs tau32  (BEST=x-axis)
    init_hist(prefix+"_BEST_t_BEST_j_"+name,100, 0.0, 1.0, 100, 0.0,   1);    // BEST(t) vs tau32  (BEST=x-axis)

    return;
}


void histogrammer::init_hists_ttbarAC(const std::string& name){
    /* Histograms for the ttbar system and AC measurement */
    init_hist("deltay_"+name,  1000,-5.0,  5.0);
    init_hist("mttbar_"+name,  5000, 0.0, 5000);
    init_hist("pTttbar_"+name,  300, 0.0,  600);
    init_hist("yttbar_"+name,   100,  0.,   10);
    init_hist("betattbar_"+name,100,  0.,    1);
    init_hist("deltaR_ttbar_"+name, 50, 0.0, 5.0);   // deltaR between leptonic top and hadronic top

    init_hist("mttbar_deltay_"+name,  5000, 0.0, 5000, 1000,-5.0,  5.0);  // m_ttbar = x-axis; deltay = y-axis
    init_hist("pTttbar_deltay_"+name,  300, 0.0,  600, 1000,-5.0,  5.0);
    init_hist("yttbar_deltay_"+name,   100,  0.,   10, 1000,-5.0,  5.0);
    init_hist("betattbar_deltay_"+name, 100, 0.0, 1.0, 1000,-5.0,  5.0);

    // resolution
    if (m_isTtbar){
        init_hist("resmat_"+name,          100,-5.0,  5.0, 100,-5.0,  5.0);  // reco=x-axis; truth=y-axis
        init_hist("deltay_dyres_"+name,   1000,-5.0,  5.0, 200,-10,10);
        init_hist("mttbar_dyres_"+name,   5000, 0.0, 5000, 200,-10,10);
        init_hist("pTttbar_dyres_"+name,   300, 0.0,  600, 200,-10,10);
        init_hist("betattbar_dyres_"+name, 100,  0.,    1, 200,-10,10);
        init_hist("yttbar_dyres_"+name,    100,  0.,   10, 200,-10,10);
    }

    return;
}


void histogrammer::bookHists( std::string name ){
    /* 
      Book histograms -- modify/inherit this function for analysis-specific hists 

      @param name   This is the string used to identify histograms for different systematics/event weights
    */
    m_names.resize(0); // append names to this to keep track of later

    if (m_useJets){
        init_hist("n_jets_"+name,   31, -0.5,  30.5);
        init_hist("n_btags_"+name,  11, -0.5,  10.5);

        init_hists_jets(name,"jet");   // all jets
        init_hists_jets(name,"jet0");  // leading jet 
        init_hists_jets(name,"jet1");  // sub-leading jet
    }

    if (m_useLjets){
        init_hist("n_ljets_"+name,  31, -0.5,   30.5);              // Number of good AK8 jets
        if (!m_isOneLeptonAnalysis) init_hists_ljets(name,"ljet");  // all AK8 jets
    }

    if (m_useLeptons){
        init_hist("n_el_"+name, 10, 0, 10);
        init_hists_leptons(name,"el");

        init_hist("n_mu_"+name, 10, 0, 10);
        init_hists_leptons(name,"mu");

        if (m_useNeutrinos){
            init_hist("deltaR_mu_nu_"+name, 50, 0.0, 5.0);
            init_hist("deltaR_el_nu_"+name, 50, 0.0, 5.0);

            init_hist("mass_mu_nu_"+name,  400, 0.0, 400);
            init_hist("mass_el_nu_"+name,  400, 0.0, 400);
        } // end if use neutrinos (& leptons)
    } // end if use leptons

    if (m_useNeutrinos){
        init_hists_4vec(name,"nu");
    }

    // kinematics
    init_hist("met_met_"+name, 500,  0.0,  500);
    init_hist("met_phi_"+name,  64, -3.2,  3.2);
    init_hist("ht_"+name,     5000,  0.0, 5000);
    init_hist("mtw_"+name,     500,  0.0,  500);


    // ttbar reconstruction + asymmetry values
    if (m_isOneLeptonAnalysis){
        init_hist("deltaR_lep_ak4_"+name, 50, 0.0, 5.0);
        init_hist("pTrel_lep_ak4_"+name, 100, 0.0, 500);
        init_hist("deltaR_lep_ak8_"+name, 50, 0.0, 5.0);
        init_hist("deltaR_ak4_ak8_"+name, 50, 0.0, 5.0);
        init_hist("deltaR_pTrel_lep_ak4_"+name, 50,0.0,5.0, 100,0.0,500);

        init_hist("el_MET_triangle_"+name,  200, 0.0, 20);
        init_hist("mu_MET_triangle_"+name,  200, 0.0, 20);
        init_hist("ak4_MET_triangle_"+name, 200, 0.0, 20);

        // sjet (selected AK4 jet) histograms
        init_hists_jets(name,"sjet");

        // leptonic top system
        init_hist("leptop_mass_"+name, 500, 0.0, 500);
        init_hist("leptop_pt_"+name,  3000, 0.0,3000);
        init_hist("leptop_eta_"+name,   60,-3.0, 3.0);

        // hadronic top system
        if (m_isTtbar){
            // different containment levels -- combine to get 'total' ttbar 
            for (const auto& c : m_containments){
                std::string cname = c+"_"+name;
                init_hists_ljets(cname,"hadtop");
            } // end loop over containments
        } // end is ttbar
        else
            init_hists_ljets(name,"hadtop");


        // plots of subjet charge/bdisc for different lepton charges
        if (m_useLeptons && m_isOneLeptonAnalysis){
            // All large-R jets that pass basic object selection
            for (const auto& q : m_lepton_charges){
                std::string qname = q+"_"+name;
                init_hist("ljet_charge_"+qname,         1000, -5.0, 5.0);

                init_hist("ljet_subjet_0_charge_"+qname,1000, -5.0, 5.0);
                init_hist("ljet_subjet_0_bdisc_"+qname,  100,  0.0, 1.0);

                init_hist("ljet_subjet_1_charge_"+qname,1000, -5.0, 5.0);
                init_hist("ljet_subjet_1_bdisc_"+qname,  100,  0.0, 1.0);

                init_hist("deltay_"+qname,  1000,-5.0,  5.0);
            }
        }


        init_hists_ttbarAC(name);    // histograms of ttbar system for AC measurement
    } // end if is single lepton analysis

    return;
}




/**** FILL HISTOGRAMS ****/
void histogrammer::fill( Event& event, const std::vector<unsigned int>& evtsel_decisions ){
    /* Fill histograms -- fill histograms based on selection, tree, or systematic weights ("nominal" but different weight)
       This is the function to modify / inherit for analysis-specific purposes
    */
    double event_weight = event.nominal_weight();

    std::vector<std::string> selections = m_config->selections();
    for (unsigned int ss=0, size=selections.size(); ss<size; ss++){
        std::string sel( selections.at(ss) );
        if (!evtsel_decisions.at(ss)) continue;
        fill( m_name+sel, event, event_weight );

        // if there are systematics stored as weights (e.g., b-tagging, pileup, etc.)
        // the following calls the fill() function with different event weights
        // to make histograms
        bool isNominal = m_config->isNominalTree( event.treeName() );
        if (m_isMC && isNominal && m_doSystWeights){
            // weight systematics
            event_weight = 1.0;
            for (const auto& syst : m_config->listOfWeightSystematics()){
                event_weight = event.getSystEventWeight( syst );
                fill( m_name+sel+"_"+syst, event, event_weight );
            } // end weight systematics

            // vector weight systematics
            event_weight = 1.0;
            for (const auto& syst : m_config->mapOfWeightVectorSystematics()){
                for (unsigned int el=0;el<syst.second;++el){
                    event_weight = event.getSystEventWeight( syst.first, el );
                    std::string weightIndex = std::to_string(el);

                    fill( m_name+sel+"_"+weightIndex+"_"+syst.first, event, event_weight );
                } // end components of vector
            } // end vector weight systematics
        } // end if nominal and doSystWeights
    } // end loop over selections

    return;
}


void histogrammer::fill( const std::string& name, Event& event, double event_weight){
    /* Fill histograms -- just use information from the event and fill histogram
       This is the function to modify / inherit for analysis-specific purposes

       //std::vector<Muon> muons = event.muons();  -- merged into "Lepton"
       //std::vector<Electron> electrons = event.electrons();
    */
    cma::DEBUG("HISTOGRAMMER : Fill histograms.");
    cma::DEBUG("HISTOGRAMMER : event weight = "+std::to_string(event_weight) );

    // physics information
    std::vector<Jet> jets = event.jets();
    std::vector<Ljet> ljets = event.ljets();
    std::vector<Lepton> leptons = event.leptons();
    std::vector<Neutrino> neutrinos = event.neutrinos();
    MET met = event.met();

    // fill histograms!
    // -- only filling the "good" objects

    if (m_useJets){
        cma::DEBUG("HISTOGRAMMER : Fill small-R jets");
        histogrammerBase::fill("n_btags_"+name, event.btag_jets().size(), event_weight );

        histogrammerBase::fill("n_jets_"+name, jets.size(), event_weight );

        for (const auto& jet : jets){
            if (!jet.isGood) continue;
            histogrammerBase::fill("jet_pt_"+name,  jet.p4.Pt(),   event_weight);
            histogrammerBase::fill("jet_eta_"+name, jet.p4.Eta(),  event_weight);
            histogrammerBase::fill("jet_phi_"+name, jet.p4.Phi(),  event_weight);
            histogrammerBase::fill("jet_bdisc_"+name, jet.bdisc,  event_weight);
        }
    }


    if (m_useLjets){
        cma::DEBUG("HISTOGRAMMER : Fill large-R jets");
        histogrammerBase::fill("n_ljets_"+name, ljets.size(), event_weight );

        for (const auto& ljet : ljets){
            if (!ljet.isGood) continue;
            histogrammerBase::fill("ljet_pt_"+name,    ljet.p4.Pt(),  event_weight);
            histogrammerBase::fill("ljet_eta_"+name,   ljet.p4.Eta(), event_weight);
            histogrammerBase::fill("ljet_phi_"+name,   ljet.p4.Phi(), event_weight);
            histogrammerBase::fill("ljet_SDmass_"+name,ljet.softDropMass, event_weight);
            histogrammerBase::fill("ljet_charge_"+name,ljet.charge,event_weight);

            histogrammerBase::fill("ljet_BEST_t_"+name,  ljet.BEST_t,  event_weight);
            histogrammerBase::fill("ljet_BEST_w_"+name,  ljet.BEST_w,  event_weight);
            histogrammerBase::fill("ljet_BEST_z_"+name,  ljet.BEST_z,  event_weight);
            histogrammerBase::fill("ljet_BEST_h_"+name,  ljet.BEST_h,  event_weight);
            histogrammerBase::fill("ljet_BEST_j_"+name,  ljet.BEST_j,  event_weight);
            histogrammerBase::fill("ljet_BEST_t_j_"+name,  ljet.BEST_t / (ljet.BEST_t+ljet.BEST_j),  event_weight);

            histogrammerBase::fill("ljet_tau1_"+name,  ljet.tau1,  event_weight);
            histogrammerBase::fill("ljet_tau2_"+name,  ljet.tau2,  event_weight);
            histogrammerBase::fill("ljet_tau3_"+name,  ljet.tau3,  event_weight);
            histogrammerBase::fill("ljet_tau21_"+name, ljet.tau21, event_weight);
            histogrammerBase::fill("ljet_tau32_"+name, ljet.tau32, event_weight);
            histogrammerBase::fill("ljet_subjet0_bdisc_"+name, ljet.subjet0_bdisc, event_weight);
            histogrammerBase::fill("ljet_subjet1_bdisc_"+name, ljet.subjet1_bdisc, event_weight);
            histogrammerBase::fill("ljet_subjet0_charge_"+name,ljet.subjet0_charge,event_weight);
            histogrammerBase::fill("ljet_subjet1_charge_"+name,ljet.subjet1_charge,event_weight);

            float subjet0_tau21 = ljet.subjet0_tau2/ljet.subjet0_tau1;
            float subjet0_tau32 = ljet.subjet0_tau3/ljet.subjet0_tau2;
            float subjet1_tau21 = ljet.subjet1_tau2/ljet.subjet1_tau1;
            float subjet1_tau32 = ljet.subjet1_tau3/ljet.subjet1_tau2;
            histogrammerBase::fill("ljet_subjet0_tau21_"+name, subjet0_tau21, event_weight);
            histogrammerBase::fill("ljet_subjet1_tau21_"+name, subjet1_tau21, event_weight);
            histogrammerBase::fill("ljet_subjet0_tau32_"+name, subjet0_tau32, event_weight);
            histogrammerBase::fill("ljet_subjet1_tau32_"+name, subjet1_tau32, event_weight);

            // fill based on which subjet has a higher b-disc value
            if (ljet.subjet0_bdisc>ljet.subjet1_bdisc)
                histogrammerBase::fill("ljet_subjet-b_subjet-w_tau21_"+name, subjet0_tau21, subjet1_tau21, event_weight);
            else
                histogrammerBase::fill("ljet_subjet-b_subjet-w_tau21_"+name, subjet1_tau21, subjet0_tau21, event_weight);

            histogrammerBase::fill("ljet_subjet0_subjet1_tau21_"+name,subjet0_tau21,subjet1_tau21,event_weight);  // subjet0 tau21 vs subjet1 tau21
            histogrammerBase::fill("ljet_subjet0_subjet1_tau32_"+name,subjet0_tau32,subjet1_tau32,event_weight);  // subjet0 tau32 vs subjet1 tau32

            histogrammerBase::fill("ljet_subjet0_charge_bdisc_"+name, ljet.subjet0_charge, ljet.subjet0_bdisc, event_weight);
            histogrammerBase::fill("ljet_subjet1_charge_bdisc_"+name, ljet.subjet1_charge, ljet.subjet1_bdisc, event_weight);

            histogrammerBase::fill("ljet_pt_eta_"+name,   ljet.p4.Pt(), ljet.p4.Eta(), event_weight);
            histogrammerBase::fill("ljet_pt_SDmass_"+name,ljet.p4.Pt(), ljet.softDropMass, event_weight);

            histogrammerBase::fill("ljet_SDmass_tau32_"+name, ljet.softDropMass, ljet.tau32,  event_weight);    // SDmass  vs tau32  (SDmass=x-axis)
            histogrammerBase::fill("ljet_BEST_t_SDmass_"+name,ljet.BEST_t, ljet.softDropMass, event_weight);    // BEST(t) vs SDmass (BEST=x-axis)
            histogrammerBase::fill("ljet_BEST_t_tau32_"+name, ljet.BEST_t, ljet.tau32,        event_weight);    // BEST(t) vs tau32  (BEST=x-axis)

            int charge(-999);
            if (m_useLeptons && leptons.size()>0) {
                // only interested in these plots for lepton+jets channel (1 lepton reconstructed)
                charge = leptons.at(0).charge;

                if (charge>0) {
                    histogrammerBase::fill("ljet_charge_Qpos_"+name, ljet.charge, event_weight);
                    histogrammerBase::fill("ljet_subjet_0_charge_Qpos_"+name, ljet.subjet0_charge,event_weight);
                    histogrammerBase::fill("ljet_subjet_0_bdisc_Qpos_"+name,  ljet.subjet0_bdisc, event_weight);
                    histogrammerBase::fill("ljet_subjet_1_charge_Qpos_"+name, ljet.subjet1_charge,event_weight);
                    histogrammerBase::fill("ljet_subjet_1_bdisc_Qpos_"+name,  ljet.subjet1_bdisc, event_weight);
                }
                else {
                    histogrammerBase::fill("ljet_charge_Qneg_"+name, ljet.charge, event_weight);
                    histogrammerBase::fill("ljet_subjet_0_charge_Qneg_"+name, ljet.subjet0_charge,event_weight);
                    histogrammerBase::fill("ljet_subjet_0_bdisc_Qneg_"+name,  ljet.subjet0_bdisc, event_weight);
                    histogrammerBase::fill("ljet_subjet_1_charge_Qneg_"+name, ljet.subjet1_charge,event_weight);
                    histogrammerBase::fill("ljet_subjet_1_bdisc_Qneg_"+name,  ljet.subjet1_bdisc, event_weight);
                }
            } // end if use leptons

            if (m_isTtbar){
                std::string cname = m_mapContainmentRev[ std::abs(ljet.containment) ]+"_"+name;
                histogrammerBase::fill("ljet_pt_"+cname,     ljet.p4.Pt(),      event_weight);
                histogrammerBase::fill("ljet_SDmass_"+cname, ljet.softDropMass, event_weight);
                histogrammerBase::fill("ljet_BEST_t_"+cname,  ljet.BEST_t,       event_weight);
                histogrammerBase::fill("ljet_BEST_w_"+cname,  ljet.BEST_w,       event_weight);
                histogrammerBase::fill("ljet_BEST_z_"+cname,  ljet.BEST_z,       event_weight);
                histogrammerBase::fill("ljet_BEST_h_"+cname,  ljet.BEST_h,       event_weight);
                histogrammerBase::fill("ljet_BEST_j_"+cname,  ljet.BEST_j,       event_weight);
                histogrammerBase::fill("ljet_BEST_t_j_"+cname,ljet.BEST_t / (ljet.BEST_t+ljet.BEST_j), event_weight);
                histogrammerBase::fill("ljet_tau21_"+cname, ljet.tau21, event_weight);
                histogrammerBase::fill("ljet_tau32_"+cname, ljet.tau32, event_weight);

                histogrammerBase::fill("ljet_pt_SDmass_"+cname, ljet.p4.Pt(),ljet.softDropMass,event_weight);

                histogrammerBase::fill("ljet_SDmass_tau32_"+cname, ljet.softDropMass, ljet.tau32,  event_weight);    // SDmass  vs tau32  (SDmass=x-axis)
                histogrammerBase::fill("ljet_BEST_t_SDmass_"+cname,ljet.BEST_t, ljet.softDropMass, event_weight);    // BEST(t) vs SDmass (BEST=x-axis)
                histogrammerBase::fill("ljet_BEST_t_tau32_"+cname, ljet.BEST_t, ljet.tau32,        event_weight);    // BEST(t) vs tau32  (BEST=x-axis)
            } // end if isttbar
        } // end loop over ljets
    } // end if use ljets

    if (m_useLeptons){
        cma::DEBUG("HISTOGRAMMER : Fill leptons");
        unsigned int n_electrons(0);
        unsigned int n_muons(0);
        for (const auto& lep : leptons){
            if (!lep.isGood) continue;

            if (lep.isElectron){
                histogrammerBase::fill("el_pt_"+name,  lep.p4.Pt(),  event_weight);
                histogrammerBase::fill("el_eta_"+name, lep.p4.Eta(), event_weight);
                histogrammerBase::fill("el_phi_"+name, lep.p4.Phi(), event_weight);
                histogrammerBase::fill("el_charge_"+name, lep.charge, event_weight);
                histogrammerBase::fill("el_ptrel_"+name,  lep.ptrel,  event_weight);
                histogrammerBase::fill("el_drmin_"+name,  lep.drmin,  event_weight);
                n_electrons++;
            }
            else if (lep.isMuon){
                histogrammerBase::fill("mu_pt_"+name,  lep.p4.Pt(),  event_weight);
                histogrammerBase::fill("mu_eta_"+name, lep.p4.Eta(), event_weight);
                histogrammerBase::fill("mu_phi_"+name, lep.p4.Phi(), event_weight);
                histogrammerBase::fill("mu_charge_"+name, lep.charge, event_weight);
                histogrammerBase::fill("mu_ptrel_"+name,  lep.ptrel,  event_weight);
                histogrammerBase::fill("mu_drmin_"+name,  lep.drmin,  event_weight);
                n_muons++;
            }
        } // end loop over leptons
        histogrammerBase::fill("n_el_"+name,  n_electrons, event_weight);
        histogrammerBase::fill("n_mu_"+name,  n_muons,     event_weight);
    } // end if use leptons

    if (m_useNeutrinos){
        cma::DEBUG("HISTOGRAMMER : Fill neutrinos");
        for (const auto& nu : neutrinos){
            histogrammerBase::fill("nu_pt_"+name,  nu.p4.Pt(),  event_weight);
            histogrammerBase::fill("nu_eta_"+name, nu.p4.Eta(), event_weight);
            histogrammerBase::fill("nu_phi_"+name, nu.p4.Phi(), event_weight);
        }
    }

    // kinematics
    cma::DEBUG("HISTOGRAMMER : Fill kinematics");
    histogrammerBase::fill("met_met_"+name, met.p4.Pt(),  event_weight);
    histogrammerBase::fill("met_phi_"+name, met.p4.Phi(), event_weight);
    histogrammerBase::fill("ht_"+name,      event.HT(),   event_weight);
    histogrammerBase::fill("mtw_"+name,     met.mtw,      event_weight);


    if (m_isOneLeptonAnalysis){
        // have the asymmetry readily available in boosted single lepton
        cma::DEBUG("HISTOGRAMMER : Fill ttbar AC values");
        Ttbar1L ttbarSL = event.ttbar1L();

        Jet tt_jet     = ttbarSL.jet;
        Ljet tt_ljet   = ttbarSL.ljet;
        Lepton tt_lep  = ttbarSL.lepton;
        Neutrino tt_nu = ttbarSL.neutrino;

        TLorentzVector top_lep;
        TLorentzVector top_had;
        TLorentzVector ttbar;
        top_lep = ttbarSL.leptop;
        top_had = tt_ljet.p4;
        ttbar   = top_had+top_lep;

        float dr_lep_ak4    = tt_jet.p4.DeltaR( tt_lep.p4 );
        float ptrel_lep_ak4 = cma::ptrel( tt_lep.p4, tt_jet.p4);

        histogrammerBase::fill("deltaR_lep_ak4_"+name, dr_lep_ak4,    event_weight);
        histogrammerBase::fill("pTrel_lep_ak4_"+name,  ptrel_lep_ak4, event_weight);
        histogrammerBase::fill("deltaR_lep_ak8_"+name, top_had.DeltaR( tt_lep.p4 ), event_weight);
        histogrammerBase::fill("deltaR_ak4_ak8_"+name, top_had.DeltaR( tt_jet.p4 ), event_weight);
        histogrammerBase::fill("deltaR_pTrel_lep_ak4_"+name, dr_lep_ak4, ptrel_lep_ak4, event_weight);

        histogrammerBase::fill("hadtop_pt_"+name,     top_had.Pt(), event_weight);
        histogrammerBase::fill("hadtop_eta_"+name,    top_had.Eta(), event_weight);
        histogrammerBase::fill("hadtop_SDmass_"+name, tt_ljet.softDropMass, event_weight);
        histogrammerBase::fill("hadtop_tau21_"+name,  tt_ljet.tau21, event_weight);
        histogrammerBase::fill("hadtop_tau32_"+name,  tt_ljet.tau32, event_weight);

        // asymmetry
        float dy     = ttbarSL.dy;
        float mtt    = ttbar.M();
        float pttt   = ttbar.Pt();
        float ytt    = std::abs(ttbar.Rapidity());
        float betatt = std::abs(top_had.Pz() + top_lep.Pz()) / (top_had.E() + top_lep.E());

        histogrammerBase::fill("deltay_"+name,  dy,   event_weight);
        histogrammerBase::fill("mttbar_"+name,  mtt,  event_weight);
        histogrammerBase::fill("pTttbar_"+name, pttt, event_weight);
        histogrammerBase::fill("yttbar_"+name,  ytt,  event_weight);
        histogrammerBase::fill("betattbar_"+name, betatt, event_weight);

        histogrammerBase::fill("mttbar_deltay_"+name,  mtt,  dy, event_weight);
        histogrammerBase::fill("pTttbar_deltay_"+name, pttt, dy, event_weight);
        histogrammerBase::fill("yttbar_deltay_"+name,  ytt,  dy, event_weight);
        histogrammerBase::fill("betattbar_deltay_"+name, betatt, dy, event_weight);

        if (m_isTtbar){
            float true_dy(0);
            std::vector<TruthTop> truth_tops = event.truth();
            std::vector<Parton> truth_partons = event.truth_partons();
            if (truth_tops.size()==2){
                TruthTop top0 = truth_tops.at(0);
                TruthTop top1 = truth_tops.at(1);
                Parton ptop0  = truth_partons.at( top0.Top );
                Parton ptop1  = truth_partons.at( top1.Top );

                true_dy = (top0.isTop && top1.isAntiTop) ? std::abs(ptop0.p4.Rapidity()) - std::abs(ptop1.p4.Rapidity()) : 
                                                           std::abs(ptop1.p4.Rapidity()) - std::abs(ptop0.p4.Rapidity());
                histogrammerBase::fill("resmat_"+name, dy, true_dy, event_weight);  // x=reco, y=truth

                float dyres = true_dy - dy;
                histogrammerBase::fill("deltay_dyres_"+name,   dy,    dyres, event_weight);
                histogrammerBase::fill("mttbar_dyres_"+name,   mtt,   dyres, event_weight);
                histogrammerBase::fill("pTttbar_dyres_"+name,  pttt,  dyres, event_weight);
                histogrammerBase::fill("betattbar_dyres_"+name,betatt,dyres, event_weight);
                histogrammerBase::fill("yttbar_dyres_"+name,   ytt,   dyres, event_weight);
            }
        } // end response matrix creation
    }

    cma::DEBUG("HISTOGRAMMER : End histograms");

    return;
}

// THE END

cet_build_plugin(SBNEventWeight art::module
  LIBRARIES
  sbnobj::Common_SBNEventWeight sbncode_SBNEventWeight_Base
  sbncode_SBNEventWeight_Calculators_CrossSection
  sbncode_SBNEventWeight_Calculators_BNBFlux nusimdata::SimulationBase
  ROOT::Geom)

cet_build_plugin(SystToolsEventWeight art::module
  LIBRARIES
  sbnobj::Common_SBNEventWeight sbncode_SBNEventWeight_Base
  systtools::interface
  systtools::interpreters
  systtools::utility
  nusyst::interface
  nusyst::systproviders
  nusyst::responsecalculators
  larnusystematics::utility
  nugen::EventGeneratorBase_GENIE
  ${GENIE_LIB_LIST}
  )

install_headers()
install_fhicl()
install_source()


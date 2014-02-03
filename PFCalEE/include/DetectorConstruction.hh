#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "SamplingSection.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <map>
#include <string>

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

/**
   @class DetectorConstruction
   @short builds a simple detector
 */
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  enum DetectorVersion { v_CALICE, v_CALICE_Pb, v_UNIFORM, v_UNIFORM_08, v_UNIFORM_05, v_UNIFORM_03, v_JV,
			 v_UNIFORM_Si60, v_UNIFORM_Si80, v_UNIFORM_Si120, v_UNIFORM_Si200, v_UNIFORM_Si300, v_UNIFORM_Si500, v_VJ,
			 v_CALICE_Pb_Si60, v_CALICE_Pb_Si80, v_CALICE_Pb_Si120, v_CALICE_Pb_Si200, v_CALICE_Pb_Si300, v_CALICE_Pb_Si500,
			 v_HGCALEE, v_HGCALEE_Si500,

			 v_HGCAL_CONCEPT=100,
			 v_HGCAL_CONCEPT_thickSi,
			 v_HGCAL_CONCEPT_thinSi,
			 v_HGCAL_CONCEPT_CALICE,
			 v_HGCAL_CONCEPT_fineSampling,
			 v_HGCAL_CONCEPT_coarseSampling
  };

  /**
     @short CTOR
   */
  DetectorConstruction(G4int ver=DetectorConstruction::v_CALICE);

  /**
     @short calorimeter structure (sampling sections)
   */
  std::vector<SamplingSection> m_caloStruct;
  std::vector<SamplingSection> *getStructure() { return &m_caloStruct; }

  /**
     @short define the calorimeter materials
   */
  void DefineMaterials(); 
  std::map<std::string, G4Material *> m_materials;
  
  /**
     @short set magnetic field
   */
  void SetMagField(G4double fieldValue);
  G4UniformMagField* m_magField;      //pointer to the magnetic field

  /**
     @short DTOR
   */
  ~DetectorConstruction();
  

  /**
     @short getters
   */
  G4double GetCalorSizeYZ() { return m_CalorSizeYZ; }
  G4double GetCalorSizeX()  { return m_CalorSizeX; }
  G4double GetWorldSizeYZ() { return m_WorldSizeYZ; }
  G4double GetWorldSizeX()  { return m_WorldSizeX; }

  /**
     @short build the detector
   */
  G4VPhysicalVolume* Construct();

private:

  //detector version
  int version_;

  /**
     @short compute the calor dimensions
   */
  void UpdateCalorSize(); 

  /**
     @short build the calorimeter
   */
  G4VPhysicalVolume* ConstructCalorimeter();     

  std::vector<G4Material* > m_SensitiveMaterial;
  
  G4double           m_CalorSizeYZ, m_CalorSizeX;
  G4double           m_WorldSizeYZ, m_WorldSizeX;
            
  G4Box*             m_solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   m_logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* m_physWorld;     //pointer to the physical World  
  
  DetectorMessenger* m_detectorMessenger;  //pointer to the Messenger
};


#endif


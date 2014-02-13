#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "HGCSSSimHit.hh"

using namespace std;

//
DetectorConstruction::DetectorConstruction(G4int ver) : version_(ver)
{
  //radiation lengths: cf. http://pdg.lbl.gov/2012/AtomicNuclearProperties/
  //W 0.3504 cm
  //Pb 0.5612 cm

  switch(version_)
    {
      //cf. http://arxiv.org/abs/0805.4833
    case v_CALICE:
      G4cout << "[DetectorConstruction] starting v_CALICE (10x0.4+10x0.8+10x1.2)X_0 with Tungsten" << G4endl;
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(0.4*3.504*mm,0.0*mm,0.525*mm,1.0*mm,2.5*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(0.8*3.504*mm,0.0*mm,0.525*mm,1.0*mm,2.5*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.2*3.504*mm,0.0*mm,0.525*mm,1.0*mm,2.5*mm) );
      break;
    case v_CALICE_Pb: case v_CALICE_Pb_Si60: case v_CALICE_Pb_Si80: case v_CALICE_Pb_Si120: case v_CALICE_Pb_Si200: case v_CALICE_Pb_Si300: case v_CALICE_Pb_Si500:
      {
	float siWidth(0.525);
	if(version_==v_CALICE_Pb_Si60)  siWidth=0.060;
	if(version_==v_CALICE_Pb_Si80)  siWidth=0.080;
	if(version_==v_CALICE_Pb_Si120) siWidth=0.120;
	if(version_==v_CALICE_Pb_Si200) siWidth=0.200;
	if(version_==v_CALICE_Pb_Si300) siWidth=0.300;
	if(version_==v_CALICE_Pb_Si500) siWidth=0.500;
	G4cout << "[DetectorConstruction] starting v_CALICE_Pb (10x0.4+10x0.8+10x1.2)X_0 with Lead Si width=" << siWidth << G4endl;
	for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(0.4*5.612*mm,0.0*mm,siWidth*mm,1.0*mm,1.0*mm) );
	for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(0.8*5.612*mm,0.0*mm,siWidth*mm,1.0*mm,1.0*mm) );
	for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.2*5.612*mm,0.0*mm,siWidth*mm,1.0*mm,1.0*mm) );
	break;
      }
    case v_HGCALEE: case v_HGCALEE_Si500:
      {
	G4cout << "[DetectorConstruction] starting v_HGCALEE" << G4endl;
	float siWidth(0.200);
	if(version_==v_HGCALEE_Si500) siWidth=0.500;
	for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.6*mm,3.0*mm,siWidth*mm,1.0*mm,1.0*mm) );
	for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(3.3*mm,3.0*mm,siWidth*mm,1.0*mm,1.0*mm) );
	for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(5.6*mm,3.0*mm,siWidth*mm,1.0*mm,1.0*mm) );
	break;
      }
    case v_UNIFORM: case v_UNIFORM_Si60 : case v_UNIFORM_Si80: case v_UNIFORM_Si120: case v_UNIFORM_Si200: case v_UNIFORM_Si300: case v_UNIFORM_Si500:
      {
	float siWidth(0.300);
	if(version_==v_UNIFORM_Si60)  siWidth=0.060;
	if(version_==v_UNIFORM_Si80)  siWidth=0.080;
	if(version_==v_UNIFORM_Si200) siWidth=0.20;
	if(version_==v_UNIFORM_Si500) siWidth=0.500;
	G4cout << "[DetectorConstruction] starting v_UNIFORM (26x1)X_0 with Lead and Si width = " << siWidth << " mm" << G4endl;
	for(int i=0; i<26; i++) m_caloStruct.push_back( SamplingSection(5.612*mm,3.0*mm,siWidth*mm,1.0*mm,1.0*mm) );
      }
      break;
    case v_UNIFORM_08:
      G4cout << "[DetectorConstruction] starting v_UNIFORM (33x0.8)X_0 with Lead" << G4endl;
      for(int i=0; i<33; i++) m_caloStruct.push_back( SamplingSection(0.8*5.612*mm,3.0*mm,0.300*mm,1.0*mm,1.0*mm) );
      break;
    case v_UNIFORM_05:
      G4cout << "[DetectorConstruction] starting v_UNIFORM (52x0.5)X_0 with Lead" << G4endl;
      for(int i=0; i<52; i++) m_caloStruct.push_back( SamplingSection(0.54*5.612*mm,3.0*mm,0.300*mm,1.0*mm,1.0*mm) );
      break;
    case v_UNIFORM_03:
      G4cout << "[DetectorConstruction] starting v_UNIFORM (87x0.3)X_0 with Lead" << G4endl;
      for(int i=0; i<87; i++) m_caloStruct.push_back( SamplingSection(0.3*5.612*mm,3.0*mm,0.300*mm,1.0*mm,1.0*mm) );
      break;
    case v_JV:
      G4cout << "[DetectorConstruction] starting v_JV" << G4endl;
      for(int i=0; i<6; i++) m_caloStruct.push_back( SamplingSection(1.0*5.612*mm,3.0*mm,0.300*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<15; i++) m_caloStruct.push_back( SamplingSection(0.66*5.612*mm,3.0*mm,0.300*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.0*5.612*mm,3.0*mm,0.300*mm,1.0*mm,1.0*mm) );
      break;
    case v_VJ:
      G4cout << "[DetectorConstruction] starting v_VJ" << G4endl;
      for(int i=0; i<15; i++) m_caloStruct.push_back( SamplingSection(0.66*5.612*mm,3.0*mm,0.300*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<16; i++) m_caloStruct.push_back( SamplingSection(1.0*5.612*mm,3.0*mm,0.300*mm,1.0*mm,1.0*mm) );
      break;

    case v_HGCAL_CONCEPT:
      G4cout << "[DetectorConstruction] starting concept" << G4endl;
      m_caloStruct.push_back( SamplingSection(0.0*mm,3.0*mm,0.200*mm,2.0*mm,2.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.6*mm,3.0*mm,0.200*mm,2.0*mm,2.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(3.3*mm,3.0*mm,0.200*mm,2.0*mm,2.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(5.6*mm,3.0*mm,0.200*mm,2.0*mm,2.0*mm) );
      break;

    case v_HGCAL_CONCEPT_thickSi:
      G4cout << "[DetectorConstruction] starting concept thick Si" << G4endl;
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.6*mm,3.0*mm,0.525*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(3.3*mm,3.0*mm,0.525*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(5.6*mm,3.0*mm,0.525*mm,1.0*mm,1.0*mm) );
      break;

    case v_HGCAL_CONCEPT_thinSi:
      G4cout << "[DetectorConstruction] starting concept thin Si" << G4endl;
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.6*mm,3.0*mm,0.100*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(3.3*mm,3.0*mm,0.100*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(5.6*mm,3.0*mm,0.100*mm,1.0*mm,1.0*mm) );
      break;

    case v_HGCAL_CONCEPT_CALICE:
      G4cout << "[DetectorConstruction] starting concept calice" << G4endl;
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(0.4*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(0.8*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.3*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      break;

    case v_HGCAL_CONCEPT_fineSampling:
      G4cout << "[DetectorConstruction] starting concept fine sampling" << G4endl;
      for(int i=0; i<6; i++) m_caloStruct.push_back( SamplingSection(0.3*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<6; i++) m_caloStruct.push_back( SamplingSection(0.5*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(0.8*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.2*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      break;

    case v_HGCAL_CONCEPT_coarseSampling:
      G4cout << "[DetectorConstruction] starting concept fine sampling" << G4endl;
      for(int i=0; i<5; i++) m_caloStruct.push_back( SamplingSection(0.5*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<13; i++) m_caloStruct.push_back( SamplingSection(0.8*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.2*5.612*mm,3.0*mm,0.200*mm,1.0*mm,1.0*mm) );
      break;




    }
  DefineMaterials();
  UpdateCalorSize();
  SetMagField(0);
  m_detectorMessenger = new DetectorMessenger(this);
}

//
DetectorConstruction::~DetectorConstruction() { delete m_detectorMessenger;}

//
void DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();
  m_materials["Abs"] = (version_== v_CALICE) ? 
    nistManager->FindOrBuildMaterial("G4_W",false) :
    nistManager->FindOrBuildMaterial("G4_Pb",false);
  m_materials["Cu"] = nistManager->FindOrBuildMaterial("G4_Cu",false); 
  m_materials["Si"] = nistManager->FindOrBuildMaterial("G4_Si",false);
  m_materials["PCB"] = new G4Material("G10",1.700 * g/cm3,4 );
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(14), 1);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(8) , 2);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(6) , 3);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(1) , 3);
  m_materials["Air"]=nistManager->FindOrBuildMaterial("G4_AIR",false);
}

//
void DetectorConstruction::UpdateCalorSize(){  

  m_CalorSizeX=0;
  for(size_t i=0; i<m_caloStruct.size(); i++)
    m_CalorSizeX=m_CalorSizeX+m_caloStruct[i].Total_thick;
  m_CalorSizeYZ=SIZE_X;

  m_WorldSizeX=m_CalorSizeX*1.1;  
  m_WorldSizeYZ=m_CalorSizeYZ*1.1;

  G4cout << "[DetectorConstruction][UpdateCalorSize] X x YZ = " << m_CalorSizeX << " x " << m_CalorSizeYZ << " mm " <<  G4endl;
}

//
G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //clean old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
   
  //detector's World
  m_solidWorld = new G4Box("Wbox",m_WorldSizeX/2,m_WorldSizeYZ/2,m_WorldSizeYZ/2);
  m_logicWorld = new G4LogicalVolume(m_solidWorld, m_materials["Air"], "Wlog");
  m_physWorld = new G4PVPlacement(0, G4ThreeVector(), m_logicWorld, "Wphys", 0, false, 0);

  //build the stack
  char nameBuf[10];
  G4double xOffset(-m_CalorSizeX/2), xOverburden(0.);
  for(size_t i=0; i<m_caloStruct.size(); i++)
    {
      sprintf(nameBuf,"Abs%d",int(i+1)); 
      std::string baseName(nameBuf);
      G4double thick = m_caloStruct[i].Pb_thick;
      if(thick>0){
	G4Box *solid          = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi = new G4LogicalVolume(solid, m_materials["Abs"], baseName+"log");
	m_caloStruct[i].Pb_vol= new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	m_caloStruct[i].Pb_X0 = m_materials["Abs"]->GetRadlen();
	G4VisAttributes *simpleBoxVisAtt= new G4VisAttributes(G4Colour::Gray());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
      }
      
      sprintf(nameBuf,"Cu%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].Cu_thick;
      if(thick>0){
	G4Box *solid = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["Cu"], baseName+"log");
	m_caloStruct[i].Cu_vol = new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	m_caloStruct[i].Cu_X0 = m_materials["Cu"]->GetRadlen();
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Black());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
      }

      sprintf(nameBuf,"Si%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].Si_thick;
      if(thick>0){
	G4Box *solid = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["Si"], baseName+"log");
	m_caloStruct[i].Si_vol = new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	m_caloStruct[i].Si_X0 = m_materials["Si"]->GetRadlen();
	G4cout << "************" << m_caloStruct[i].Si_X0 << " " << m_caloStruct[i].Si_thick << G4endl;
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::White());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
      }

      sprintf(nameBuf,"PCB%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].PCB_thick;
      if(thick>0){
	G4Box *solid = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["PCB"], baseName+"log");
	m_caloStruct[i].PCB_vol = new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	m_caloStruct[i].PCB_X0 = m_materials["PCB"]->GetRadlen();
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Blue());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
      }

      sprintf(nameBuf,"Air%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].Air_thick;
      if(thick>0){
	G4Box *solid = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["Air"], baseName+"log");
	m_caloStruct[i].Air_vol = new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Cyan());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
      }
    }
 
  //                                        
  // Visualization attributes
  //
  m_logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  return m_physWorld;
}

//
void DetectorConstruction::SetMagField(G4double fieldValue)
{

  if(fieldValue<=0) return; 

  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if(m_magField) delete m_magField;                //delete the existing magn field
  m_magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
  fieldMgr->SetDetectorField(m_magField);
  fieldMgr->CreateChordFinder(m_magField);
  fieldMgr->SetDetectorField(m_magField);  
}


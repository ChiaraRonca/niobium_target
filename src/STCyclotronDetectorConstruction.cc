
#include "STCyclotronDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh" 

#include "G4Element.hh" 
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Tubs.hh" 
#include "G4Cons.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Region.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4VisAttributes.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "STCyclotronRun.hh"
#include "STCyclotronSensitiveTarget.hh"
#include "STCyclotronSensitiveFoil.hh"

#include "STCyclotronRunAction.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4LogicalVolumeStore.hh"
STCyclotronDetectorConstruction::STCyclotronDetectorConstruction()
 
{ 
  
}

STCyclotronDetectorConstruction::~STCyclotronDetectorConstruction()
{
 
}

G4VPhysicalVolume* STCyclotronDetectorConstruction::Construct()
{  
  //Initialization of messenger parameters
  fTarget_diameter = 0.95*cm;
  fDensity_target = 1.005 * g/cm3;
  fTarget_thickness = 1.23*cm;
 
  fFoil_thickness = 25.*micrometer;
  fDensity_foil = 8.3*g/cm3;
   

  //Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;


  //Create the world
  G4double world_hx = 1.*m;
  G4double world_hy = 1.*m;
  G4double world_hz = 1.*m;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld 
    = new G4Box("World",
		world_hx, 
		world_hy, 
		world_hz);
  
  fLogicWorld
    = new G4LogicalVolume(solidWorld,
			  world_mat,
			  "World");

  G4VPhysicalVolume* physWorld 
    = new G4PVPlacement(0,                     //no rotation
			G4ThreeVector(),       //at (0,0,0)
			fLogicWorld,            //its logical volume
			"World",               //its name
			0,                     //its mother  volume
			false,                 //no boolean operation
			0);                    //copy number



  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////Create the detector////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  

  //Overall parameters
  G4double startAngle     = 0.*deg;
  G4double spanningAngle  = 360.*deg;
  
  ////////////////////////////////////
  /////////Define materials///////////
  ////////////////////////////////////


  //ALUMINIUM//
 
  G4Material* al = nist->FindOrBuildMaterial("G4_Al");



  



  //HELIUM
  G4double helium_Z, helium_A, helium_density, helium_pressure, helium_temperature;
  helium_Z = 2.;
  helium_A = 4*g/mole;
  helium_density = 0.004002*g/cm3; // with T=298.15 K To modify !
  helium_pressure = 24.82*bar;
  helium_temperature = 298.15*kelvin; //15 to 20°
  G4Material* helium = new G4Material("helium", helium_Z, helium_A, helium_density, kStateGas, helium_temperature, helium_pressure);


  //PLATINIUM
  G4Material* pt =  nist->FindOrBuildMaterial("G4_Pt");
  
  
  //Oxygen 18 water 
   G4Isotope* O16 = new G4Isotope("O16", 8, 16, 15.994 * g/mole);
   G4Isotope* O18 = new G4Isotope("O18", 8, 18, 17.999 * g/mole);
   
   G4Element* O18Element = new G4Element("Oxygen18", "O18", 2);
   O18Element->AddIsotope(O16, 2. * perCent);
   O18Element->AddIsotope(O18, 98. * perCent);
   
   G4Element* H = nist->FindOrBuildElement("H");
   
   G4Material* H2O18 = new G4Material("Water18O", 1.005 * g/cm3, 2);
   H2O18->AddElement(H,2);
   H2O18->AddElement(O18Element,1);
   
   
   //Create vacuum around the beam//
  G4double vacuum_density = 1.e-24*g/cm3;
  
  G4Material* vacuum_beam = new G4Material("vacuumBeam", vacuum_density, 1);
  vacuum_beam->AddElement(H,2);
	
  //fTarget_Material =  nist->FindOrBuildMaterial("G4_Ni");
  fTarget_Material =  H2O18;
  
  G4Material* niobium = nist->FindOrBuildMaterial("G4_Nb");
  

 G4int ncomponents;
 
  //FOIL MATERIAL
  //fFoil_Material =  nist->FindOrBuildMaterial("G4_Co");
   G4Element* Co = nist->FindOrBuildElement("Co");
   G4Element* Cr = nist->FindOrBuildElement("Cr");
   G4Element* Fe = nist->FindOrBuildElement("Fe");
   G4Element* Ni = nist->FindOrBuildElement("Ni");
   G4Element* W = nist->FindOrBuildElement("W");
   G4Element* Mo = nist->FindOrBuildElement("Mo");
   G4Element* Mn = nist->FindOrBuildElement("Mn");
   G4Element* C = nist->FindOrBuildElement("C");
   
 
   
    G4Material* Havar = new G4Material("Havar", fDensity_foil, ncomponents=8);
    Havar->AddElement(Co, 41.78*perCent);
    Havar->AddElement(Cr, 22.29*perCent);
    Havar->AddElement(Fe, 18.11*perCent);
    Havar->AddElement(Ni, 12.83*perCent);
    Havar->AddElement(W, 0.88*perCent);
    Havar->AddElement(Mo, 1.45*perCent);
    Havar->AddElement(Mn, 1.69*perCent);
    Havar->AddElement(C, 0.96*perCent);
    
    fFoil_Material = Havar;
     //PVC
 G4Element* Cl = nist->FindOrBuildElement("Cl");

 G4Material* PVC = new G4Material("PVC", 1.4*g/cm3, ncomponents=3);
 PVC->AddElement(C,2);
 PVC->AddElement(H,3);
  PVC->AddElement(Cl,1);
  
  ////////////////////////////////////////////////////////////////////
  /////////////////////////////PART1//////////////////////////////////
  ////////////////////////////////////////////////////////////////////

//box
 G4Box* rearF1 = new G4Box("aluminumBox", 5.9263 / 2.*cm, 5.9263 / 2.*cm, 8.45 / 2. * cm);
	
//cylindrical hole
G4double hole_radius = 0.875 * cm;
G4double hole_height = 0.5*15.95*cm;
G4double ext_radius = 1.25 * cm;
G4double ext_height = 0.5*15.95*cm;

G4Tubs* solidHole = new G4Tubs("Hole", 0, hole_radius, hole_height, startAngle, spanningAngle);
G4Tubs* extTube = new G4Tubs("Hole", hole_radius, ext_radius, ext_height, startAngle, spanningAngle);

G4LogicalVolume* logicHole = new G4LogicalVolume(solidHole, vacuum_beam, "Hole");
G4LogicalVolume* logicExtTube = new G4LogicalVolume(extTube, al, "tube");
G4VisAttributes* ExtTubeVisAtt = new G4VisAttributes(G4Colour(0.4, 0.5, 0.5));
	
logicExtTube->SetVisAttributes(ExtTubeVisAtt);

G4Tubs* buco = new G4Tubs("buco", 0, ext_radius, 0.5*8.46*cm, startAngle, spanningAngle);


//box with hole
G4SubtractionSolid* boxHole = new G4SubtractionSolid("boxHole", rearF1, buco, 0, G4ThreeVector(0,0,0));


G4LogicalVolume* boxHoleLogical = new G4LogicalVolume(boxHole, PVC, "boxHoleLogical");


G4VPhysicalVolume* boxHolePhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 4.25*cm + 8.45/2.*cm), boxHoleLogical, "boxHolePhysical", fLogicWorld, false, 0, true);

G4VPhysicalVolume* physicalHole = new G4PVPlacement(0, G4ThreeVector(0., 0., 15.95/2.*cm), logicHole, "physicalHole", fLogicWorld, false, 0, true);
G4VPhysicalVolume* physicalTube = new G4PVPlacement(0, G4ThreeVector(0., 0., 15.95/2.*cm), logicExtTube, "physicalTube", fLogicWorld, false, 0, true);

G4Tubs* PVC_cyl = new G4Tubs("PVC_cyl", ext_radius, 2.95*cm, 2.07/2.*cm, startAngle, spanningAngle);
G4LogicalVolume* PVC_cylLogical = new G4LogicalVolume(PVC_cyl, PVC, "PVC_cylLogical");
G4VPhysicalVolume* PVC_cylPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 2.18*cm + 2.07/2.*cm), PVC_cylLogical, "PVC_cylPhysical", fLogicWorld, false, 0, true);


G4Tubs* initial = new G4Tubs("init_cyl", ext_radius, 2.95*cm, 0.3/2.*cm, startAngle, spanningAngle);
G4LogicalVolume* initialLogical = new G4LogicalVolume(initial, al, "init_cylLog");
G4VPhysicalVolume* initialPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 2.18*cm+0.15*cm), initialLogical, "initial_physical", fLogicWorld, false, 0, true);

G4Tubs* ending = new G4Tubs("ending_cyl", ext_radius, 2.95*cm, 0.6/2.*cm, startAngle, spanningAngle);
G4LogicalVolume* endingLogical = new G4LogicalVolume(ending, al, "ending_cylLog");
G4VPhysicalVolume* endingPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 15.95*cm-2*1.65*cm-0.15*cm), initialLogical, "initial_physical", fLogicWorld, false, 0, true);


  ////////////////////////////////////////////////////////////////////
  //////////////////////////aluminum FRONT FLANGE//////////////////////////////
  ////////////////////////////////////////////////////////////////////

G4Box* alFlange = new G4Box("alFlange", 4.637 / 2.*cm, 4.637 / 2.*cm, 1.625 /2. * cm);

G4Tubs* bucoF1 = new G4Tubs("bucoF1", 0, ext_radius, 1.8*cm, startAngle, spanningAngle);

   G4SubtractionSolid* AlFlange = new G4SubtractionSolid("AlFlange", alFlange, bucoF1, 0, G4ThreeVector(0,0,0));
G4LogicalVolume* alFlangeLogical = new G4LogicalVolume(AlFlange, al, "alFlangeLogical");
G4VisAttributes* alFlangeVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
	
alFlangeLogical->SetVisAttributes(alFlangeVisAtt);
G4VPhysicalVolume* alFlangePhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 8.45*cm+2.18*cm+2.07*cm+1.625/2*cm), alFlangeLogical, "alFlangePhysical", fLogicWorld, false, 0, true);

G4Box* alFlange2 = new G4Box("alFlange", 4.637 / 2.*cm, 6.173 / 2.*cm, 1.625 /2. * cm);
G4SubtractionSolid* AlFlange2 = new G4SubtractionSolid("AlFlange2", alFlange2, bucoF1, 0, G4ThreeVector(0.,-0.768*cm,0));
G4LogicalVolume* alFlange2Logical = new G4LogicalVolume(AlFlange2, al, "alFlange2Logical");

G4VPhysicalVolume* alFlange2Physical = new G4PVPlacement(0, G4ThreeVector(0., 0.768*cm, 8.45*cm+2.18*cm+2.07*cm+1.625*cm + 1.625/2.*cm), alFlange2Logical, "alFlange2Physical", fLogicWorld, false, 0, true);


  ////////////////////////////////////////////////////////////////////
  //////////////////////////CHAMBER IN NIOBIUM//////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //Create the foil
G4double havar_radius = 2.185*cm;
G4double havar_thickness = 0.000025*m;

  fSolidFoil = new G4Tubs("Foil", 0, havar_radius, 0.5* havar_thickness, startAngle, spanningAngle);
  fLogicFoil = new G4LogicalVolume(fSolidFoil, Havar, "Foil");
 
  fPhysFoil = new G4PVPlacement(0, G4ThreeVector(0., 0., 15.95*cm + havar_thickness/2.), fLogicFoil, "Foil", fLogicWorld, false, 0, true);
  
  fZ_foil_position =  15.95*cm + havar_thickness/2.;
  fFoilVolume = 3.14*havar_radius*havar_radius*havar_thickness;
  //fPhysFoil
  //  = new G4PVPlacement(0,                                    
	//		G4ThreeVector(0.65*cm,0.*mm, fZ_foil_position),   
	//		fLogicFoil,                           
	//		"Foil",                              
	//		fLogicWorld,                          
	//		false,                               
	//		0,                                    
	//		checkOverlaps);                    

  
  //////////////////////////////////////////////////////////////////
  ////////////////Create the 18O water target ///////////////////////
  //////////////////////////////////////////////////////////////////


  G4double tube_length = 8.45*cm;
  
  //G4double target_innerRadius   = 0.*mm;
  //G4double target_outerRadius   = 2.331*cm;
  
  
  G4double boxX = 0.5*0.95*cm;  
  G4double boxY = 0.5*1.45*cm;
  G4double boxZ = 0.5*1.23*cm;
  
  //G4double boxZ = 0.5*3.43*cm;
  
  G4double target_hz  = 2.*boxZ;
  
  G4Box* box = new G4Box("Box", boxX, boxY, boxZ);
  
  G4double radius = boxX;
  G4double semiZ = boxZ;
  G4Tubs* semicircle = new G4Tubs("Semicircle", 0, radius, semiZ, 0, 180*deg);
  G4ThreeVector posSemiUp(0, boxY, 0);
  G4ThreeVector posSemiDown(0, -boxY, 0);
  G4RotationMatrix* rotSemi = new G4RotationMatrix();
  rotSemi->rotateX(180*deg);
  
  
  G4UnionSolid* waterTarget = new G4UnionSolid("waterTarget", box, semicircle, rotSemi, posSemiDown);
  //G4UnionSolid* heliumTarget = new G4UnionSolid("heliumTarget", Hbox, semicircle, 0, posSemiUp);
  
  G4UnionSolid* completeGeometry = new G4UnionSolid("CompleteGeometry", waterTarget, semicircle, 0, posSemiUp);
  //fSolidTarget = completeGeometry;
  
  //fSolidTarget = new G4Tubs("Tube",
//		 target_innerRadius,
	//	 target_outerRadius,
	//	 target_hz,
	//	 startAngle,
	//	 spanningAngle);
		 
  fTargetVolume = 2.5*cm3;
  
  fLogicTarget
    = new G4LogicalVolume(waterTarget,
			  fTarget_Material,
			  "Target");
			  
  G4LogicalVolume* logicHeliumTarget = new G4LogicalVolume(semicircle, helium, "hTarget");
  

  /*  G4VPhysicalVolume* physTube_PART1 = */
  new G4PVPlacement(0,                                               
		    G4ThreeVector(0.,0.,15.95*cm+ boxZ + havar_thickness),
		    fLogicTarget,                                 
		    "Target",                                    
		    fLogicWorld,                               
		    false,                                           
		    0,                                              
		    checkOverlaps);      
		    
new G4PVPlacement(0,                                               
		    G4ThreeVector(0.,boxY,15.95*cm+ boxZ + havar_thickness),
		    logicHeliumTarget,                                 
		    "heliumTarget",                                    
		    fLogicWorld,                               
		    false,                                           
		    0,                                              
		    checkOverlaps); 
		    
  fTarget_z_position=15.95*cm+ havar_thickness;  
  

              



  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////// PART 5 = NIOBIUM intorno target /////////////////////////////////  
  ///////////////////////////////////////////////////////////////////////////////////////
  
   G4Box* rearF = new G4Box("niobiumBox", 4.637 / 2.*cm, 6.173 / 2.*cm, 1.75 / 2. * cm);
     
   //G4Tubs* buco2 = new G4Tubs("buco2", 0, 0.75*cm, 1.8*cm, startAngle, spanningAngle);
   G4SubtractionSolid* boxHole2 = new G4SubtractionSolid("boxHole2", rearF,completeGeometry , 0, G4ThreeVector(0.,-0.768*cm,-0.875*cm));

   G4LogicalVolume* rearFLogical = new G4LogicalVolume(boxHole2, niobium, "rearFLogical");
	
	G4VisAttributes* rearFVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
	
rearFLogical->SetVisAttributes(rearFVisAtt);
	
   G4VPhysicalVolume* rearFPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0.768*cm, 15.95*cm+1.75/2*cm + havar_thickness), rearFLogical, "rearFPhysical", fLogicWorld, false, 0, true);
	
	//posterior Flange 
	
   G4Box* postF = new G4Box("postFlange", 4.637 / 2.*cm, 6.173 / 2.*cm, 2.8 / 2. * cm);
   G4LogicalVolume* postFLogical = new G4LogicalVolume(postF, al, "postFLogical");
   G4VPhysicalVolume* postFPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0.768*cm, 15.95*cm+1.75*cm+1.4*cm + havar_thickness), postFLogical, "postFPhysical", fLogicWorld, false, 0, true);



  //////////////////////////////////////////////
  /////////// Set sensitive region /////////////
  //////////////////////////////////////////////


  if(!fRegionTarget)
    {
      fRegionTarget = new G4Region("Target");
      fLogicTarget -> SetRegion(fRegionTarget);
     //logicHeliumTarget->SetRegion(fRegionTarget);
      fRegionTarget -> AddRootLogicalVolume(fLogicTarget);
    }

 
  if(!fRegionFoil&&fPhysFoil!=nullptr)
    {
      fRegionFoil = new G4Region("Foil");
      fLogicFoil -> SetRegion(fRegionFoil);
      fRegionFoil -> AddRootLogicalVolume(fLogicFoil);
    }
  
 

  //Always return physical world//

  return physWorld;


}

void STCyclotronDetectorConstruction::ConstructSDandField()
{

  if(fLogicTarget != nullptr){
    STCyclotronSensitiveTarget* TargetSD = new STCyclotronSensitiveTarget("TargetSD", this);
    SetSensitiveDetector("Target",TargetSD);
  }
  
 if(fLogicFoil != nullptr){
    STCyclotronSensitiveFoil* FoilSD = new STCyclotronSensitiveFoil("FoilSD", this);
    SetSensitiveDetector("Foil",FoilSD);
  }

}

void DefineScoring() {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    // Create a new multi-functional detector
    G4MultiFunctionalDetector* detector = new G4MultiFunctionalDetector("ScoringMesh");
    
    // Define a primitive scorer for energy deposition
    G4VPrimitiveScorer* scorer = new G4PSEnergyDeposit("Edep");
    detector->RegisterPrimitive(scorer);

    // Attach the detector to the logical volume of the target
    G4LogicalVolume* fLogicTarget = G4LogicalVolumeStore::GetInstance()->GetVolume("Target");
    sdManager->AddNewDetector(detector);
    fLogicTarget->SetSensitiveDetector(detector);
}






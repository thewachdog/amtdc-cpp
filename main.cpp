/*
C++ Version: C++14
*/

#include <iostream> // basic input output
#include <fstream> // for impoting file
#include <vector> // for vector in c++

// External libraries
#include <json.hpp> // to parse json file
#include <libInterpolate/Interpolate.hpp>
// #include "package/spline.h"

// Read value from JSON file
using json = nlohmann::json;
std::ifstream f("input.json");
json valueDict = json::parse(f);

int main(int argc, char **argv) {

  // Contact wire
  float rho_C = valueDict["contact_linearDensity"];
  float A_C = valueDict["contact_crossSec"];
  A_C *= pow(10, -6);
  float E_C = valueDict["contact_youngs"];
  E_C *= pow(10, 9);
  float I_C = valueDict["contact_inertia"];
  I_C *= pow(10, -9);
  int T_C = valueDict["contact_tension"];
  T_C = -T_C;
  long int EA_C = round(E_C * A_C);
  double EI_C = E_C * I_C;

  // Messenger wire
  float rho_M = valueDict["messenger_linearDensity"];
  float A_M = valueDict["messenger_crossSec"];
  A_M *= pow(10, -6);
  float E_M = valueDict["messenger_youngs"];
  E_M *= pow(10, 9);
  float I_M = valueDict["messenger_inertia"];
  I_M *= pow(10, -9);
  int T_M = valueDict["messenger_tension"];
  T_M = -T_M;
  long int EA_M = round(E_M * A_M);
  double EI_M = E_M * A_M;

  // Droppers
  float rho_D = valueDict["dropper_youngsModulus"];
  double EA_D_temp = valueDict["dropper_youngsModulus"];
  double EA_D = valueDict["dropper_crossSec"];
  EA_D *= EA_D_temp * pow(10, 3);
  float EI_D = valueDict["dropper_inertia"];
  EI_D *= EA_D_temp;

  int No_span = valueDict["miscellaneous_spanCount"]; // Number of Spans

  // Contact wire support
  long int EA_S_CW_temp = valueDict["miscellaneous_youngsModulusRegArm"];
  long int EA_S_CW = valueDict["miscellaneous_crossSecRegArm"];
  EA_S_CW *= pow(10, 9) * EA_S_CW_temp * pow(10, -6);
  long int EI_S_CW = valueDict["miscellaneous_inertiaRegArm"];
  EI_S_CW *= pow(10, 9) * EA_S_CW_temp * pow(10, -6);

  // Messenger wire support
  long int EA_S_MW_temp = valueDict["miscellaneous_youngsModulusSteadyArm"];
  long int EA_S_MW = valueDict["miscellaneous_crossSecSteadyArm"];
  EA_S_MW *= pow(10, 9) * EA_S_MW_temp * pow(10, -6);
  int EI_S_MW = valueDict["miscellaneous_inertiaSteadyArm"];
  EI_S_MW *= pow(10, 9) * EA_S_MW_temp * pow(10, -6);

  int lengthSpan = valueDict["miscellaneous_spanLength"]; // Length of one span
  int T_c = - T_C;
  int T_m = - T_M;
  float w_c = rho_C * 9.8; // # % weight density of the contact wire
  float w_m = rho_M * 9.8; // # % weight density of the messenger wire
  float w_dm = valueDict["miscellaneous_clampMassMessenger"]; // # % weight of dropper clamp on messenger wire
  w_dm *= 9.8;
  float w_dc = valueDict["miscellaneous_clampMassContact"]; // # % weight of dropper clamp on contact wire
  w_dc *= 9.8;

  int No_droppers = valueDict["dropper_count"];

  std::vector<float> dropper_location  =  {};
  for (auto& val : valueDict.items()) {
    // std::cout << val.key() << std::endl;
		if(val.key().substr(0, strlen("dropper_location")) == "dropper_location"){
      dropper_location.push_back(val.value());
    }
  }
  float d_c = valueDict["contact_preSag"]; // # % Maximum presag on the contact wire
  double encumbrance = valueDict["miscellaneous_encumbrance"];

  //   # %% Dropper force and length
  // # % presag calculated at each droppper location
  
  // Need to implement D_dropper
  std::vector<double> D_dropper(No_droppers, 0);

  for(int x = 0; x < No_droppers; x++){
    double a = 4 * d_c / (pow((dropper_location[No_droppers - 1] - dropper_location[0]), 2));
    double b = dropper_location[x] - dropper_location[0];
    double c = dropper_location[No_droppers - 1] - dropper_location[0] - dropper_location[x] + dropper_location[0];
    D_dropper[x] = a * b * c;
    // std::cout << D_dropper[x] << std::endl;
  }

  // # % Forces at the droppers
  std::vector<float> F(No_droppers, 0);
  float f1 = dropper_location[0];
  float f2 = (dropper_location[1] - dropper_location[0])/2;
  float f3 = (T_c * (D_dropper[1] - D_dropper[0])) / (w_c * (dropper_location[1] - dropper_location[0]));
  float F_1_9 = (f1+f2+f3) * w_c + w_dc;
  F[0] = F_1_9;
  F[8] = F_1_9;

  for(int i = 1; i < No_droppers - 1; i++){
    double f4 = (dropper_location[i + 1] - dropper_location[i-1]) / 2;
    double f5 = (T_c * (D_dropper[i] - D_dropper[i-1])) / (w_c * (dropper_location[i] - dropper_location[i-1]));
    double f6 = (T_c * (D_dropper[i+1] - D_dropper[i])) / (w_c * (dropper_location[i+1] - dropper_location[i]));
    F[i] = (f4 - f5 + f6) * w_c + w_dc;
  }

  std::vector<double> r1(No_droppers, 0);
  double r2 = 0.0;
  float R_a = ((w_m * lengthSpan) / 2) + r2;
  
  // # % Reaction force at the messenger wire support
  for(int i = 0; i < No_droppers; i++){
    r1[i] = (F[i] * (lengthSpan - dropper_location[i])) / lengthSpan;
    r2 += r1[i];
  }

  std::vector<double> F_droppers_m(No_droppers, 0);
  std::vector<double> F_droppers_m1(No_droppers, 0);

  for(int i = 0; i < No_droppers; i++){
    for(int j = 0; j < i; j++){
      F_droppers_m1[j] = F[j] * 5.75 * (i - j);
    }
    for(int j = 0; j < No_droppers; j++){
      F_droppers_m[i] += F_droppers_m1[j];
    }
    if(i > No_droppers/2){
      F_droppers_m[i] = F_droppers_m[No_droppers - i - 1];
    }
  }

  std::vector<double> c_m(No_droppers, 0);
  // # # % Sag calculation at the messenger wire at the dropper locations
  for(int i = 0; i < No_droppers; i++){
    float c1 = R_a * dropper_location[i];
    float c2 = pow(w_m * dropper_location[i], 2) / 2;
    float c3 = F_droppers_m[i];
    c_m.push_back((c1 - c2 - c3) / T_m);
    if(i > No_droppers / 2){
      c_m[i] = c_m[No_droppers - i - 1];
    }
  }

  // # % Calculation of dropper lengths 
  std::vector<double> length_dropper(No_droppers, 0);
  for(int i = 0; i < No_droppers; i++){
    length_dropper.push_back(encumbrance - c_m[i] + D_dropper[i]);
    // std::cout << length_dropper[i] << std::endl;
  }

  // # %% Node coordinates

  float EL = 0.25; // # % length of each element in the wire
  int numberNodes_w = round(lengthSpan * 4 + 1); // # % number of elements for contact wire and messenger wire for a single span
  int numberElements = No_span * ((numberNodes_w  - 1) * 2 + No_droppers) + 2 * No_span + 2; // # % total number of elements (220+220+9) contact wire, messenger wire, droppers
  int numberElements_W = (No_span * lengthSpan / EL); // # % Total Number of elements for contact wire and messneger wire
  int numberNodes_W = (No_span * lengthSpan / EL) + 1; // # % Total Number of nodes for contact wire and messneger wire
  
  // std::vector < std::vector<int> > matrix (rows, std::vector<int>(columns, value) );
  std::vector < std::vector<int> > NC_CW (numberNodes_W, std::vector<int>(2, 0) ); // #  % nodal coordinates of the contact wire
  std::vector < std::vector<int> > NC_MW (numberNodes_W, std::vector<int>(2, 0) ); // # % nodal coordinates of the messenger wire
  std::vector<double> mw_x;
  for(double i = 0; i <= lengthSpan; i += EL){
    mw_x.push_back(i); // # % Discretization of the wire into finite elements for sigle span
  }

  int Nodes_effective = (2 * numberNodes_W + 2 * No_span + 2); // # % Total number of nodes of the OHE system

  // # %dropper position and support on the contact wire (x & y coordinates)
  std::vector<double> cW_x(No_droppers+2);
  std::vector<double> cW_y(No_droppers+2);

  // # %dropper position and support on the contact wire (x coordinates)
  for(int i = 1; i <= No_droppers; i++){
    cW_x[i] = dropper_location[i - 1];
    cW_x[No_droppers + 1] = lengthSpan;

    cW_y[i] = D_dropper[i - 1];
  }

  // # % dropper position and support on the contact wire (y coordinates)
  std::vector<double> CW_Y;
  _1D::LinearInterpolator<double> interp1_cy;
  interp1_cy.setData(cW_x, cW_y);
  for(int i = 0; i < mw_x.size(); i++){
    CW_Y.push_back( - interp1_cy(mw_x[i]));
  }

  // # %Interpolation of the dropper point to find out the configuration of the entire contact wire
  std::vector<double> MW_x = cW_x;
  std::vector<double> MW_y(No_droppers + 2);

  for(int i = 0; i < No_droppers + 2; i++){
    if(i == 0 || i == No_droppers + 1){
      MW_y[i] = 1.2;
    }
    else{
      MW_y[i] = length_dropper[i - 1];
    }
  }



}
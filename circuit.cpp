#include "circuit.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <ctime>

using namespace std;

ostream& operator<<(ostream& os, const Fault& fault) {
    os << "Net " << fault.netID << " s-a-" << fault.stuckAt;
    return os;
}

Circuit::Circuit(){
  gateHead = nullptr;
  gateTail = nullptr;
  netHead = nullptr;
  nextGateID = 0;
  originalFaultCount = 0; 
}

Circuit::~Circuit(){
  Gate* currentGate = gateHead;
  while(currentGate != nullptr) {
    Gate* toDelete = currentGate;
    currentGate = currentGate->next;
    delete toDelete;
  }

  Net* currentNet = netHead;
  while(currentNet != nullptr) {
    Net* toDelete = currentNet;
    currentNet = currentNet->next;
    delete toDelete;
  }
}

// === Public Methods ===

int Circuit::getPIcount() const {
  return inputList.size();
}

double Circuit::getCoverage() const {
  if (originalFaultCount == 0) {
    return 0.0;
  } 
  return 100.0 * (double)detectedFaults.size() / originalFaultCount;
}

void Circuit::parseNetList(const string& path){
  ifstream file(path);
  string line;

  if (file.is_open()){
    while (getline(file, line)){
      stringstream ss(line);
      string gate;
      ss >> gate;

      if(gate == "INPUT"){
        int id;
        while(ss >> id && id != -1){
          inputList.push_back(id);
          findNet(id)->sourceGate = -1; 
        }
      } else if (gate == "OUTPUT"){
        int id;
        while(ss >> id && id != -1){
          outputList.push_back(id);
          findNet(id);
        }
      } else {  // Gates
        Gate* newGate = new Gate();
        newGate->gateType = gate;
        newGate->gateID = nextGateID++;

        vector<int> nums;
        int x;
        while(ss >> x){
          nums.push_back(x);
        }

        if(gate == "INV" || gate == "BUF"){ // single input
          newGate->inputNets.push_back(nums[0]);
          newGate->outputNet = nums[1];
        } else { // two inputs
          newGate->inputNets.push_back(nums[0]);
          newGate->inputNets.push_back(nums[1]);
          newGate->outputNet = nums[2];
        }

        for(size_t i = 0; i < newGate->inputNets.size(); i++){
          Net* net = findNet(newGate->inputNets[i]);
          net->connectedGates.push_back(newGate->gateID);
        }
        findNet(newGate->outputNet)->sourceGate = newGate->gateID;
        appendGate(newGate);
      }        
    }
    file.close();
  } else {
    cout << "Unable to open file" << endl;
  }
}

void Circuit::reset(){
  for(Net* n = netHead ; n != nullptr ; n = n->next) {
    n->logicValue = -1; // Reset to unknown
  }
}

vector<int> Circuit::run_fault_free_sim(string primaryInput){
  reset();
  size_t want = inputList.size();
  size_t have = primaryInput.size();
  if (have < want) {
    cout << "[warn] vector shorter than PI count: " << have << " < " << want << "\n";
  }
  for (size_t i = 0; i < want; ++i) {
    int b = (i < have && (primaryInput[i] == '0' || primaryInput[i] == '1')) ? primaryInput[i] - '0' : 0;
    findNet(inputList[i])->logicValue = b;
  }

  bool progress;
  do {
    progress = false;
    for(Gate* g = gateHead ; g != nullptr ; g = g->next) {
      if(g->outputNet != -1 && findNet(g->outputNet)->logicValue == -1) {
        if(inputAssigned(g)) {
          vector<int> inputValues;
          for(int netID : g->inputNets){
            inputValues.push_back(findNet(netID)->logicValue);
          }
          int outputValue = evalLogic(g->gateType, inputValues);
          findNet(g->outputNet)->logicValue = outputValue;
          progress = true;
        }
      }
    }
  } while(progress);

  vector<int> outputValues;
  for(int netID : outputList){
    Net* net = findNet(netID);
    outputValues.push_back(net->logicValue);
  }
  return outputValues;
}

void Circuit::initializeFaults(const string& faultFile){
  allFaults.clear();
  detectedFaults.clear();
  originalFaultCount = 0;

  if(faultFile.empty()){
    cout << "No fault file provided. Generating all faults...\n";
    for (Net* n = netHead; n != nullptr; n = n->next) {
        allFaults.push_back({n->netID, 0}); // s-a-0
        allFaults.push_back({n->netID, 1}); // s-a-1
    }
  } else {
    ifstream r_faultFile(faultFile);

    string line;
    int netID, stuckAt;
    while(getline(r_faultFile, line)){
      stringstream ss(line);
      if(ss >> netID >> stuckAt){
        findNet(netID); // create new net if not exist
        allFaults.push_back({netID, stuckAt});
      }
    }
    r_faultFile.close();
  }
  originalFaultCount = allFaults.size();
}

vector<int> Circuit::run_deductive_sim(string primaryInput){
  vector<int> faultFreeOutput = run_fault_free_sim(primaryInput);

  for (Net* n = netHead; n != nullptr; n = n->next){
    n->faultList.clear();
  }

  // initialize PI faults
  for(size_t i = 0 ; i < inputList.size() ; i++){
    Net* piNet = findNet(inputList[i]);
    int netID = piNet->netID;
    int goodValue = piNet->logicValue;
    Fault selfFault;

    if (goodValue == 0) {
        selfFault = {netID, 1}; // s-a-1
    } else if (goodValue == 1) {
        selfFault = {netID, 0}; // s-a-0
    } else {
        continue;
    }
    
    // Only add fault if it's in our simulation 
    if (check_fault_contain(allFaults, selfFault)) {
        add_fault_2_list(piNet->faultList, selfFault);
    }
  }

  bool progress;
  int iterations = 0;
  int gateCount = 0;
  for(Gate* g = gateHead ; g != nullptr ; g = g->next){
    gateCount++;
  }

  for(Net* n = netHead ; n != nullptr ; n=  n->next){
    if(n->sourceGate != -1){ // means n is NOT PI
      n->logicValue = -1;
    }
  }

  do {
    progress = false;
    for(Gate* g = gateHead ; g != nullptr ;  g = g->next){
      if(g->outputNet != -1 && findNet(g->outputNet)->logicValue == -1){
        if(inputAssigned(g)){
          vector<int> inputValues;

          for(int netID : g->inputNets){
            inputValues.push_back(findNet(netID)->logicValue);
          }

          findNet(g->outputNet)->logicValue = evalLogic(g->gateType, inputValues);
          
          calculate_fault_list(g);

          progress = true;
        }
      }
    }

    if(iterations++ > gateCount + 10){
      cout << "ERROR: Fault simulation time out. \n";
      break;
    }
  }while(progress);

  for(int poNetID : outputList){
    Net* poNet = findNet(poNetID);
    for(const auto& fault : poNet->faultList){
      add_fault_2_list(detectedFaults, fault); 
    }
  }

  vector<Fault> remainingFaults;
  for (const auto& fault : allFaults){
    if(!check_fault_contain(detectedFaults, fault)){
      remainingFaults.push_back(fault);
    }
  }

  allFaults = remainingFaults;  // leave the faults that are not detectable in allFaults

  return faultFreeOutput;
}

void Circuit::reportDetectedFaults(ostream& os) const {
  os << "\n" << string(60, '=') << "\n";
  os << "FAULT SIMULATION REPORT\n";
  os << string(60, '=') << "\n";

  if (originalFaultCount == 0) {
      os << "No faults were simulated.\n";
      return;
  }

  os << "Total Faults Simulated: " << originalFaultCount << "\n";
  os << "Total Faults Detected:  " << detectedFaults.size() << "\n";
  os << "Remaining Faults:       " << allFaults.size() << "\n";
  double coverage = 100.0 * (double)detectedFaults.size() / originalFaultCount;
  os << "Fault Coverage:         " << fixed << setprecision(2) << coverage << "%\n";

  os << "\n--- Detected Faults List ---\n";
  if (detectedFaults.empty()) {
      os << "(None)\n";
  } else {
      for (const auto& fault : detectedFaults) {
          os << "  " << fault << "\n";
      }
  }
}


Circuit::Net* Circuit::findNet(int netID) {
  for(Net* n = netHead ; n != nullptr ; n = n->next) {
    if(n->netID == netID) {
      return n;
    }
  }
  Net* newNet = new Net();
  newNet->netID = netID;
  newNet->next = netHead;
  netHead = newNet;
  return newNet;
}

Circuit::Gate* Circuit::findGate(int gateID) {
  for(Gate* g = gateHead ; g != nullptr ; g = g->next) {
    if(g->gateID == gateID) {
      return g;
    }
  }
  return nullptr;
}

void Circuit::appendGate(Gate* newGate) {
  if(gateHead == nullptr) {
    gateHead = newGate;
    gateTail = newGate;
  } else {
    gateTail->next = newGate;
    gateTail = newGate;
  }
}

int Circuit::evalLogic(const string& gateType, const vector<int>& inputValues) const {
  if(gateType == "AND") {
    return inputValues[0] & inputValues[1];
  } else if(gateType == "OR") {
    return inputValues[0] | inputValues[1];
  } else if(gateType == "NAND") {
    return !(inputValues[0] & inputValues[1]);
  } else if(gateType == "NOR") {
    return !(inputValues[0] | inputValues[1]);
  } else if(gateType == "XOR") {
    return inputValues[0] ^ inputValues[1];
  } else if(gateType == "XNOR") {
    return !(inputValues[0] ^ inputValues[1]);
  } else if(gateType == "INV") {
    return !inputValues[0];
  } else if(gateType == "BUF") {
    return inputValues[0];
  }
  return -1; // Unknown gate type
}

bool Circuit::inputAssigned(Gate* g) {
  for(int netID : g->inputNets){
    Net* net = findNet(netID);
    if(net->logicValue == -1) {
      return false;
    }
  }
  return true;
}


bool Circuit::check_fault_contain(const vector<Fault>& list, const Fault& fault) const {
  for(const auto& f : list ){
    if(f.netID == fault.netID && f.stuckAt == fault.stuckAt){
      return true;
    }
  }
  return false;
}

void Circuit::add_fault_2_list(vector<Fault>& list, const Fault& fault) const { 
  if(!check_fault_contain(list, fault)){ // if already exist in the list then dont add
    list.push_back(fault);
  }
}

// ========== set operation =================

vector<Fault> Circuit::union_op(const vector<Fault>& listA, const vector<Fault>& listB) const {
  vector<Fault> result = listA;
  for(const auto& f : listB){
    add_fault_2_list(result, f);
  }
  return result;
}

vector<Fault> Circuit::intersection_op(const vector<Fault>& listA, const vector<Fault>& listB) const {
  vector<Fault> result;
  for(const auto& f : listA){
    if(check_fault_contain(listB, f)){
      result.push_back(f);
    }
  }
  return result;
}

// L_A - L_B
// *** THIS IS THE FIX ***
vector<Fault> Circuit::difference_op(const vector<Fault>& listA, const vector<Fault>& listB) const {
  vector<Fault> result;
  for(const auto& f : listA){
    if(!check_fault_contain(listB, f)){
      result.push_back(f);
    }
  }
  return result;
} 

void Circuit::calculate_fault_list(Gate* g){
  Net* outNet = findNet(g->outputNet);
  outNet->faultList.clear();

  vector<Fault> propagatedList;

  vector<vector<Fault>> inputFaultList;
  vector<int> inputValues;
  for(int inNetID : g->inputNets){
    Net* n = findNet(inNetID);
    inputFaultList.push_back(n->faultList); 
    inputValues.push_back(n->logicValue);
  }

  const string& gType = g->gateType;  // for simplification, copied to gType

  // ---- evaluate faults to include in the list from input faults ------

  if(gType == "INV" || gType == "BUF"){  // always propagate the input fault
    //L_out = L_in
    if(!inputFaultList.empty()){
      propagatedList = inputFaultList[0];
    }
  }
  else if(gType == "AND" || gType == "NAND"){
    int inputA = inputValues[0];
    int inputB = inputValues[1];

    const auto& L_a = inputFaultList[0];
    const auto& L_b = inputFaultList[1];

    if(inputA == 0 && inputB == 0){ // intersection
      propagatedList = intersection_op(L_a, L_b);
    } else if (inputA == 0 && inputB == 1){ // L_a is controlling value
      propagatedList = difference_op(L_a, L_b);
    } else if (inputA == 1 && inputB == 0){
      propagatedList = difference_op(L_b, L_a);
    } else if (inputA == 1 && inputB == 1){ // union
      propagatedList = union_op(L_a, L_b);
    }
  }
  else if (gType == "OR" || gType == "NOR"){
    int inputA = inputValues[0];
    int inputB = inputValues[1];

    const auto& L_a = inputFaultList[0];
    const auto& L_b = inputFaultList[1];

    if(inputA == 0 && inputB == 0){ // intersection
      propagatedList = union_op(L_a, L_b);
    } else if (inputA == 0 && inputB == 1){ // L_a is controlling value
      propagatedList = difference_op(L_b, L_a);
    } else if (inputA == 1 && inputB == 0){
      propagatedList = difference_op(L_a, L_b);
    } else if (inputA == 1 && inputB == 1){ // union
      propagatedList = intersection_op(L_a, L_b);
    }
  }
  else if (gType == "XOR" || gType == "XNOR"){
    const auto& L_a = inputFaultList[0];
    const auto& L_b = inputFaultList[1];
    
    vector<Fault> u = union_op(L_a, L_b);
    vector<Fault> i = intersection_op(L_a, L_b);
    propagatedList = difference_op(u, i);  // include all the faults EXCEPT faults that changes both input
  }

  outNet->faultList = propagatedList;

  int NetID = outNet->netID;
  int goodValue = outNet->logicValue;
  Fault selfFault;
  if(goodValue == 0){
    selfFault = {NetID, 1};
  } else if (goodValue == 1){
    selfFault = {NetID, 0};
  } else {
    return;
  }

  if(check_fault_contain(allFaults, selfFault)){
    add_fault_2_list(outNet->faultList, selfFault);
  }
}

// ============== new functions for PODEM ===============

int Circuit::logic_NOT(int v) const {
  switch (v) {
    case V_0: return V_1;
    case V_1: return V_0;
    case V_D: return V_D_BAR;
    case V_D_BAR: return V_D;
    case V_X: return V_X;
    default: return V_X;
  }
}

int Circuit::logic_AND(int a, int b) const {
  if (a == V_0 || b == V_0) return V_0; // Controlling value
  if (a == V_1) return b;
  if (b == V_1) return a;
  if (a == V_X || b == V_X) return V_X; // X propagation (Corrected this line)
  if (a == b) return a; // D&D=D, D_BAR&D_BAR=D_BAR
  return V_0; // D & D_BAR = 0
}

int Circuit::logic_OR(int a, int b) const {
    if (a == V_1 || b == V_1) return V_1; // Controlling value
    if (a == V_0) return b;
    if (b == V_0) return a;
    if (a == V_X || b == V_X) return V_X; // X propagation
    if (a == b) return a; // D|D=D, D_BAR|D_BAR=D_BAR
    return V_1; // D | D_BAR = 1
}

int Circuit::logic_NAND(int a, int b) const {
    return logic_NOT(logic_AND(a, b));
}

int Circuit::logic_NOR(int a, int b) const {
    return logic_NOT(logic_OR(a, b));
}

int Circuit::logic_XOR(int a, int b) const {
    if (a == V_X || b == V_X) return V_X;
    if (a == b) return V_0; // D^D = 0
    if (a == V_1) return logic_NOT(b);
    if (b == V_1) return logic_NOT(a);
    if (a == V_0) return b;
    if (b == V_0) return a;
    // Only D and D_BAR left
    return V_1; // D ^ D_BAR = 1
}

int Circuit::evalGate5Val(Gate* g) {
    vector<int> inputs;
    // copy input nets value into vector
    for (int inNetID : g->inputNets) {
        inputs.push_back(findNet(inNetID)->logicValue);
    }

    const string& gType = g->gateType;
    if (gType == "AND") return logic_AND(inputs[0], inputs[1]);
    if (gType == "OR") return logic_OR(inputs[0], inputs[1]);
    if (gType == "NAND") return logic_NAND(inputs[0], inputs[1]);
    if (gType == "NOR") return logic_NOR(inputs[0], inputs[1]);
    if (gType == "XOR") return logic_XOR(inputs[0], inputs[1]);
    if (gType == "XNOR") return logic_NOT(logic_XOR(inputs[0], inputs[1]));
    if (gType == "INV") return logic_NOT(inputs[0]);
    if (gType == "BUF") return inputs[0];

    return V_X;
}

// === PODEM Public Functions ===

std::string Circuit::run_podem(const Fault& fault) {
    targetFault = fault;
    reset5ValLogic(); // Set all nets to V_X

    if (podem_recursive()) {
        return getPIValues(); // SUCCESS
    } else {
        return "UNDETECTABLE"; // FAILURE
    }
}

bool Circuit::isFaultDetected(const Fault& fault) const {
    return check_fault_contain(detectedFaults, fault);
}


// === PODEM Private Helpers (Propagation) ===

void Circuit::reset5ValLogic() {
    for (Net* n = netHead; n != nullptr; n = n->next) {
        n->logicValue = V_X; // Reset all to unknown
    }
    D_frontier.clear(); // Clear the D-frontier
}

// basically converting PI value into string
std::string Circuit::getPIValues() {
    std::string vec = "";
    for (int piNetID : inputList) {
        Net* piNet = findNet(piNetID);
        switch (piNet->logicValue) {
            case V_0: vec += '0'; break;
            case V_1: vec += '1'; break;
            case V_D: vec += '1'; break;
            case V_D_BAR: vec += '0'; break;
            default:  vec += 'X'; break;
        }
    }
    return vec;
}

void Circuit::imply_and_propagate(int piNetID, int value) {
    int finalVal = value;

    // Check if this PI is the specific fault location
    if (piNetID == targetFault.netID) {
        // If we try to set 1, but it's Stuck-at-0 -> It becomes D (1/0)
        if (value == V_1 && targetFault.stuckAt == 0) {
            finalVal = V_D;
        } 
        // If we try to set 0, but it's Stuck-at-1 -> It becomes D_BAR (0/1)
        else if (value == V_0 && targetFault.stuckAt == 1) {
            finalVal = V_D_BAR;
        }
        // If value matches stuck-at, it stays simple 0 or 1 (Fault not excited)
    }

    // Set the PI value
    findNet(piNetID)->logicValue = finalVal;
    
    // Propagate the consequences
    propagate_5val_logic();
}

void Circuit::propagate_5val_logic() {
    bool progress;

    do {
        progress = false;
        for (Gate* g = gateHead; g != nullptr; g = g->next) {
          Net* outNet = findNet(g->outputNet);
          int oldVal = outNet->logicValue;
          
          // Calculate the new output value using 5-val logic
          int newVal = evalGate5Val(g); 

          if(outNet->netID == targetFault.netID){
            int goodVal = newVal; // newVal = where the fault net's value should have been

            if(goodVal == V_0 || goodVal == V_1){ // converts value to D or D-bar to propagate
              if(goodVal == V_1 && targetFault.stuckAt == 0){
                newVal = V_D;
              } else if (goodVal == V_0 && targetFault.stuckAt == 1){
                newVal = V_D_BAR;
              }
            }
          }
          // when any of the net value changes due to ckt simulation, keep going
          if(newVal != oldVal){
            outNet->logicValue = newVal;
            progress = true;
          }
        }

    } while (progress);

    D_frontier.clear();
    for (Gate* g = gateHead; g != nullptr; g = g->next) {
        bool hasD = false;
        for(int inNetID : g->inputNets) {
            int val = findNet(inNetID)->logicValue;
            if (val == V_D || val == V_D_BAR) { // checking D or D_bar is still being propagated and waiting at the input
                hasD = true;
                break; // once D is found, then stop looking
            }
        }

        int outVal = findNet(g->outputNet)->logicValue;
        if (hasD && outVal == V_X) {
            D_frontier.push_back(g->gateID); // if input of the gate is D or Dbar and output of the gate is 'X' then added to D-F
        }
    }
}


// --- PODEM Private Helpers  ---

int Circuit::getControllingValue(const std::string& gateType) const {
    if (gateType == "AND" || gateType == "NAND") return V_0;
    if (gateType == "OR" || gateType == "NOR") return V_1;
    if (gateType == "XOR" || gateType == "XNOR") return V_0; // for XOR choose any value
    return V_X; // BUF, INV, XOR have no controlling value
}

int Circuit::getInversion(const std::string& gateType) const {
    // 1 (true) if inverting, 0 (false) if not
    if (gateType == "INV" || gateType == "NAND" || gateType == "NOR") return V_1;
    return V_0; // No inversion
}

bool Circuit::checkFaultAtPO() { // if true, then there is fault at PO
  for(int poNetID : outputList){
    Net* poNet = findNet(poNetID);
    int val = poNet->logicValue;
    if(val == V_D || val == V_D_BAR){
      return true;  // if any of the PO is either D or D_bar then no error at PO
    }
  }
  return false;
}

bool Circuit::testNotPossible() { // if true, then test not possible
  Net* faultNet = findNet(targetFault.netID);
  int faultVal = faultNet->logicValue;

  if(faultVal != V_D && faultVal != V_D_BAR){
    if(faultVal == V_0 && targetFault.stuckAt == 0) return true;
    if(faultVal == V_1 && targetFault.stuckAt == 1) return true;
  }

  if((faultVal == V_D || faultVal == V_D_BAR) && D_frontier.empty() && !checkFaultAtPO()){
    return true;
  }
  return false;
}

Circuit::PodemTarget Circuit::objective(){
  Net* faultNet = findNet(targetFault.netID);
  int faultVal = faultNet->logicValue;

  if(faultVal == V_X){ // if the fault has not been excited yet, set it as the opposite of s-a value
    int requiredValue = (targetFault.stuckAt == 0) ? V_1 : V_0;
    return {targetFault.netID, requiredValue};
  }

  if(faultVal == V_D || faultVal == V_D_BAR){
    if(D_frontier.empty()){
      return {-1, -1};
    }
    int gateG_ID = D_frontier.front();
    Gate* G = findGate(gateG_ID);
    if(G == nullptr){
      return {-1, -1};
    }

    for(int inputNetID : G->inputNets){ // go over two input Nets
      if(findNet(inputNetID)->logicValue == V_X){ // find the input that has not been assigned yet
        int ctrlVal = getControllingValue(G->gateType); // find the controlling value of the gate
        if(ctrlVal == V_X) {
            ctrlVal = V_0; 
        }
        int nonCtrlVal = logic_NOT(ctrlVal);  // give the non-controlling value to the unassigned input
        return {inputNetID, nonCtrlVal};      // so that D or D_bar can propagate
      }
    }
  }
  return {-1, -1};

}

Circuit::PodemTarget Circuit::backtrace(int k_netID, int vk_value){
  int k = k_netID;
  int v = vk_value;

  while(true){
    Net* k_net = findNet(k);
    if(k_net->sourceGate == -1){ // k is a PI
      break; // we exit from while loop only when k went all the way back to PI
    }

    Gate* G = findGate(k_net->sourceGate); // find the gate that drives this net
    if(G == nullptr){
      break;
    }
    int i = getInversion(G->gateType);

    int j = -1;
    for(int inputNetID : G->inputNets){
      if(findNet(inputNetID)->logicValue == V_X){
        j = inputNetID; // if i find not assigned input value then break
        break;
      }
    }

    if(j == -1){
      j = G->inputNets[0]; // if no input is X, pick the first one
    }

    if(i == V_1){
      v = logic_NOT(v); // if the source gate is inversion gate, then we want v_bar for the input
    }

    k = j; // update the goal, now net(j) is what we have to acheive
  }

  return {k, v}; // finally we return what value v, PI (net k) needs to be 
}

bool Circuit::podem_recursive() {
  // main PODEM loop

  // 1. check for fault not able to propagate to PO
  if(checkFaultAtPO()) return true;

  // 2. check if fault cannot be excited
  if(testNotPossible()) return false;

  // 3. get objective
  PodemTarget objective_pair = objective();
  int k_netID = objective_pair.NetID;
  int vk_value = objective_pair.value;

  if(k_netID == -1){
    return false; // if objective fails, this path is failure
  }

  // 4. backtrace
  PodemTarget pi_assignment = backtrace(k_netID, vk_value);
  int j_piNetID = pi_assignment.NetID;
  int vj_piValue = pi_assignment.value;

  std::vector<int> netValuesBackup;
  for(Net* n = netHead; n != nullptr; n = n->next){
  netValuesBackup.push_back(n->logicValue);
  }

  std::vector<int> dFrontierBackup = D_frontier;

  // 5. imply(j, vj) and recurse
  imply_and_propagate(j_piNetID, vj_piValue);
  if(podem_recursive()){
    return true;
  }


  int i = 0;
  for (Net* n = netHead; n != nullptr; n = n->next, ++i) {
    n->logicValue = netValuesBackup[i];
  }
  D_frontier = dFrontierBackup;

  // 6. imply(j, !vj) reverse decision
  imply_and_propagate(j_piNetID, logic_NOT(vj_piValue));
  if(podem_recursive()){
    return true;
  }

  int j = 0;
  for (Net* n = netHead ; n != nullptr ; n = n->next, ++j){
    n->logicValue = netValuesBackup[j];
  }
  D_frontier = dFrontierBackup;

  return false;

}
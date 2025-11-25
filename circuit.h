#ifndef CIRCUIT_H_
#define CIRCUIT_H_

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>

// --- Logic Values for 5-Valued Simulation ---
// -1 = X (Unknown)
//  0 = 0
//  1 = 1
#define V_X -1
#define V_0 0
#define V_1 1
#define V_D 2
#define V_D_BAR 3

// Forward-declare the ostream operator for Fault
struct Fault;
std::ostream& operator<<(std::ostream& os, const Fault& fault);

struct Fault{
  int netID;
  int stuckAt; // 0 or 1
};

class Circuit
{
  struct Gate {
    std::string gateType;
    int gateID;
    std::vector <int> inputNets;
    int outputNet;
    Gate* next;
    Gate() : gateID (-1), outputNet(-1), next(nullptr) {}
  };

  struct Net {
    int netID;
    std::vector <int> connectedGates;
    int sourceGate;
    int logicValue; 
    std::vector<Fault> faultList;
    Net* next;
    Net() : netID(-1), sourceGate(-1), logicValue(V_X), next(nullptr) {}
  };

  struct PodemTarget {
    int NetID;
    int value;
  };

  private:
    Gate* gateHead;
    Gate* gateTail;
    Net* netHead;

    std::vector<int> inputList;
    std::vector<int> outputList;
    int nextGateID;

    std::vector<Fault> allFaults;
    std::vector<Fault> detectedFaults;

    int originalFaultCount;

    // --- Deductive_FS Helper functions ---
    Net* findNet(int netID);
    Gate* findGate(int gateID);
    void appendGate(Gate* newGate);
    int evalLogic(const std::string& gateType, const std::vector<int>& inputValues) const;
    bool inputAssigned(Gate* g);
    bool check_fault_contain(const std::vector<Fault>& list, const Fault& fault) const;
    void add_fault_2_list(std::vector<Fault>& list, const Fault& fault) const;
    std::vector<Fault> union_op(const std::vector<Fault>& listA, const std::vector<Fault>& listB) const;
    std::vector<Fault> intersection_op(const std::vector<Fault>& listA, const std::vector<Fault>& listB) const;
    std::vector<Fault> difference_op(const std::vector<Fault>& listA, const std::vector<Fault>& listB) const;
    void calculate_fault_list(Gate* g);

    // ------- PODEM ----------
    Fault targetFault;
    std::vector<int> D_frontier; // stores gateID

    void reset5ValLogic();
    std::string getPIValues();
    bool podem_recursive();
    
    PodemTarget objective();
    PodemTarget backtrace(int k, int vk);
    
    void imply_and_propagate(int piNetID, int value);
    void propagate_5val_logic();
    bool checkFaultAtPO();
    bool testNotPossible();
    int getControllingValue(const std::string& gateType) const;
    int getInversion(const std::string& gateType) const;
    int logic_NOT(int v) const;
    int logic_AND(int a, int b) const;
    int logic_OR(int a, int b) const;
    int logic_NAND(int a, int b) const;
    int logic_NOR(int a, int b) const;
    int logic_XOR(int a, int b) const;
    int evalGate5Val(Gate* g);

  public:
    Circuit();  
    ~Circuit(); 

    int getPIcount() const;
    double getCoverage() const;
    void parseNetList(const std::string& path);
    void reset();
    std::vector<int> run_fault_free_sim(std::string primaryInput);
    void initializeFaults(const std::string& faultFile);
    std::vector<int> run_deductive_sim(std::string primaryInput);
    void reportDetectedFaults(std::ostream& os) const;

    // ----- PODEM -----
    std::string run_podem(const Fault& fault);
    bool isFaultDetected(const Fault& fault) const;
};

#endif // CIRCUIT_H_
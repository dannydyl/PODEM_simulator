#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <iomanip> 
#include <map>     
#include <algorithm>

#include "circuit.h"

using namespace std;

bool parseFaultInput(const string& input, Fault& outFault) {
    stringstream ss(input);
    int netID, stuckAt;
    
    if (ss >> netID && ss >> stuckAt) {
        if (stuckAt == 0 || stuckAt == 1) {
            outFault = {netID, stuckAt};
            return true;
        }
    }
    return false;
}


int main(int argc, char** argv) {
    if (argc < 2) {
        cout << "usage: " << argv[0] << " <netlist.txt>\n";
        cout << "This program runs PODEM ATPG interactively." << endl;
        return 0;
    }

    string netlistFile = argv[1];
    
    // --- Load the circuit ---
    Circuit ckt;
    try {
        ckt.parseNetList(netlistFile);
    } catch (...) {
        cout << "Error: Could not parse netlist file: " << netlistFile << endl;
        return 1;
    }

    cout << "Circuit loaded successfully." << endl;
    cout << "Enter a fault to test (e.g., '16 0' for Net 16 s-a-0)." << endl;
    cout << "Type 'exit' or 'quit' to end." << endl;

    string line;
    while (true) {
        cout << "\n> Enter fault (netID saVal): ";
        if (!getline(cin, line)) {
            break; // End of input
        }

        if (line == "exit" || line == "quit") {
            break;
        }

        Fault targetFault;
        if (!parseFaultInput(line, targetFault)) {
            cout << "  Error: Invalid format. Please use: <netID> <0 or 1>" << endl;
            continue;
        }

        cout << "  Targeting fault: " << targetFault << endl;

        // Run PODEM 
        Circuit podemCkt;
        podemCkt.parseNetList(netlistFile);

        string testVector = podemCkt.run_podem(targetFault);

        cout << "  PODEM result:    " << testVector << endl;

        // Verify the test vector with deductive sim 
        if (testVector != "UNDETECTABLE") { // if the test is not undetectable
            string simVector = testVector;
            // Replace all 'X' with '0'
            std::replace(simVector.begin(), simVector.end(), 'X', '0'); 
            cout << "  Verifying with:  " << simVector << endl;

            Circuit simCkt;
            simCkt.parseNetList(netlistFile);

            // Initialize the fault simulator with this one fault
            string tempFaultFile = "temp_single_fault.txt";
            ofstream faultOut(tempFaultFile);
            faultOut << targetFault.netID << " " << targetFault.stuckAt << "\n";
            faultOut.close();
            
            simCkt.initializeFaults(tempFaultFile);

            // Run the deductive simulator
            simCkt.run_deductive_sim(simVector);

            // Check if the fault was detected
            if (simCkt.isFaultDetected(targetFault)) {
                cout << "  Verification:    SUCCESS" << endl;
            } else {
                cout << "  Verification:    *** FAILED ***" << endl;
            }

            remove(tempFaultFile.c_str());
            
        } else {
            cout << "  Verification:    N/A (fault is undetectable)" << endl;
        }
    }

    cout << "Exiting." << endl;
    return 0;
}

// how to run
// ./main file/s27.txt
// user input format
// netID saVal
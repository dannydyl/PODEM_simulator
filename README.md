# PODEM Test Vector Generation and Deductive Fault Simulation

## Overview

This project provides a C++ implementation of two essential tools for digital circuit testing.

1. **PODEM Path Oriented Decision Making**  
   An Automatic Test Pattern Generation method that finds an input vector able to detect a specific stuck at fault.

2. **Deductive Fault Simulator**  
   A method to evaluate fault coverage for a sequence of input vectors across an entire fault list.

Both programs support standard benchmark circuits such as ISCAS eighty five netlists.

## Key Features

* Automatic Test Pattern Generation using the PODEM method
* Efficient fault coverage evaluation using deductive fault simulation
* Flexible simulation modes for fault lists and input vectors
* C plus plus implementation using STL containers for circuit modeling

## Implementation Details

The circuit is represented with clear data structures that allow efficient simulation.

### Data Structures

| Component | Description | Implementation |
| :--- | :--- | :--- |
| Circuit topology | Gate and net connectivity | Linked lists |
| Core sets | D Frontier, Primary Input list, Primary Output list | std::vector |
| Logic states | X, zero, one, D, D bar | Macro constants |

### PODEM Method Summary

The main logic is implemented inside the function podem_recursive.  
It works similarly to a depth first search.

Steps:

1. **Check observability**  
   Determine whether the fault effect D or D bar reaches a Primary Output. A successful detection ends the search.

2. **Check excitability**  
   Confirm that the target stuck at fault is activated.

3. **Determine objective**  
   Select the required target net and logic value to either activate the fault or propagate the D Frontier.

4. **Backtrace**  
   Find a Primary Input assignment that satisfies the objective.

5. **Simulate and backtrack**  
   Apply the assignment, simulate implications, and restore previous states when a conflict occurs.

## Getting Started

### Requirements

* C plus plus compiler with C plus plus twenty support  
* Make utility

### Build

Run:

make

Clean generated files:

make clean

## Running the Program

The application provides two main functions.

* PODEM for interactive ATPG
* Deductive fault simulation in three modes

---

## PODEM Simulation

Interactive ATPG mode loads a circuit and allows repeated fault entry.

### Command

./main netlist_file.txt

### Example

$ ./main netlist/s27.txt
Circuit loaded successfully.
Enter a fault to test (for example 16 0)
Type exit or quit to end.

> Enter fault (netID saVal): 16 0
Targeting fault: Net 16 s a 0
PODEM result: X0X10X0
Verifying with: 0001000
Verification: SUCCESS

---

## Deductive Fault Simulation

Three operation modes are supported.

### Mode One  
File based simulation using an input vector file.

Command:

./main vector_file optional_fault_list

Output:  
result.txt

### Mode Two  
Random vector simulation with progressive coverage logging.

Command:

./main random

Output:  
coverage.csv

### Mode Three  
Interactive simulation where the user enters faults and input vectors.

Command:

./main interactive

---

## Test Results Summary

The implementation was verified on several ISCAS circuits.

| Circuit | Target Fault | PODEM Result | Verification |
| :--- | :--- | :--- | :--- |
| s27 | Net sixteen s a zero | X0X10X0 | SUCCESS |
| s27 | Net ten s a one | X00XXX0 | SUCCESS |
| s298f two | Net seventy s a one | 01X1XXXXXXXXXX0XX | SUCCESS |
| s344 two | Net one hundred sixty six s a zero | 01X00XXXXX011XX0XXXX XXXXX | SUCCESS |

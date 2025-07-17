#include "slbe.hpp"

int main(int argc, char* argv[])
{
  bigend = fBigEndian();
  fDefineSystem();
  fSetSerialDomain();
  fStartDLMESO();
  fMemoryAllocation();
  fGetModel();
  fPrintSystemInfo();
  fsPrintDomainInfo();
  fInputParameters();
  fPrintParameters();
  fReadSpaceParameter();
  fsCreateIOGroups();
  fNeighbourBoundary();
  fInitializeSystem();               // assuming no restart for this simulation
  fReadInitialState();
  if(outformat==2)                   // Grid for Plot3D output files
    fsOutputGrid(); 
//  if(interact==1 || interact==2)
//    fCalcPotential_ShanChen();
  timetotal=fCheckTimeSerial();
  for(lbcurstep=0; lbcurstep<= lbtotstep; lbcurstep++) {
    if(lbcurstep%lbsave == 0 || lbcurstep==lbequstep) {
      if(lbequstep>0 && lbcurstep==lbequstep) {
        fPrintEndEquilibration();
        postequil = 1;
        fReadSpaceParameter();
        fNeighbourBoundary();
      }
      if(lbcurstep>=lbequstep)
        fOutput();
      cout << left << setw(12) << lbcurstep << " ";
      fPrintDomainMass();
      double timeelapsed = fCheckTimeSerial();
      cout << left << setw(12) << timeelapsed << " ";
      fPrintDomainMomentum();
    }
    fInteractionForceZero();         // Remove if no Boussinesq forces or mesophase interactions used
    fsCalcPhaseIndex_Lishchuk();     // Substitute with required potential/phase index/density/concentration gradients
                                     // (or remove if no mesophase interactions required)
    fsInteractionForceLishchuk();    // Substitute with required mesophase interactions and/or add additional forces
    fCollisionBGKGuoLishchuk();      // Substitute with required collision and forcing algorithm
    fPostCollBoundary();
    fPropagationSwap();
    fPostPropBoundary();
  }
  timetotal=fCheckTimeSerial();
  fFreeMemory();
  fFinishDLMESO();
  return 0;
}

    

//compute_LOS_profiles.cpp -- compute LOS tau profiles for a specified
//spacecraft geometry

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "definitions.h"
#include "nr3.h"
#include "iph_sim.h"

int main() {
  using std::cout; 
  using std::endl;
  
  LOSprofinterp LOS_prof(losproffname);
  
  //coronal parameters (number density and temperature)
  std::cout << "Enter nexo: ";
  double nexo;
  std::cin >> nexo;
  std::cout << "nexo = " << nexo << endl;
  std::cout << "Enter Texo: ";
  double Texo;
  std::cin >> Texo;
  std::cout << "Texo = " << Texo << endl;
  //get the most-probable thermal velocity
  double mpvel = sqrt(2*kB*Texo/(mH));//cm per sec
  
  //new file-based computation
  std::cout << "Input filename: ";
  string infilename;
  std::cin >> infilename;
  std::cout << "Reading geometry values from file: " << infilename << endl;
  ifstream infile;
  infile.open(infilename.c_str());
  
  //prepare for read in
  int nprof=0; //number of profiles to compute
  infile >> nprof;
  //setup vectors for file reading
  VecDoub sc_alt_km;
  sc_alt_km.resize(nprof);
  VecDoub lc;
  lc.resize(nprof);
  VecDoub tanpt_alt_km;
  tanpt_alt_km.resize(nprof);
  VecDoub lb;
  lb.resize(nprof);
  VecDoub sza_deg;
  sza_deg.resize(nprof);
  VecDoub t0;
  t0.resize(nprof);
  //and for output
  VecDoub *avec;
  avec = new VecDoub[nprof];
  VecDoub *tvec;
  tvec = new VecDoub[nprof];
  for (int iprof=0; iprof<nprof; iprof++) {
    std::cout << "For point " << iprof << endl;

    infile >> sc_alt_km[iprof];
    //the file has these in km altitude, need to convert to lambda= GMm/kTr
    lc[iprof] = (G*mMars*mH)/(kB*Texo*(rMars+200e5));
    // std::cout << "   lc = " << lc[iprof] << endl;

    infile >> tanpt_alt_km[iprof];
    lb[iprof] = (G*mMars*mH)/(kB*Texo*(rMars+tanpt_alt_km[iprof]*1e5));
    // std::cout << "   lb = " << lb[iprof] << endl;

    //add some error trapping logic in case alt_lb > alt_lc, in which
    //case the code would fail:
    if (lb[iprof] > lc[iprof]) {
      std::cout << "Tangent point below critical level in index "
		<< iprof << "of input file. Errors will result. Aborting!";
      return 1;
    }
    infile >> sza_deg[iprof];
    //the angle here is a Solar Zenith Angle, but as these are Solar
    //occultations the transformation to theta_0 is easy (as long as
    //the tangent point altitude is being computed correctly, and the
    //tangent point is between the spacecraft and the Sun):
    t0[iprof] = pi/180.*sza_deg[iprof];
    t0[iprof] = pi/2-t0[iprof];//0->pi, pi/2->0, pi->-pi/2
    // std::cout << "   t0 = " << t0[iprof] << endl;

    // std::cin.get();

    //compute the distributions
    LOS_prof.interp(lc[iprof],lb[iprof],t0[iprof],avec[iprof],tvec[iprof]);

    //multiply by most probable velocity to convert to cm/s and
    //         by user specified density to obtain optical depth
    for (int itau=0;itau<tvec[iprof].size();itau++) {
      avec[iprof][itau] *= mpvel;
      avec[iprof][itau] /= 1e5;//convert from cm/s to km/s
      tvec[iprof][itau] *= nexo;
      //print these values to console
      // std::cout << endl;
      // std::cout << "   avec["<< iprof << "][" << itau << "] = " << avec[iprof][itau] << endl;
      // std::cout << "   tvec["<< iprof << "][" << itau << "] = " << tvec[iprof][itau] << endl;
      // std::cout << endl;
    }
    //    std::cin.get();
  }
  infile.close();
  
  //prepare for write out
  std::cout << "Output filename: ";
  string outfilename;
  std::cin >> outfilename;
  std::cout << "Writing profile values to file: " << outfilename << endl;
  ofstream outfile;
  outfile.open(outfilename.c_str());
  //find the largest vector of taus
  int maxsize = 0;
  int maxpos = 0;
  for (int iprof=0; iprof<nprof; iprof++) {
    maxsize = tvec[iprof].size() > maxsize ? tvec[iprof].size() : maxsize;
    maxpos = tvec[iprof].size() > maxsize ? iprof : maxpos;
  }
  //now write out to file
  outfile << "Each line of this file gives optical depths as a function of Doppler shift.\n";
  outfile << "The first row gives the names of the parameters and velocities in km/s at \n";
  outfile << "which the optical depths are measured, in the rest frame of the planet.\n";
  outfile.width(10);
  outfile << "\"i\"";
  outfile << ",";
  outfile.width(10);
  outfile << "\"sc_alt\"";
  outfile << ",";
  outfile.width(10);
  outfile << "\"tanpt_alt\"";
  outfile << ",";
  outfile.width(10);
  outfile << "\"sza_deg\"";
  outfile << ",";
  for (int ivel=0; ivel < maxsize-1; ivel++) {
    outfile.width(10);
    outfile << avec[maxpos][ivel];
    outfile << ",";
  }
  outfile.width(10);
  outfile << avec[maxpos][maxsize-1];
  outfile << "\n";
  int thissize=0;
  int padzeroes=0;
  for (int iprof=0; iprof<nprof; iprof++) {
    //how big is this array?
    thissize = tvec[iprof].size();
    padzeroes= (maxsize-thissize)/2;
    
    //print the identifiers
    outfile.width(10);
    outfile << iprof;
    outfile << ",";
    outfile.width(10);
    outfile << sc_alt_km[iprof];
    outfile << ",";
    outfile.width(10);
    outfile << tanpt_alt_km[iprof];
    outfile << ",";
    outfile.width(10);
    outfile << sza_deg[iprof];
    outfile << ",";
    
    //now print the taus, with the zero padding for arrays smaller than the max size
    for (int ipad=0;ipad<padzeroes;ipad++) {
      outfile.width(10);
      outfile << 0.0;
      outfile << ",";
    }
    for (int ivel=0; ivel < thissize-1; ivel++) {
      outfile.width(10);
      outfile << tvec[iprof][ivel];
      outfile << ",";
    }
    outfile.width(10);
    outfile << tvec[iprof][thissize-1];
    if (padzeroes>0) {
      outfile << ",";
    }
    for (int ipad=0;ipad<padzeroes-1;ipad++) {
      outfile.width(10);
      outfile << 0.0;
      outfile << ",";
    }
    if (padzeroes > 1) {
      outfile.width(10);
      outfile << 0.0;
    }
    outfile << std::endl;
  }
  outfile.close();
  
  return 0;
}


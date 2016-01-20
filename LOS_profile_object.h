//LOS_profile_object.cpp -- object to compute LOS tau profiles for a specified
//spacecraft geometry, mostly to pass to Python.

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "definitions.h"
#include "nr3.h"
#include "iph_sim.h"

using std::cout; 
using std::endl;

struct H_LOS_simulator {
  //LOS interpolation object
  LOSprofinterp LOS_prof;
  //coronal parameters
  double nexo, Texo;
  double lc, mpvel;
  //spacecraft observation parameters
  int nprof;
  VecDoub sc_alt_km;
  VecDoub tanpt_alt_km;
  VecDoub lb;
  VecDoub sc_sza_deg;
  VecDoub t0;
  bool scposloaded;
  //output parameters (arrays containing LOS profiles for each spacecraft position)
  VecDoub* avec;
  VecDoub* tvec;
  bool simcomplete;

  H_LOS_simulator() {
    scposloaded=FALSE;
    simcomplete=FALSE;
  };

  ~H_LOS_simulator() {
    if (simcomplete) {
      delete [] avec;
      delete [] tvec;
    }
  };
    
  void set_scpos(double* tanpt_alt_kmm, double* sc_sza_degg, int nproff) {
    VecDoub tanpt_vec(nproff,tanpt_alt_kmm);
    VecDoub sc_sza_vec(nproff,sc_sza_degg);
    set_scpos(tanpt_vec,sc_sza_vec);
  }

  void set_scpos(VecDoub tanpt_alt_kmm, VecDoub sc_sza_degg) {
    nprof=tanpt_alt_kmm.size();
    tanpt_alt_km=tanpt_alt_kmm;
    sc_sza_deg=sc_sza_degg;
    lb.resize(tanpt_alt_km.size());
    t0.resize(tanpt_alt_km.size());
    for (int iprof=0;iprof<tanpt_alt_km.size();iprof++) {
      t0[iprof] = pi/180.*sc_sza_deg[iprof];
      t0[iprof] = pi/2-t0[iprof];//0->pi, pi/2->0, pi->-pi/2
    }
    scposloaded=TRUE;
  }

  void read_scpos_from_file(string infilename) {
    std::cout << "Reading geometry values from file: " << infilename << endl;
    ifstream infile;
    infile.open(infilename.c_str());
    
    //prepare for read in
    infile >> nprof;
    //setup vectors for file reading
    sc_alt_km.resize(nprof);
    tanpt_alt_km.resize(nprof);
    lb.resize(nprof);
    sc_sza_deg.resize(nprof);
    t0.resize(nprof);
    for (int iprof=0; iprof<nprof; iprof++) {
      //std::cout << "For point " << iprof << endl;
      
      infile >> sc_alt_km[iprof];

      infile >> tanpt_alt_km[iprof];
      
      infile >> sc_sza_deg[iprof];
      //the angle here is a Solar Zenith Angle, but as these are Solar
      //occultations the transformation to theta_0 is easy (as long as
      //the tangent point altitude is being computed correctly, and the
      //tangent point is between the spacecraft and the Sun):
      t0[iprof] = pi/180.*sc_sza_deg[iprof];
      t0[iprof] = pi/2-t0[iprof];//0->pi, pi/2->0, pi->-pi/2
      // std::cout << "   t0 = " << t0[iprof] << endl;
      
      // std::cin.get();
    }
    scposloaded=TRUE;
  }
  
  void simulate(double nexoo, double Texoo) {
    if (!scposloaded) {
      cout << "Load spacecraft positions before calling simulate()!";
      return;
    }
    //set coronal parameters (number density and temperature)
    nexo=nexoo;
    Texo=Texoo;
    //get the most-probable thermal velocity
    mpvel = sqrt(2*kB*Texo/(mH));//cm per sec
    //and the Jeans parameter
    lc = (G*mMars*mH)/(kB*Texo*(rMars+200e5));

    //set up parameters for output
    avec = new VecDoub[nprof];
    tvec = new VecDoub[nprof];
    for (int iprof=0; iprof<nprof; iprof++) {
      lb[iprof] = (G*mMars*mH)/(kB*Texo*(rMars+tanpt_alt_km[iprof]*1e5));
      //    std::cout << "   lb = " << lb[iprof] << endl;
      
      //add some error trapping logic in case alt_lb > alt_lc, in which
      //case the code would fail:
      if (lb[iprof] > lc) {
	cout << "Tangent point below critical level in index "
	     << iprof << "of requested profile. Errors will result. Aborting!";
	throw("bad alatitudes in read_scpos_from_file");
      }
      

      //compute the distributions
      LOS_prof.interp(lc,lb[iprof],t0[iprof],avec[iprof],tvec[iprof]);
      
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
    pad_profiles();
    simcomplete=TRUE;
  }
  
  void pad_profiles() {
    //pads all avec,tvec profiles to a common length, adding zeroes at the ends

    //find the largest vector of taus
    int maxsize = 0;
    int maxpos = 0;
    for (int iprof=0; iprof<nprof; iprof++) {
      maxsize = tvec[iprof].size() > maxsize ? tvec[iprof].size() : maxsize;
      maxpos = tvec[iprof].size() > maxsize ? iprof : maxpos;
    }
    VecDoub ttvec;

    int thissize=0;
    int padzeroes=0;
    for (int iprof=0; iprof<nprof; iprof++) {
      //how big is this array?
      thissize = tvec[iprof].size();
      padzeroes= (maxsize-thissize)/2;
      ttvec=tvec[iprof];

      tvec[iprof].resize(maxsize);
      for (int ivec=0; ivec<maxsize; ivec++) {
	if (ivec<padzeroes || ivec > maxsize-padzeroes-1) {
	  tvec[iprof][ivec]=0.0;
	} else {
	  tvec[iprof][ivec]=ttvec[ivec-padzeroes];
	}
      }
      avec[iprof]=avec[maxpos];
    }
  }

  void write_to_file(string outfilename) {
    if (!simcomplete) {
      cout << "Run simlate() before calling write_to_file()!";
      return;
    }

    //prepare for write out
    std::cout << "Writing profile values to file: " << outfilename << endl;
    ofstream outfile;
    outfile.open(outfilename.c_str());
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

    int nvel=avec[0].size();
    for (int ivel=0; ivel < nvel-1; ivel++) {
      outfile.width(10);
      outfile << avec[0][ivel];
      outfile << ",";
    }
    outfile.width(10);
    outfile << avec[0][nvel-1];
    outfile << "\n";
    for (int iprof=0; iprof<nprof; iprof++) {
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
      outfile << sc_sza_deg[iprof];
      outfile << ",";
    
      //now print the taus, with the zero padding for arrays smaller than the max size
      for (int ivel=0; ivel < nvel-1; ivel++) {
	outfile.width(10);
	outfile << tvec[iprof][ivel];
	outfile << ",";
      }
      outfile.width(10);
      outfile << tvec[iprof][nvel-1];
      outfile << std::endl;
    }
    outfile.close();
  }
  
};

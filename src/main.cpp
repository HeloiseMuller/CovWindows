/*
 * Copyright 2021, Héloïse Muller
 *
 * heloise.muller@egce.cnrs-gif.fr
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *
 * */

#include <vector>
#include <string>
#include <array>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stdexcept>
#include <cmath> 

#include <unistd.h>
#include <omp.h>

class Window {
  public:
     Window(std::string chrom, uint64_t begin, uint64_t end, std::string wnames="NA"): chrom_(chrom), begin_(begin), end_(end), wnames_(wnames), average_(0.) {}
    ~Window() = default;

    const std::string& chrom() const {return chrom_;}
    uint64_t begin() const {return begin_;}
    uint64_t end() const {return end_;}
    const std::string& wnames() const {return wnames_;}
    
    double average() const {return average_;}
    void set_average(double a) {average_ = a;}

  private:
    std::string chrom_;
    uint64_t begin_;
    uint64_t end_;
    std::string wnames_;
    double average_;
};

class SCF {
  public:
    SCF(): chrom_(), positions() {}
    SCF(std::string chrom): chrom_(chrom), positions() {}
    ~SCF() = default;

    void add_position(std::pair<uint64_t, double> pos) {
      positions.push_back(pos);
    }

    size_t size() const {return positions.size();}

    const std::pair<uint64_t, double>& operator[](size_t i) const {
      if(positions.size() == 0 || i >= positions.size()) {
	std::string mssg = "Invalid index to SCF of " + std::to_string(i);
        throw std::out_of_range(mssg);
      }

      return positions[i];
    }

    std::pair<uint64_t, double>& operator[](size_t i) {
      if(positions.size() == 0 || i >= positions.size()) {
	std::string mssg = "Invalid index to SCF of " + std::to_string(i);
        throw std::out_of_range(mssg);
      }

      return positions[i];
    }

  private:
    std::string chrom_;

    // Vector of all positions (pair.first) and how many times
    // they were read (pair.second)
    std::vector<std::pair<uint64_t, double>> positions;
};

std::vector<Window> read_windows(const std::string& fname, bool windowNames) {
  std::vector<Window> windows;

  std::ifstream file(fname);

  std::string chrom;
  uint64_t begin;
  uint64_t end;
  std::string wnames;
  int i = 0;

  std::string element;
  while(file >> element) { //Go through file word by word
    if(i == 0) {
      chrom = element;
      i++;
    } else if(i == 1) {
      if(std::stoll(element) < 0) {
        std::cerr << "ERROR: begin position negative\n";
        std::exit(1);
      }
      begin = std::stoull(element) -1; //scf file begin at 1 whereas index begin at 0 in c++
      i++;
    } else if(i == 2) {
      if(std::stoull(element) <= begin) {
       std::cerr << "ERROR: Found end equal or inferior to begin\n"; 
       std::exit(1);
      }
      end = std::stoull(element) -1; //scf file begin at 1 whereas index begin at 0 in c++
      if(windowNames==false){
        // Save map
        windows.push_back(Window(chrom, begin, end));
        // Reset i
        i = 0;
      } else {
        i++;
      }
    } else if(i == 3 && windowNames==true){
      wnames = element;
      // Save map
      windows.push_back(Window(chrom, begin, end, wnames));
      // Reset i
      i = 0;
    }  
  }

  if(i != 0 && element != "\n") {
    std::cerr << " ERROR: Unexpected dangling values in windows file\n";
    std::exit(1);
  }

  file.close();

  return windows;
}

std::unordered_map<std::string, SCF> read_scfs(const std::string& fname) {
  std::unordered_map<std::string, SCF> scfs;

  std::ifstream file(fname);

  std::string chrom;
  uint64_t pos;
  double count;
  int i = 0;

  std::string element;
  while(file >> element) {
    if(i == 0) {
      // Ensure same SCF
      if(element != chrom) {
        chrom = element;
        scfs[chrom] = SCF(chrom);
      }
      i++;
    } else if(i == 1) {
      pos = std::stoull(element) - 1;
      i++;
    } else if(i == 2) {
      count = std::stod(element); 

      // Save position
      scfs[chrom].add_position({pos, count});

      // Reset i
      i = 0;
    }
  }

  if(i != 0 && element != "\n") {
    std::cerr << " ERROR: Unexpected dangling values in coverage file\n";
    std::exit(1);
  }

  file.close();

  return scfs;
}

void get_window_average(Window& window, const SCF& scf) {
  double sum = 0.;
  double count = 0.;

  //Begin of window not included so +1 / End included so <=
  for(size_t i = window.begin()+1; i <= window.end(); i++) {
    sum += scf[i].second;
    count += 1.;
  }

  window.set_average(sum / count);
}

double get_whole_average(std::vector<Window> windows){
  double sum = 0;
  double weight = 0;
    for(size_t i = 0; i < windows.size(); i++) {
      weight += windows[i].end()-windows[i].begin();
      sum += (windows[i].end()-windows[i].begin()) * windows[i].average();
    }  
  return sum/weight;
}


const char* OPTIONS = "fdht:w:c:o:";

void usage() {
  std::cout << " Avebed -w <window file> -c <coverage file> [-o <output file>] [-t <num threads>] [-f] [-d] [-h]\n\n";

  std::cout << "   -w <window file>   : <window file> is name of file containing windows in a bed format: ]n:m]\n";
  std::cout << "   -f                 : To set when <window file> has a fourth colomn containing their names\n";
  std::cout << "   -c <coverage file> : <coverage file> is name of file containing coverage at each positions\n";
  std::cout << "   -d                 : Give the average coverage by window, instead of the whole average\n"; 
  std::cout << "   -o <output file>   : <output file> is the name of the file to save results in, when -d is set\n";
  std::cout << "   -t <num threads>   : If not used, default number is 1\n";
  std::cout << "   -h                 : Show this message and exit\n\n";

  std::cout << " Writen by Héloïse Muller (heloise.muller@egce.cnrs-gif.fr)\n";
  std::cout << " Copyright (c) 2021, Héloïse Muller\n";
  std::cout << " All rights reserved.\n\n";

  std::cout << " Released under the terms and conditions of the CeCILL-v2.1 license.\n";
}

int main(int argc, char** argv) {

  // Parameters
  std::string windows_file;
  std::string coverage_file;
  std::string out_file;
  bool details = false;
  int nthreads = 1;
  bool windowNames = false;

  // Parse command line arguments
  int opt = 0;
  while((opt = getopt(argc, argv, OPTIONS)) != -1) {
    switch (opt) {
      case 't':
        nthreads = std::stoi(optarg);
        break;

      case 'w':
        windows_file = std::string(optarg);
        break;

      case 'f':
        windowNames = true;
        break;

      case 'c':
        coverage_file = std::string(optarg);
        break;

      case 'd':
        details = true;
        break;

      case 'o':
        out_file = std::string(optarg);
        break;

      case 'h':
        usage();
        return 0;
        break;
      
      default:
        usage();
        return 1;
        break;
    }
  }

  if(windows_file.size() == 0 || coverage_file.size() == 0 || (out_file.size() == 0 && details==true)) {
    std::cerr << " ERROR: Improper arguments provided.\n\n";
    usage();
    return 1;
  }

  // Maps to hold windows and SCFs
  std::cout << " Reading windows file...\n";
  std::vector<Window> windows = read_windows(windows_file, windowNames);

  std::cout << " Reading coverage file...\n";
  std::unordered_map<std::string, SCF> scfs = read_scfs(coverage_file);

  // Set number of threads
  int nmax_threads = omp_get_max_threads();
  if(nmax_threads < nthreads) nthreads = nmax_threads;
  omp_set_num_threads(nthreads);

  // Go through each window, and get the average for the window
  std::cout << " Calculating window averages...\n";
  #pragma omp parallel for schedule(static)
  for(size_t i = 0; i < windows.size(); i++) {
    if (scfs.find(windows[i].chrom()) == scfs.end()) {
      #pragma omp critical
      {
        std::cerr << "Chrom " << windows[i].chrom() << " not present in coverage file.\n";
        //If windows from BUSCO, check that no contig of named "contig_1:4-1933"
      }
      std::exit(1);
    } else {
      get_window_average(windows[i], scfs.at(windows[i].chrom()));
    }
    
  }

  //Get the whole average
  double average = get_whole_average(windows);

//Calculate the whole standart deviation !! TO FINISH !!
std::cout << " Calculating the standart deviation...\n";
double sum = 0;
double count = 0;
  for(size_t i = 0; i < windows.size(); i++) { //Go through windows
    //+1 to reset positions ; another +1 because don't count 1st position:
    for(size_t j = windows[i].begin()+1; j <= windows[i].end(); j++) { 
      sum += pow(scfs.at(windows[i].chrom())[j].second-average, 2);
      //std::cout << scfs.at(windows[i].chrom())[j].second << "\n";
      count += 1.;    
    }  
  }
  //std::cout << count << "\n";
  //std::cout << sum << "\n\n";
  double stdev = sqrt(sum/count);

  //Print the whole average with standart deviation
  std::cout << "The whole average through these windows is: " << average << " +- " << stdev << "\n";

  // Write output file
  if(details){
    std::cout << " Writing detailed output file...\n";
    std::ofstream file(out_file);
    for(size_t i = 0; i < windows.size(); i++) {
      file << windows[i].chrom() << "  ";
      file << windows[i].begin()+1 << "  "; //+1 to reset like in input
      file << windows[i].end()+1 << "  "; //+1 to reset like in input
      if(windowNames==true){
        file << windows[i].wnames() << " ";
      }
      file << windows[i].average() << "\n";
    }
    file.close();
   }
  return 0;
}

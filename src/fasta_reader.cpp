#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <Rcpp.h>
#include <ctime>

// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace Rcpp;

std::ifstream::pos_type filesize(const char* filename)
{
  std::ifstream input(filename, std::ifstream::ate | std::ifstream::binary);
  return input.tellg();
}
// [[Rcpp::export]]
DataFrame read_fasta(const vector<string> & files, const vector<string> & proteins)
{
  set<string> set;
  for(auto protein : proteins)
  {
    set.insert(protein);
  }
  time_t interval = time(nullptr);
  vector<string> accessions;
  vector<string> descriptors;
  vector<string> sequences;
  int fsize = 0;
  int bytes = 0;
  for (auto f : files)
  {
    fsize += filesize(f.c_str());
  }
  for (auto f : files)
  {
    ifstream file(f);
    string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      if(set.find(entry.substr(0, entry.find_first_of(' '))) != set.end())
      {
        while (entry.find('\n') == std::string::npos) {
          string tmp = entry;
          getline(file, entry, '>');
          entry = tmp + '>' + entry;
        }
        accessions.push_back(entry.substr(0, entry.find_first_of(' ')));
        descriptors.push_back(entry.substr(entry.find_first_of(' ') + 1, entry.find_first_of('\n') - entry.find_first_of(' ') - 1));
        sequences.push_back(entry.substr(entry.find_first_of('\n') + 1, entry.size() - 1));
      }
      bytes += entry.length();
      if (interval != time(nullptr))
      {
        interval = time(nullptr);
        Rcout << "reading fasta file(s) ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
  Rcout << "reading fasta file(s) ... 100%" << std::endl;
  return DataFrame::create(_["accession"] = accessions, _["description"] = descriptors, _["sequence"] = sequences,
                           _["stringsAsFactors"] = false);
}


// [[Rcpp::export]]
void modify_fasta(const vector<string> & db, const vector<string> & proteins, std::string filename)
{
  set<string> set;
  for(auto protein : proteins)
  {
    set.insert(protein);
  }
  time_t interval = time(nullptr);
  int fsize = 0;
  int bytes = 0;
  for (auto f : db)
  {
    fsize += filesize(f.c_str());
  }
  ofstream ofile(filename.c_str());
  for (auto f : db)
  {
    ifstream file(f);
    string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      if(set.find(entry.substr(0, entry.find_first_of(' '))) != set.end())
      {
        while (entry.find('\n') == std::string::npos) {
          string tmp = entry;
          getline(file, entry, '>');
          entry = tmp + '>' + entry;
        }
        ofile << '>' << entry;
      }
      bytes += entry.length();
      if (interval != time(nullptr))
      {
        interval = time(nullptr);
        Rcout << "reading fasta file(s) ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
  ofile.close();
  Rcout << "reading fasta file(s) ... 100%" << std::endl;
}

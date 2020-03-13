#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <Rcpp.h>
#include <ctime>
#include <regex>

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

std::vector<std::string> split(const std::string & s, std::string & delimiter)
{
  std::regex rgx(delimiter);
  std::sregex_token_iterator iter(s.begin(),
                                  s.end(),
                                  rgx);
  std::vector<std::string> res;
  std::sregex_token_iterator end;
  
  for ( ; iter != end; ++iter) {
    std::cout << *iter << '\n';
    res.push_back(*iter);
  }
  return res;
}

// [[Rcpp::export]]
void trypsin_digestion(const vector<string> & files, int missed_cleavage, int min_length, int max_length) {
  time_t interval = time(nullptr);
  vector<string> accession;
  vector<string> name;
  vector<string> sequence;
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
      while (entry.find('\n') == std::string::npos) {
        string tmp = entry;
        getline(file, entry, '>');
        entry = tmp + '>' + entry;
      }
      std::string acc = entry.substr(0, entry.find_first_of(' '));
      accession.push_back(entry.substr(0, entry.find_first_of(' ')));
      std::string nm = entry.substr(entry.find_first_of(' ') + 1, entry.find_first_of('\n') - entry.find_first_of(' ') - 1);
      name.push_back(nm);
      std::string seq = entry.substr(entry.find_first_of('\n') + 1, entry.size() - 1);
      sequence.push_back(seq);
      for(int i = 1; i <= (missed_cleavage + 1); i++) {
        std::stringstream ss;
        ss << "(\\w{" << min_length << "," << max_length << "}?(?!P)[KR]){" << i << "}";
        std::string delim = ss.str();
        std::vector<std::string> sseq = split(seq, delim);
        std::cout << "i:" << i << std::endl;
        for(auto s : sseq) {
          sequence.push_back(s);
        }
      }
      Rcpp::stop("");
      bytes += entry.length();
      if (interval != time(nullptr))
      {
        interval = time(nullptr);
        Rcout << "reading fasta file(s) ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
}
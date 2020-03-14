#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <Rcpp.h>
#include <ctime>
#include <regex>
// [[Rcpp::plugins(cpp11)]]

std::ifstream::pos_type filesize(const char* filename)
{
  std::ifstream input(filename, std::ifstream::ate | std::ifstream::binary);
  return input.tellg();
}
// [[Rcpp::export]]
Rcpp::DataFrame read_fasta(const std::vector<std::string> & files, const std::vector<std::string> & proteins)
{
  std::set<std::string> set;
  for(auto protein : proteins)
  {
    set.insert(protein);
  }
  time_t interval = time(nullptr);
  std::vector<std::string> accessions;
  std::vector<std::string> descriptors;
  std::vector<std::string> sequences;
  int fsize = 0;
  int bytes = 0;
  for (auto f : files)
  {
    fsize += filesize(f.c_str());
  }
  for (auto f : files)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      if(set.find(entry.substr(0, entry.find_first_of(' '))) != set.end())
      {
        while (entry.find('\n') == std::string::npos) {
          std::string tmp = entry;
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
        Rcpp::Rcout << "reading fasta file(s) ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
  Rcpp::Rcout << "reading fasta file(s) ... 100%" << std::endl;
  return Rcpp::DataFrame::create(Rcpp::_["accession"] = accessions,
                                 Rcpp::_["description"] = descriptors,
                                 Rcpp::_["sequence"] = sequences,
                                 Rcpp::_["stringsAsFactors"] = false);
}


// [[Rcpp::export]]
void modify_fasta(const std::vector<std::string> & db, const std::vector<std::string> & proteins, std::string filename)
{
  std::set<std::string> set;
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
  std::ofstream ofile(filename.c_str());
  for (auto f : db)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      if(set.find(entry.substr(0, entry.find_first_of(' '))) != set.end())
      {
        while (entry.find('\n') == std::string::npos) {
          std::string tmp = entry;
          getline(file, entry, '>');
          entry = tmp + '>' + entry;
        }
        ofile << '>' << entry;
      }
      bytes += entry.length();
      if (interval != time(nullptr))
      {
        interval = time(nullptr);
        Rcpp::Rcout << "reading fasta file(s) ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
  ofile.close();
  Rcpp::Rcout << "reading fasta file(s) ... 100%" << std::endl;
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
void trypsin_digestion(const std::vector<std::string> & files, int missed_cleavage, int min_length, int max_length) {
  time_t interval = time(nullptr);
  std::vector<std::string> accession;
  std::vector<std::string> name;
  std::vector<std::string> sequence;
  int fsize = 0;
  int bytes = 0;
  for (auto f : files)
  {
    fsize += filesize(f.c_str());
  }
  for (auto f : files)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      while (entry.find('\n') == std::string::npos) {
        std::string tmp = entry;
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
        Rcpp::Rcout << "reading fasta file(s) ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
}
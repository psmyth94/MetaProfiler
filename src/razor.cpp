#include <Rcpp.h>
#include <vector>
#include <string>

// [[Rcpp::plugins(cpp11)]]

// void quicksort_indices(std::vector<int> &array_indices,const std::vector<std::vector<std::string>> peptides, int left ,int right)
// {
//   if(left<right)
//   {
//     int middle;
//     int x=peptides[array_indices[left]].size();
//     int l=left;
//     int r=right;
//     while(l<r)
//     {
//       while((peptides[array_indices[l]].size() <= x)&&(l<right)) l++ ;
//       while((peptides[array_indices[r]].middle_costs.objective>x)&&(r>=left)) r-- ;
//       if(l<r)
//       {
//         int temp = array_indices[l];
//         array_indices[l]=array_indices[r];
//         array_indices[r]=temp ;
//       }
//     }
//     middle=r;
//     int temp=array_indices[left];
//     array_indices[left]=array_indices[middle];
//     array_indices[middle]=temp;
//     
//     quicksort_indices_SO(array_indices,gen,left,middle-1);
//     quicksort_indices_SO(array_indices,gen,middle+1,right);
//   }
// }

// [[Rcpp::export]]
Rcpp::List razor(std::vector<std::vector<std::string>> x, std::vector<std::string> id, const bool verbose = true) 
{
  std::vector<std::vector<std::string>> group(1);
  group[0].push_back(id[0]);
  std::vector<std::vector<std::string>> cmp;
  cmp.push_back(x[0]);
  std::vector<std::vector<int>> count;
  count.push_back(std::vector<int>(x[0].size(), 1));
  std::vector<std::vector<int>> unique;
  unique.push_back(std::vector<int>(x[0].size(), 0));
  time_t last_invoke_ = time(nullptr);
  for(size_t i = 1; i < x.size(); i++)
  {
    for(size_t j = 0; j < cmp.size(); j++)
    {
      std::vector<int> c = count[j];
      std::vector<int> u = unique[j];
      for(std::vector<std::string>::iterator it = x[i].begin(); it != x[i].end(); it++)
      {
        auto found = find(cmp[j].begin(), cmp[j].end(), *it);
        if(found == cmp[j].end())
        {
          goto ctn1;
        }
        ++c[distance(cmp[j].begin(), found)];
        if(x[i].size() == 1) ++u[distance(cmp[j].begin(), found)];
      }
      group[j].push_back(id[i]);
      count[j] = c;
      unique[j] = u;
      goto ctn2;
      ctn1:;
    }
    cmp.push_back(x[i]);
    group.resize(group.size() + 1);
    group[group.size() - 1].push_back(id[i]);
    count.push_back(std::vector<int>(x[i].size(), 1));
    unique.push_back(std::vector<int>(x[i].size(), 0));
    
    ctn2:;
    if (verbose && last_invoke_ != time(nullptr))
    {
      last_invoke_ = time(nullptr);
      Rcpp::Rcout << round((double)(i + 2)/(double)x.size() * 100) << "%                     \r";
    }
  }
  if (verbose)
  {
    Rcpp::Rcout << "100%                     \n";
  }
  for(size_t i = 0; i < cmp.size(); i++)
  {
    std::vector<int> sorted_indices(cmp[i].size());
    for(size_t j=0;j<cmp[i].size();j++)
      sorted_indices[j] = j;
    std::sort(sorted_indices.begin(), sorted_indices.end(),
              [&](const int& l, const int& r) {
                int lhs = count[i][l] + unique[i][l];
                int rhs = count[i][r] + unique[i][r];
                return lhs > rhs;
              }
    );
    std::vector<std::string> tmp1(cmp[i].size());
    std::vector<int> tmp2(count[i].size());
    std::vector<int> tmp3(unique[i].size());
    for(size_t j=0;j<sorted_indices.size();j++)
    {
      tmp1[j] = cmp[i][sorted_indices[j]];
      tmp2[j] = count[i][sorted_indices[j]];
      tmp3[j] = unique[i][sorted_indices[j]];
    }
    cmp[i] = tmp1;
    count[i] = tmp2;
    unique[i] = tmp3;
  }

  return Rcpp::List::create(Rcpp::_["x"] = cmp, Rcpp::_["id"] = group, Rcpp::_["count"] = count, Rcpp::_["unique"] = unique);
}
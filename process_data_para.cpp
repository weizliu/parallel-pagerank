#include <iostream>
#include "FileNode.h"
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <queue>
#include <cmath>
#include <climits>
#include <omp.h>
using namespace std;

int total_nodes;

struct M{
  omp_lock_t l;
  M(){omp_init_lock(&l);}
  ~M(){omp_destroy_lock(&l);}
  void Lock(){omp_set_lock(&l);}
  void Unlock(){omp_unset_lock(&l);}
};

M m1;

void read_data_to_adjancy_list(string filename, map<int, vector<int> >& file_nodes){
  cout << "start" << endl;
    ifstream file(filename.c_str());
    cout << "open file over" << endl;
    string line;
    int count_for_debug = 0;
      while(getline(file, line) && count_for_debug < 100000){
        count_for_debug++;
        if(line == "") break;
        if(line[0] == '#') continue;
        int first_space_idx = line.find_first_of('\t');
        string file_id_str = line.substr(0, first_space_idx);
        string dest_id_str = line.substr(first_space_idx+1);
        //cout << "22222222222222222222222:" << file_id_str << ", " << dest_id_str << endl;
        int file_id = stoll(file_id_str);
        int dest_id = stoll(dest_id_str);
        if(total_nodes < file_id) total_nodes = file_id;
        if(total_nodes < dest_id) total_nodes = dest_id;
        if(dest_id >= 1 && file_id >= 1)
          file_nodes[file_id-1].push_back(dest_id-1);
        else
          cout << "wrong!!! id is < 0 " << dest_id << ", " << file_id << endl;
      }

      //for(int i = 0; i < total_nodes; i++){
      //  FileNode* fn = new FileNode(i);
      //  id_filenode[i] = fn;
      //}
    file.close();
    cout << "finish read file" << endl;
    cout << "total: " << total_nodes << endl;
    /*ofstream file_out("adj_list.txt");
    map<int, vector<int> >::iterator it = file_nodes.begin();
    while(it != file_nodes.end()){
      file_out << it->first << ": ";
      vector<int> dests = (it->second);
      for(int i = 0; i < dests.size(); i++){
        file_out << dests[i] << " ";
      }
      file_out << "\n";
      it++;
    }
    cout << "finish write file" << endl;*/
}

double distance(const vector<double>& v1, const vector<double>& v2){
  double res = 0;
  for(int i = 0; i < v1.size(); i++){
    res += abs(v1[i]-v2[i]);
  }
  return res;
}

//can parallel
double distance_approaching(const vector<double>& v1, const vector<double>& v2){
  double res = INT_MAX;
  for(int i = 0; i < v1.size(); i++){
    if(abs(v1[i]-v2[i]) < res){
      res = abs(v1[i]-v2[i]);
    }
  }
  return res;
}

void write(ofstream& file, vector<double> v){
  for(int i = 0; i < v.size(); i++){
    if(v[i] != 0)
    file << i << ": " << v[i] << "\n";
  }
}

vector<double> serial_pr(map<int, vector<int> >& file_nodes){
  //threshold, alpha
  vector<double> pr_new(total_nodes, 0);
  vector<double> pr_old(total_nodes, 1.0/sqrt(total_nodes));
  //total_nodes = 100;
  //ofstream file("r_single.txt");
  double threshold = 10e-6;
  double dangling_value = 0;
  double alpha = 0.85;
  double error;
  int c = 1;
  while(1){
    //c--;
    //can parallel
    //fill_n(pr_new.begin(), pr_new.end(), 0);
    int nthreads, th_id;
    int count = 0;
    #pragma omp parallel private(nthreads, th_id, count)
    {
      th_id = omp_get_thread_num();
      nthreads = omp_get_num_threads();
      //cout << nthreads << endl;
      //cout << "th_id: " << th_id << endl;
      int num_per_thread = total_nodes/nthreads;
      count = 0;
      int end = (th_id == nthreads-1) ? total_nodes : (th_id+1)*num_per_thread;
      for(int i = th_id * num_per_thread; i < end; i++){
        count++;
        pr_new[i] = 0;
      }
    }
    #pragma omp barrier

    //can parallel
   // #pragma omp parallel private(nthreads, th_id, count)
   // {
      th_id = omp_get_thread_num();
      nthreads = omp_get_num_threads();
      //cout << "000000000000000000000000000000000000: " << th_id << ", " << nthreads << endl;
      th_id = 0;
      nthreads = 1;
      int num_per_thread = total_nodes/nthreads;
      count = 0;
      int end = (th_id == nthreads-1) ? total_nodes : (th_id+1)*num_per_thread;
      vector<int> dests;
      for(int i = th_id * num_per_thread; i < end; i++){
        if(file_nodes.count(i) != 0){
          dests = file_nodes[i];
          for(int j = 0; j < dests.size(); j++){
            int d = dests[j];
            if(d >= total_nodes || d < 0){
              cout << "omg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
              continue;
            }
            //m1.Lock();
            //cout << d << ", " << i << ", " << j << endl;
            pr_new[d] += alpha * pr_old[i] / double(file_nodes[i].size());
            //m1.Unlock();
          }
        }
        else{
          //i+1 is a dangling node(no outlink)
          //m1.Lock();
          dangling_value += alpha * pr_old[i];
          //m1.Unlock();
        }
      }
  //  }
  //  #pragma omp barrier

    //can parallel
    double sum = 0;
    #pragma omp parallel private(nthreads, th_id, count)
    {
      th_id = omp_get_thread_num();
      nthreads = omp_get_num_threads();
      //cout << nthreads << endl;
      //cout << "th_id: " << th_id << endl;
      int num_per_thread = total_nodes/nthreads;
      count = 0;
      int end = (th_id == nthreads-1) ? total_nodes : (th_id+1)*num_per_thread;
      for(int i = th_id * num_per_thread; i < end; i++){
        pr_new[i] += (dangling_value + 1 - alpha) / double(total_nodes);
        //m1.Lock();
        //sum += pr_new[i]*pr_new[i];
        //m1.Unlock();
      }
    }
    #pragma omp barrier
  
    for(int i = 0; i < total_nodes; i++){
      sum += pr_new[i] * pr_new[i];
    }

    #pragma omp parallel private(nthreads, th_id, count)
    {
      th_id = omp_get_thread_num();
      nthreads = omp_get_num_threads();
      //cout << nthreads << endl;
      //cout << "th_id: " << th_id << endl;
      int num_per_thread = total_nodes/nthreads;
      count = 0;
      int end = (th_id == nthreads-1) ? total_nodes : (th_id+1)*num_per_thread;
      for(int i = th_id * num_per_thread; i < end; ++i){
        pr_new[i] /= sqrt(sum);
      }
    }
    #pragma omp barrier
    
    error = 0;
    for(int i = 0; i < pr_old.size(); i++){
      error += abs(pr_new[i]-pr_old[i]);
    }

    if(error < threshold) break;
    //cout << "error: " << error << endl;
    pr_old = pr_new;
  }

  return pr_new;
}

class com{
public:
  bool operator()(const pair<double, int>& p1, const pair<double, int>& p2){
    return p1.first < p2.first;
  }
};

int main(int argc, char* argv[]){
  map<int, vector<int> > file_nodes;
  read_data_to_adjancy_list("com-amazon.ungraph.txt", file_nodes);

  double start = omp_get_wtime();
  vector<double> pr = serial_pr(file_nodes);
  cout << "time: " << omp_get_wtime()-start << endl;

  priority_queue<pair<double, int>, vector<pair<double, int> >, com > pq;
  for(int i = 0; i < pr.size(); i++){
    pq.push(make_pair(pr[i], i+1));
  }
  cout << "result: " << endl;
  for(int i = 0; i < 10; i++){
    cout << pq.top().second << ": " << pq.top().first << endl;
    pq.pop();
  }
}


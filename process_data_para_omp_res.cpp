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

struct M{
  omp_lock_t l;
  M(){omp_init_lock(&l);}
  ~M(){omp_destroy_lock(&l);}
  void Lock(){omp_set_lock(&l);}
  void Unlock(){omp_unset_lock(&l);}
};

M m1, c;

void read_data_to_adjancy_list(string filename, map<int, vector<int> >& file_nodes, int& total_nodes){
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
        int file_id = stoll(file_id_str);
        int dest_id = stoll(dest_id_str);
        if(total_nodes < file_id) total_nodes = file_id;
        if(total_nodes < dest_id) total_nodes = dest_id;
        if(dest_id >= 1 && file_id >= 1)
          file_nodes[file_id-1].push_back(dest_id-1);
        else
          cout << "wrong!!! id is < 0 " << dest_id << ", " << file_id << endl;
      }

    file.close();
    cout << "finish read file" << endl;
    cout << "total: " << total_nodes << endl;
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

vector<double> serial_pr(map<int, vector<int> >& file_nodes, double alp, int num_threads, int total_nodes){
  vector<double> pr_new(total_nodes, 0);
  vector<double> pr_old(total_nodes, 1.0/sqrt(total_nodes));
  double threshold = 10e-9;
  double dangling_value = 0;
  double alpha = alp;
  int c = 1;
  int nthreads, th_id;
  double sum = 0;
  double error = 0;
  omp_set_num_threads(num_threads);
  #pragma omp parallel private(nthreads, th_id)
  {
	  th_id = omp_get_thread_num();
	  nthreads = omp_get_num_threads();
	  int num_per_thread = total_nodes/nthreads;
	  int count = 0;
	  int end = (th_id == nthreads-1) ? total_nodes : (th_id+1)*num_per_thread;
	  while(1){
		error = 0;
		sum = 0;
	        dangling_value = 0;
	      for(int i = th_id * num_per_thread; i < end; i++){
		count++;
		pr_new[i] = 0;
	      }
	      #pragma omp barrier

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
		    m1.Lock();
		    pr_new[d] += alpha * pr_old[i] / double(file_nodes[i].size());
		    m1.Unlock();
		  }
		}
		else{
		  //i+1 is a dangling node(no outlink)
		  m1.Lock();
		  dangling_value += alpha * pr_old[i];
		  m1.Unlock();
		}
	      }
	
	    #pragma omp barrier
	    
	    double sum_p = 0;
	    for(int i = th_id * num_per_thread; i < end; i++){
		pr_new[i] += (dangling_value + 1 - alpha) / double(total_nodes);
		sum_p += pr_new[i]*pr_new[i];
	    }
	    m1.Lock();
	    sum += sum_p;
	    m1.Unlock();
	    #pragma omp barrier
	  
	    th_id = omp_get_thread_num();
	    nthreads = omp_get_num_threads();
	    for(int i = th_id * num_per_thread; i < end; ++i){
		pr_new[i] /= sqrt(sum);
	    }
	    
	   
	    
	    double error_p = 0;
	    for(int i = th_id * num_per_thread; i < end; ++i){
	      error_p += abs(pr_new[i]-pr_old[i]);
	    }
	    m1.Lock();
	    error += error_p;
	    m1.Unlock();

	    #pragma omp barrier

	   if(error < threshold) break;
	   
	   if(th_id == 0){
	   	cout << "error: " << error << " in " << th_id << endl;
	   }
	   #pragma omp barrier
	   
	  for(int i = th_id * num_per_thread; i < end; i++){
		 pr_old[i] = pr_new[i];
	  }
	}
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
  double alp = 0.85;
  int num_threads = 1;
  if(argc == 3){
	alp = atof(argv[1]);
	num_threads = atoi(argv[2]);
  }
  else if(argc != 1){
	cout << "you don't input appropriate parameters, you can" << endl;
	cout << "1) input nothing(default alpha is 0.85, and default num of threads is 1);" << endl;
	cout << "2) input alpha and num of threads" << endl;
	cout << "alpha should be a double number in range (0, 1), and num of threads should be a number in range [1, 4]" << endl;
	return -1;
  }
  if(alp <= 0 || alp >= 1 || num_threads < 0){
	cout << "you don't input appropriate parameters, you can" << endl;
	cout << "1) input nothing(default alpha is 0.85, and default num of threads is 1);" << endl;
	cout << "2) input alpha and num of threads" << endl;
	cout << "alpha should be a double number in range (0, 1), and num of threads should be a number in range [1, 4]" << endl;
	return -1;
  }

  int total_nodes;
  map<int, vector<int> > file_nodes;
  read_data_to_adjancy_list("com-amazon.ungraph.txt", file_nodes, total_nodes);
  for(alp = 0.6; alp < 1; alp += 0.05){
    for(num_threads = 1; num_threads <= 8; num_threads++){
	cout << "alpha: " << alp << ", num_threads: " << num_threads << endl;
  	double start = omp_get_wtime();
  	vector<double> pr = serial_pr(file_nodes, alp, num_threads, total_nodes);
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
  }
}


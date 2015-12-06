#include <iostream>
#include "FileNode.h"
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <queue>
#include <cmath>
#include <climits>
using namespace std;

map<int, vector<int> > file_nodes;
map<int, FileNode*> id_filenode;
int total_nodes;
vector<double> pr_old;
vector<double> pr_new;

void read_data_to_adjancy_list(string filename){
  cout << "start" << endl;
    ifstream file(filename.c_str());
    cout << "open file over" << endl;
    string line;
    int count_for_debug = 0;
    while(getline(file, line) && count_for_debug < 100){
      //count_for_debug++;
      if(line == "") break;
      if(line[0] == '#'){
        if(line[2] == 'N'){
          string get_total = line.substr(9);
          int first_space_idx = get_total.find_first_of(' ');
          string total_str = get_total.substr(0, first_space_idx);
          total_nodes = stoll(total_str);
          cout << "there are " << total_nodes << " nodes" << endl; 
        }  
        continue;
      }
      int first_space_idx = line.find_first_of('\t');
      string file_id_str = line.substr(0, first_space_idx);
      string dest_id_str = line.substr(first_space_idx+1);
      int file_id = stoll(file_id_str);
      int dest_id = stoll(dest_id_str);
      file_nodes[file_id-1].push_back(dest_id-1);
    }

    for(int i = 0; i < total_nodes; i++){
      FileNode* fn = new FileNode(i);
      id_filenode[i] = fn;
    }
    
    file.close();
    cout << "finish read file" << endl;

    ofstream file_out("adj_list.txt");
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
    cout << "finish write file" << endl;
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

vector<double>& serial_pr(){
  //threshold, alpha
  pr_old.resize(total_nodes);
  pr_new.resize(total_nodes);
  //fill_n(pr_old.begin(), pr_old.end(), 1.0/total_nodes);
  for(int i = 0; i < total_nodes; i++){
    pr_old[i] = 1.0/sqrt(total_nodes);
  }
  double threshold = 10e-2;
  double dangling_value = 0;
  double alpha = 0.85;
  double error;
  while(1){

    //can parallel
    //fill_n(pr_new.begin(), pr_new.end(), 0);
    for(int i = 0; i < total_nodes; i++){
      pr_new[i] = 0;
    }
   
    //can parallel
    vector<int> dests;
    for(int i = 0; i < total_nodes; i++){
      if(file_nodes.count(i) != 0){
        dests = file_nodes[i];
        for(int j = 0; j < dests.size(); j++){
          int d = dests[j];
          pr_new[d] += alpha * pr_old[i] / double(file_nodes[i].size());
        }
      }
      else{
        //i+1 is a dangling node(no outlink)
        dangling_value += alpha * pr_old[i];
      }
    }

    //can parallel
    double sum = 0;
    for(int i = 0; i < total_nodes; i++){
      pr_new[i] += (dangling_value + 1 - alpha) / double(total_nodes);
      sum += pr_new[i]*pr_new[i];
    }

    for(int i = 0; i < total_nodes; ++i){
      pr_new[i] /= sqrt(sum);
    }
    error = 0;
    for(int i = 0; i < pr_old.size(); i++){
      error += abs(pr_new[i]-pr_old[i]);
    }

    if(error < threshold) break;
    cout << "error: " << error << endl;
    pr_old = pr_new;
  }
  cout << "finish main loop" << endl;
  return pr_new;
}

class com{
public:
  bool operator()(pair<double, int> p1, pair<double, int> p2){
    return p1.first < p2.first;
  }
};

int main(int argc, char* argv[]){
  read_data_to_adjancy_list("com-amazon.ungraph.txt");
  vector<double> pr = serial_pr();
  cout << "return from loop" << endl;
  cout << pr.size() << endl;
  priority_queue<pair<double, int>, vector<pair<double, int> >, com > pq;
  for(int i = 0; i < pr.size(); i++){
    pq.push(make_pair(pr[i], i+1));
  }
  cout << "finish push into pq" << endl;
  cout << pq.size() << endl;
  for(int i = 0; i < 10; i++){
    cout << pq.top().second << ": " << pq.top().first << endl;
    pq.pop();
  }
}


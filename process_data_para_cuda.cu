#include <iostream>
#include "FileNode.h"
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <queue>
#include <cmath>
#include <climits>
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

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

float distance(const vector<float>& v1, const vector<float>& v2){
  float res = 0;
  for(int i = 0; i < v1.size(); i++){
    res += abs(v1[i]-v2[i]);
  }
  return res;
}

//can parallel
float distance_approaching(const vector<float>& v1, const vector<float>& v2){
  float res = INT_MAX;
  for(int i = 0; i < v1.size(); i++){
    if(abs(v1[i]-v2[i]) < res){
      res = abs(v1[i]-v2[i]);
    }
  }
  return res;
}

void write(ofstream& file, vector<float> v){
  for(int i = 0; i < v.size(); i++){
    if(v[i] != 0)
    file << i << ": " << v[i] << "\n";
  }
}

int find_nearest_power_of_2(int num){
	int tmp = 1;
	while(tmp < num){
		tmp = tmp << 1;	
	}
	return tmp;
}

/*vector<float> serial_pr(map<int, vector<int> >& file_nodes, float alp, int num_threads, int total_nodes){
  vector<float> pr_new(total_nodes, 0);
  vector<float> pr_old(total_nodes, 1.0/sqrt(total_nodes));
  vector<M*> locks(total_nodes);
  for(int i = 0; i < total_nodes; i++){
	locks[i] = new M;
  }
  float threshold = 10e-9;
  float dangling_value = 0;
  float alpha = alp;
  int c = 1;
  int nthreads, th_id;
  float sum = 0;
  float error = 0;
  omp_set_num_threads(num_threads);
  int remain_idx = num_threads * (total_nodes/num_threads);
  int ri;
  #pragma omp parallel private(nthreads, th_id)
  {
	  th_id = omp_get_thread_num();
	  nthreads = omp_get_num_threads();
	  int num_per_thread = total_nodes/nthreads;
	  int count = 0;
	  ri = remain_idx;
	  //int end = (th_id == nthreads-1) ? total_nodes : (th_id+1)*num_per_thread;
	  int end = (th_id+1)*num_per_thread;
	  while(1){
		error = 0;
		sum = 0;
	        dangling_value = 0;
		ri = remain_idx;
		#pragma omp barrier
	      for(int i = th_id * num_per_thread; i < end; i++){
		count++;
		pr_new[i] = 0;
	      }
	      while(ri < total_nodes){
		mri.Lock();
		int i = ri;
		ri++;
		mri.Unlock();
		pr_new[i] = 0;
	      }
	      #pragma omp barrier
	      ri = remain_idx;
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
		//cout << locks.size() << ", " << d << endl;
		    locks[d]->Lock();
		    pr_new[d] += alpha * pr_old[i] / float(file_nodes[i].size());
		    locks[d]->Unlock();
		  }
		}
		else{
		  //i+1 is a dangling node(no outlink)
		  md.Lock();
		  dangling_value += alpha * pr_old[i];
		  md.Unlock();
		}
	      }
		
	     while(ri < total_nodes){
		mri.Lock();		
		int i = ri++;
		mri.Unlock();
		if(file_nodes.count(i) != 0){
		  dests = file_nodes[i];
		  for(int j = 0; j < dests.size(); j++){
		    int d = dests[j];
		    if(d >= total_nodes || d < 0){
		      cout << "omg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		      continue;
		    }
		    locks[d]->Lock();
		    pr_new[d] += alpha * pr_old[i] / float(file_nodes[i].size());
		    locks[d]->Unlock();
		  }
		}
		else{
		  //i+1 is a dangling node(no outlink)
		  md.Lock();
		  dangling_value += alpha * pr_old[i];
		  md.Unlock();
		}
	      }
	    #pragma omp barrier
	    ri = remain_idx;
	    #pragma omp barrier
	    
	    float sum_p = 0;
	    for(int i = th_id * num_per_thread; i < end; i++){
		pr_new[i] += (dangling_value + 1 - alpha) / float(total_nodes);
		sum_p += pr_new[i]*pr_new[i];
	    }
	    while(ri < total_nodes){
		mri.Lock();
		int i = ri;
		ri++;
		mri.Unlock();
		pr_new[i] += (dangling_value + 1 - alpha) / float(total_nodes);
		sum_p += pr_new[i]*pr_new[i];
	    }
	    msum.Lock();
	    sum += sum_p;
	    msum.Unlock();
	    #pragma omp barrier
	    ri = remain_idx;
	    #pragma omp barrier
	  
	    for(int i = th_id * num_per_thread; i < end; ++i){
		pr_new[i] /= sqrt(sum);
	    }
	    while(ri < total_nodes){
		mri.Lock();
		int i = ri;
		ri++;
		mri.Unlock();
		pr_new[i] /= sqrt(sum);
	    }
	    #pragma omp barrier
	    ri = remain_idx;
	    #pragma omp barrier

	    float error_p = 0;
	    for(int i = th_id * num_per_thread; i < end; ++i){
	      error_p += abs(pr_new[i]-pr_old[i]);
	    }
	    while(ri < total_nodes){
		mri.Lock();
		int i = ri++;
		mri.Unlock();
		error_p += abs(pr_new[i]-pr_old[i]);
	    }
	    merror.Lock();
	    error += error_p;
	    merror.Unlock();

	    #pragma omp barrier

	   if(error < threshold) break;
	   
	   if(th_id == 0){
	   	cout << "error: " << error << " in " << th_id << endl;
	   }
	   #pragma omp barrier
	   ri = remain_idx;
	   #pragma omp barrier
	  for(int i = th_id * num_per_thread; i < end; i++){
		 pr_old[i] = pr_new[i];
	  }
 	  while(ri < total_nodes){
		mri.Lock();
		int i = ri++;
		mri.Unlock();
		pr_old[i] = pr_new[i];
	  }
	  
	}
  }
  for(int i = 0; i < total_nodes; ++i){
	delete locks[i];
  }
  return pr_new;
}*/

//set pr_old to pr_new
__global__ void calculatePR_part0(float* pr_new, float* pr_old){
	printf("%s\n", "begin of calculate PR4");
	//printf("%d\n", filenodes->size);
	unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
	pr_old[id] = pr_new[id];
}

//basic calculate
__global__ void calculatePR_part1(float* pr_new, float* pr_old, int* nodes_range, int* nodes_dest, const size_t n, const size_t num_sides, float* dv, const float alpha){
	printf("%s\n", "begin of calculate PR1");
	//printf("%d\n", filenodes->size);
	unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(id < n){	
		int start = nodes_range[id];
		int end = nodes_range[id+1];
		
		if(end == start){ //dangling nodes
			printf("%s, %d, %d, %d\n", "this is in", id, start, end);
			atomicAdd(dv, alpha * pr_old[id]);	
			/*for(int i = 0; i < n; i++){
				printf("%s, %d, %f\n", "pr_old", i, pr_old[i]);
			}*/		
		}
		else{
			//printf("%s, %d, %d, %d\n", "this is in", id, start, end);
			for(int i = start; i < end; ++i){
				int d = nodes_dest[i];
				float tmpold = pr_old[id];
				float tmp = alpha * tmpold / (end - start);
				//printf("%s, %d, %d, %d, %s, %d, %s, %d, %s, %f, %s, %f\n", "this is in", id, start, end, "d: ", d, "pr_new + d: ", pr_new + d, "alpha * pr_old[id] / (end - start): ", tmp, "tmpold: ", tmpold);
				atomicAdd(pr_new + d, tmp);
			}
			/*for(int i = 0; i < n; i++){
				printf("%s, %d, %f\n", "pr_old", i, pr_old[i]);
			}*/
		}
		//__syncthreads();
	}
	else{
		pr_new[id] = 0;
	}
}

//add dangling value
__global__ void calculatePR_part2(float* pr_new, const size_t n, float* dv, const float alpha){
	printf("%s\n", "begin of calculate PR2");
	//printf("%d\n", filenodes->size);
	unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
	if(id < n){
		pr_new[id] += (*dv + 1 - alpha) / float(n);
	}
	else{
		pr_new[id] = 0;	
	}
}

//calculate square sum of all pr_new
__global__ void calculatePR_part3(float* pr_new, float* per_block_res, const size_t n){
	printf("%s\n", "begin of calculate PR3");	
	extern __shared__ float sdata[];
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	float x = 0;
	if(i < n) x = pr_new[i];
	sdata[threadIdx.x] = x;
	int idx = ceil(blockDim.x/2.0);	
	__syncthreads();

	while(idx != 0){
		if(i + idx < n && threadIdx.x < idx) sdata[threadIdx.x] = sdata[threadIdx.x] * sdata[threadIdx.x] + sdata[threadIdx.x + idx] * sdata[threadIdx.x + idx];
		idx /= 2;
		__syncthreads();
	}

	if(threadIdx.x == 0){
		per_block_res[blockIdx.x] = sdata[0];
		printf("%d, %f\n", blockIdx.x, per_block_res[blockIdx.x]);
	}
}

__global__ void calculatePR_part4(float* pr_new, const size_t n, const float num){
	printf("%s\n", "begin of calculate PR4");
	//printf("%d\n", filenodes->size);
	unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
	if(id < n){
		pr_new[id] /= num;
	}
	else{
		pr_new[id] = 0;	
	}
}



__global__ void block_sum(const float *input, float *per_block_res, const size_t n){
	extern __shared__ float sdata[];
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	float x = 0;
	if(i < n) x = input[i];
	sdata[threadIdx.x] = x;
	int idx = ceil(blockDim.x/2.0);	
	__syncthreads();

	while(idx != 0){
		if(i + idx < n && threadIdx.x < idx) sdata[threadIdx.x] += sdata[threadIdx.x + idx];
		idx /= 2;
		__syncthreads();
	}

	if(threadIdx.x == 0){
		per_block_res[blockIdx.x] = sdata[0];
		printf("%d, %f\n", blockIdx.x, per_block_res[blockIdx.x]);
	}
}

__global__ void setValue(float* input, const size_t n, float num){
	unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
	if(id < n) input[id] = num;
	else input[id] = 0;
}

__global__ void diffPr(float* input1, float* input2, float* diff, const size_t n){
	unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
	if(id < n) diff[id] = input1[id] > input2[id] ? (input1[id] - input2[id]) : (input2[id] - input1[id]);
	else diff[id] = 0;
}

class com{
public:
  bool operator()(const pair<float, int>& p1, const pair<float, int>& p2){
    return p1.first < p2.first;
  }
};

int main(int argc, char* argv[]){
	float alp = 0.85;
	int num_threads = 1;
	if(argc == 3){
	alp = atof(argv[1]);
	num_threads = atoi(argv[2]);
	}
	else if(argc != 1){
	cout << "you don't input appropriate parameters, you can" << endl;
	cout << "1) input nothing(default alpha is 0.85, and default num of threads is 1);" << endl;
	cout << "2) input alpha and num of threads" << endl;
	cout << "alpha should be a float number in range (0, 1), and num of threads should be a number in range [1, 4]" << endl;
	return -1;
	}
	if(alp <= 0 || alp >= 1 || num_threads < 0){
	cout << "you don't input appropriate parameters, you can" << endl;
	cout << "1) input nothing(default alpha is 0.85, and default num of threads is 1);" << endl;
	cout << "2) input alpha and num of threads" << endl;
	cout << "alpha should be a float number in range (0, 1), and num of threads should be a number in range [1, 4]" << endl;
	return -1;
	}

	int total_nodes = 6;
	float threshold = 10e-6;
	float alpha = 0.85;
	//map<int, vector<int> >* file_nodes;
	//(*file_nodes)[1].push_back(2);	
	//vector<FileNode*> file_nodes;
	
	//read_data_to_adjancy_list("com-amazon.ungraph.txt", file_nodes, total_nodes);	
	
	
	/*FileNode* file_nodes[total_nodes];
		

	FileNode* a0 = new FileNode(0);
	a0->dest.push_back(1);
	a0->size = 10;
	file_nodes[0] = a0;
	FileNode* a1 = new FileNode(1);
	file_nodes[1] = a1;	
	a1->size = 20;
	cout << "file_nodes[1].id: " << file_nodes[1]->id << endl;*/

	//copy file_nodes to GPU
	/*FileNode** d_file_nodes;
	for(int i = 0; i < total_nodes; ++i){
		file_nodes[i] = new FileNode(i);		
		FileNode* fn = file_nodes[i];
		FileNode* d_fn;
		cudaMalloc((void**)&d_fn, sizeof(FileNode*));
		cudaMemcpy(d_fn, fn, sizeof(FileNode*), cudaMemcpyHostToDevice);
		file_nodes[i] = d_fn;
	}
	cudaMalloc((void**)&d_file_nodes, sizeof(FileNode*) * total_nodes);
	cudaMemcpy(d_file_nodes, file_nodes, sizeof(FileNode*) * total_nodes, cudaMemcpyHostToDevice);*/
	
	int num_sides = 10;
	int nodes_range[total_nodes + 1]; //nodes_range[idx] = start_idx means dest of node idx starts from nodes_dest[start_idx] to nodes_dest[nodes_range[idx+1]-1]; last ele should be num_sides
	int nodes_dest[num_sides]; //dest nodes
	int *d_nodes_range, *d_nodes_dest;
	nodes_range[0] = 0; nodes_range[1] = 3; nodes_range[2] = 6; nodes_range[3] = 9; nodes_range[4] = 9; nodes_range[5] = 10; nodes_range[6] = 10; 
	nodes_dest[0] = 1; nodes_dest[1] = 3; nodes_dest[2] = 5; nodes_dest[3] = 0; nodes_dest[4] = 2; nodes_dest[5] = 4; nodes_dest[6] = 1; nodes_dest[7] = 4; nodes_dest[8] = 5; nodes_dest[9] = 2; //nodes_dest[10] = 1; nodes_dest[11] = 2; nodes_dest[12] = 3; nodes_dest[13] = 4; 
	cudaMalloc((void**)&d_nodes_range, sizeof(int) * (1 + total_nodes));
	cudaMemcpy(d_nodes_range, nodes_range, sizeof(int) * (1 + total_nodes), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&d_nodes_dest, sizeof(int) * num_sides);
	cudaMemcpy(d_nodes_dest, nodes_dest, sizeof(int) * num_sides, cudaMemcpyHostToDevice);
	
	int block_size = 2;
	int block_num = (int)ceil(total_nodes/float(block_size));
	cout << "block_num is " << block_num << endl; 

	//allocate d_pr_old, d_pr_new on GPU
	float *d_pr_old, *d_diff, *d_pr_new;
	float pr_old[block_size * block_num], pr_new[block_size * block_num];
	cudaMalloc((void**)&d_pr_old, sizeof(float) * block_num * block_size); 
	cudaMalloc((void**)&d_pr_new, sizeof(float) * block_num * block_size); 
	cudaMalloc((void**)&d_diff, sizeof(float) * block_num * block_size);

	//set pr_old to 1/sqrt(n) using GPU and pr_new to 0
	float num = 1.0/sqrt(total_nodes);
	setValue<<<block_num, block_size>>>(d_pr_new, total_nodes, num);
	
	float debug_pr_old[block_size * block_num], debug_pr_new[block_size * block_num];
	cudaMemcpy(debug_pr_new, d_pr_new, sizeof(float) * block_num * block_size, cudaMemcpyDeviceToHost);
	cout << "pr_pre: " << endl;
	for(int i = 0; i < block_num * block_size; ++i){
		cout << i << ": " << debug_pr_new[i] << endl;
	}
	num = 0;
	setValue<<<block_num, block_size>>>(d_pr_old, total_nodes, num);
	
	while(1){
		//calculate if the condition now can exit
		diffPr<<<block_num, block_size>>>(d_pr_old, d_pr_new, d_diff, total_nodes);
		
		float debug_diff[block_num * block_size];
		cudaMemcpy(debug_diff, d_diff, sizeof(float) * block_num * block_size, cudaMemcpyDeviceToHost);
		cout << "diff: " << endl;
		for(int i = 0; i < block_num * block_size; ++i){
			cout << i << ": " << debug_diff[i] << endl;
		}
		cout << endl;

		float *d_sum, sum, *d_partial_sum;
		cudaMalloc((void**)&d_sum, sizeof(float));	
		cudaMalloc((void**)&d_partial_sum, sizeof(float) * (block_num));	
		block_sum<<<block_num, block_size, block_size * sizeof(float)>>>(d_diff, d_partial_sum, total_nodes);

		float debug_ps[block_num];
		cudaMemcpy(debug_ps, d_partial_sum, sizeof(float) * block_num, cudaMemcpyDeviceToHost);
		cout << "partial sum: " << endl;
		for(int i = 0; i < block_num ; ++i){
			cout << i << ": " << debug_ps[i] << endl;
		}

		cout << "block_num and find_nearest_power_of_2(block_num): " << block_num << ", " << find_nearest_power_of_2(block_num) << endl;
		block_sum<<<1, find_nearest_power_of_2(block_num), find_nearest_power_of_2(block_num) * sizeof(float)>>>(d_partial_sum, d_sum, block_num);
		cudaMemcpy(&sum, d_sum, sizeof(float), cudaMemcpyDeviceToHost);
		cout << "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeerror:" << sum << endl;
		if(sum < threshold) break; //if distance smaller than threshold, done!
		//break;	
		
		//finish one iterations
		cout << "before cPR0" << endl;
		calculatePR_part0<<<block_num, block_size>>>(d_pr_new, d_pr_old);
		cout << "after cPR0" << endl;

		cudaMemcpy(debug_pr_old, d_pr_old, sizeof(float) * block_num * block_size, cudaMemcpyDeviceToHost);
		cout << "pr_old: " << endl;
		for(int i = 0; i < block_num * block_size; ++i){
			cout << i << ": " << debug_pr_old[i] << endl;
		}

		//set pr_new to 0
		num = 0;
		setValue<<<block_num, block_size>>>(d_pr_new, total_nodes, num);

		//set dangling value
		float dangling_value = 0;
		float* d_dv;
		cudaMalloc((void**)&d_dv, sizeof(float));
		cudaMemcpy(d_dv, &dangling_value, sizeof(float), cudaMemcpyHostToDevice);
		
		cudaMemcpy(debug_pr_old, d_pr_old, sizeof(float) * block_num * block_size, cudaMemcpyDeviceToHost);
		cout << "pr_old: " << endl;
		for(int i = 0; i < block_num * block_size; ++i){
			cout << i << ": " << debug_pr_old[i] << endl;
		}
	
		cout << "before cPR1" << endl;
		
		calculatePR_part1<<<block_num, block_size, block_size * sizeof(float)>>>(d_pr_new, d_pr_old, d_nodes_range, d_nodes_dest, total_nodes, num_sides, d_dv, alpha);
		
		cout << "after cPR1" << endl;
		
		cudaMemcpy(&dangling_value, d_dv, sizeof(float), cudaMemcpyDeviceToHost);
		cout << "dangling_value: " << dangling_value << endl;
		cudaMemcpy(debug_pr_new, d_pr_new, sizeof(float) * block_num * block_size, cudaMemcpyDeviceToHost);
		cout << "pr_new: " << endl;
		for(int i = 0; i < block_num * block_size; ++i){
			cout << i << ": " << debug_pr_new[i] << endl;
		}
		
		cout << "before cPR2" << endl;
		calculatePR_part2<<<block_num, block_size>>>(d_pr_new, total_nodes, d_dv, alpha);
		cout << "after cPR2" << endl;
	
		cudaMemcpy(debug_pr_new, d_pr_new, sizeof(float) * block_num * block_size, cudaMemcpyDeviceToHost);
		cout << "pr_new: " << endl;
		for(int i = 0; i < block_num * block_size; ++i){
			cout << i << ": " << debug_pr_new[i] << endl;
		}
		
		float *d_square_sum, square_sum, *d_partial_square_sum;
		cudaMalloc((void**)&d_square_sum, sizeof(float));	
		cudaMalloc((void**)&d_partial_square_sum, sizeof(float) * (block_num));
	
		cout << "before cPR3" << endl;
		calculatePR_part3<<<block_num, block_size, block_size * sizeof(float)>>>(d_pr_new, d_partial_square_sum, total_nodes);
		block_sum<<<1, find_nearest_power_of_2(block_num), find_nearest_power_of_2(block_num) * sizeof(float)>>>(d_partial_square_sum, d_square_sum, block_num);
		cout << "after cPR3" << endl;

		cudaMemcpy(&square_sum, d_square_sum, sizeof(float), cudaMemcpyDeviceToHost);
		cout << "square_sum: "  << square_sum << endl;
		
		cout << "before cPR4" << endl;
		float num = sqrt(square_sum);
		calculatePR_part4<<<block_num, block_size>>>(d_pr_new, total_nodes, num);
		cout << "after cPR4" << endl;
		
		cudaMemcpy(debug_pr_new, d_pr_new, sizeof(float) * block_num * block_size, cudaMemcpyDeviceToHost);
		cout << "pr_new: " << endl;
		for(int i = 0; i < block_num * block_size; ++i){
			cout << i << ": " << debug_pr_new[i] << endl;
		}

		

		//break;
	}
	
	//cudaMalloc((void**)&d_diff, sizeof(float) * block_num * block_size);
	//cudaMemcpy(d_diff, pr_old, sizeof(float) * block_num * block_size, cudaMemcpyHostToDevice);

	

	//float threshold = 10e-9;
	
	//pr<<<block_num, block_size, block_size * sizeof(float)>>>(d_file_nodes, total_nodes, threshold, d_pr_old, d_diff);

	//cudaMemcpy(file_nodes, d_file_nodes, sizeof(file_nodes) * total_nodes, cudaMemcpyDeviceToHost);
	//cout << file_nodes[0].getId() << endl;















	/*priority_queue<pair<float, int>, vector<pair<float, int> >, com > pq;
	for(int i = 0; i < pr.size(); i++){
	pq.push(make_pair(pr[i], i+1));
	}
	cout << "result: " << endl;
	for(int i = 0; i < 10; i++){
	cout << pq.top().second << ": " << pq.top().first << endl;
	pq.pop();
	}*/
	cudaFree(d_pr_new);
	cudaFree(d_pr_old);
	cudaFree(d_diff);
	//cudaFree(d_file_nodes);
}


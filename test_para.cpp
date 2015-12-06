#include "iostream"
#include <omp.h>
#include <vector>
using namespace std;

int main(){
  int nthreads, tid;
  int total = 12500;
  vector<int> a(total, 0);
  double start = omp_get_wtime();
int c = 100;
#pragma omp parallel
  {
	while(c){
		c--;
    		int nthreads = omp_get_num_threads();
    		int tid = omp_get_thread_num();
    		int num_per_t = total/nthreads;
    		int end = (tid == nthreads-1) ? total : ((tid+1) * num_per_t);
		int count = 0;
    		for(int i = tid * num_per_t; i < end; ++i){
      		a[i] = i;
		count++;
    	} 
	cout << tid << ": " << count << endl;
  	//#pragma omp barrier
	}
}
  cout << "time: " << omp_get_wtime()-start << endl;
}

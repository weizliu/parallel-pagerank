#include "iostream"
#include <vector>

using namespace std; 

class FileNode{
public:
  int id;  
  vector<int> dest;
  int size;
  FileNode() : size(0){}
  FileNode(int i) : id(i), size(0){}
};

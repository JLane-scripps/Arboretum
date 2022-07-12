// KD TREE in C++

#include <cmath>
using namespace std;


// Adjust K to have as many points as needed. Note: Capital K
const float K = 3.0; 

// A structure for PSMs, which will form the KD Tree
struct PSM {
    float datapoints[K];       // Each PSM will have K datapoints
                        // Each PSM is made of an int array with K slots
    PSM* left, * right;  // Pointers to connect children PSM's
};

// Creating PSMs
// Pass into this function an array of values to assign to a new PSM
struct PSM* newPSM(float arr[]) {           // A struct for a PSM **Pointer** that points to a new float array *YOU PASS IN*
    struct PSM* temp = new PSM;             // Another struct for a PSM Pointer called "temp" that creates a new PSM.
    for (int i = 0; i < K; i++) {           // Start at 0, iterate up to (technically 1 less than) the number of datapoints
        temp->datapoints[i] = arr[i];       // "temp"'s datapoints all take on the value in corresponding array YOU PASSED IN
    }
    temp->left = temp->right = NULL;        // Direct new PSM (named "temp")'s pointers to point to NULL for now
    return temp;                            // return the new PSM
}
// QUESTION: Double check this is not returning an array, or expecting an array back
// Is an array the best way to pass this information in / out? How will a PSM be 


// Adds PSM's to the tree
// Special Case: Tree is empty, there are no Nodes.
// Base Case: Add new PSM as either left or right child
// NOTE: ALWAYS ADDS. Do not call unless certain a PSM is ready to be added here.
PSM* add(PSM* root, float datapoint[], float Depth) {
    
    // IF TREE IS EMPTY
    if (root == NULL){
        return newPSM(datapoint[]);
    }

    // Calculate current dimension of comparison
    unsigned current_dimension = Depth % K;

    // Compare the new data with root on current dimension. Decide left or right     
    if (datapoint[current_dimension] < (root->datapoint[current_dimension])){
        root->left = newPSM(root->left, datapoint, Depth+1);
    } else {
        root->right = newPSM(root->right, datapoint, Depth+1);
    }

    return root;
}

// Are these 2 Datapoints the same?
// Parameters: established float datapoint array, new float datapoint array. (Flexible order)
bool arePointsSame(float datapoint1[], float datapoint2[]) {
    for (int i = 0; i < K; i++){
        if (datapoint1[i] != datapoint2[i]) {
            return false;
        }
    }
    return true;
}

// Searches within a single datapoint array.
// Parameters: set a Node pointer to be root, a datapoint array, and the depth of the current position.
// Special Case: Tree is empty.
// Base Case: Tree has 1 PSM.
// Normal Case: Multiple PSM's.
bool searchPSM(PSM* root, float datapoint[], float depth){

    // Special Case
    if (root == NULL) {
        return false;
    }
    // Base Case
    if (arePointsSame(root->datapoints, datapoint)) {
        return true;
    }

    // Current dimension: which value we are comparing
    unsigned current_dimension = depth % K;

    // Compare root's & searched-for values for the correct dimension
    if (datapoint[current_dimension] < root->datapoint[current_dimension]) {
        return searchPSM(root->left, datapoint, depth+1);
    }

    return searchPSM(root->right, datapoint, depth+1);
}

//Search for a PSM in a KD Tree
bool search(PSM* root, float datapoint[]) {
    
    // Use searchPSM. Pass current depth as 0
    return searchPSM(root, datapoint, 0);
}

// Driver program to test above functions
int main() {
    struct PSM *root = NULL;
    float points[][K] = {{3, 6}, {17, 15}. {13, 15}, {6, 12}, {9, 1}, {2, 7}, {10, 19}};

    int n = sizeof(datapoints)/sizeof(datapoints[0]);
    
    for (int i=0; i<n; i++){
        root = insert(root, datapoints[i]);
    }

    float datapoint1[] = {10, 19};
    (search(root, datapoint1[i]))? cout << "Found.\n": cout << "Not Found.\n";

    float datapoint2[] = {12, 19};
    (search(root, datapoint2[i]))? cout << "Found.\n": cout << "Not Found.\n";

    return 0;
}
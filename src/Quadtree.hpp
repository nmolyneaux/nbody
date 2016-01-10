/* ------------------------------------------------

Author: Nicholas Molyneaux
Date : 28 November 2015

Quad tree approach to the n-body problem. Implemented
for the course of program parallelization for PC clusters

 ------------------------------------------------*/

#include <vector>

// Class defining the body, which are stored in the leaf nodes
class Body
{
public:
    // Constructor and destructor
    Body();
    Body(double x, double y, double m, double vel_x, double vel_y);
    ~Body();
    
    // Attributes
    double pos_x, pos_y, mass, vel_x, vel_y, acc_x, acc_y;   
    int proc;

    // Member functions
    void printBody(std::ostream & os);
};

// Definition of the node class. Nodes can either be simple nodes or leaves
class Node
{
public:
    // Constructor and destructor, with deletion method
    Node();
    Node(double x, double y, double h, double w);
    Node(Node *parent_ptr, double x, double y, double h, double w);
    ~Node();
    void deleteNode();

    // Attributes   
    double center_x, center_y, width, height;
    double pos_x, pos_y, mass, vel_x, vel_y, nb_bodies;
    int depth;
    bool is_leaf;
    bool contains_body;
    Body body_local;
    Node* parent;
    Node* tr;
    Node* tl;
    Node* bl;
    Node* br;    

    // Member functions
    void setBody(Body &body);
    void updateNodeOnInsert(Body &body);

};

// Quad Tree definition. Stores nodes and leaves, and contains the methods for traversing the tree.
class Quadtree
{
public:
    // Constructor and destructor, with empty method
    Quadtree();
    Quadtree(double x, double y, double w, double h, double time_step);
    Quadtree(double x, double y, double w, double h, double time_step, int procs, double precision);
    ~Quadtree();
    void empty();

    // Attributes
    double dt, theta;
    Node root;

    int min_bodies_per_node, max_bodies_per_node;   
    int min_block_size, max_block_size, nb_proc;   
    
    // Member functions
    void updateRoot(double w, double h);

    void findQuadrant(Body &body, Node &node);
    void insertBody(Body &body, Node &local_root);
    void createLeaves(Node &node);
    
    double calculateDistance(Body &body, Node &node);
    double calculateDistance(Body &body, Body &body_effect); 

    void calculateForce(Body &body, Node & node);
    void calculateForce(Body &body, Body &body_effect);
    
    //void calculateAccelerations();

    void calculateAllForcesBody(Body & body, Node & node);

    void calculateForcesInBranch(Node & node);
    void calculateForcesInBranch(Node & node, int my_rank);

    void printPositions(Node &node, double time, std::ofstream &file);
    //void updateNodesAfterMove(Node &node, Body & body);

    void moveBodies(Node &node);
    void collectBodies(std::vector<double> &bodies, Node & node, bool with_proc);

    std::vector<std::vector<Node*> > findLocalNodes(int *bodies_per_node);
    void assignNode(int *bodies_per_node, std::vector<std::vector<Node*> > &node_assignment, Node &node);

};

// Overloading the << operator to print trees, bodies and nodes
std::ostream& operator<< (std::ostream & out, Quadtree const& qt);
std::ostream& operator<< (std::ostream & out, Node const& node);
std::ostream& operator<< (std::ostream & out, Body const& body);



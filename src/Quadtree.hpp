/* ------------------------------------------------

Author: Nicholas Molyneaux
Date : 28 November 2015

Quad tree approach to the n-body problem. Implemented
or the course of program parallelization for PC clusters

 ------------------------------------------------*/

#include <vector>

class Body
{
public:
    Body();
    Body(double x, double y, double m, double vel_x, double vel_y);

    ~Body();
    double pos_x, pos_y, mass, vel_x, vel_y, acc_x, acc_y;   

    void printBody(std::ostream & os);
};

class Node
{
public:
    Node();
    Node(double x, double y, double h, double w);
    Node(Node *parent_ptr, double x, double y, double h, double w);
    ~Node();
    void deleteNode();

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

    void setBody(Body &body);
    void updateNodeOnInsert(Body &body);

};

class Quadtree
{
public:
    Quadtree();
    Quadtree(double x, double y, double w, double h, double time_step);
    ~Quadtree();
   
    void empty();
    void updateRoot(double w, double h);


    double dt;
    Node root;

    int min_bodies_per_node, max_bodies_per_node;   
    int min_block_size, max_block_size;   
    
    void findQuadrant(Body &body, Node &node);
    void insertBody(Body &body, Node &local_root);
    void createLeaves(Node &node);
    
    double calculateDistance(Body &body, Node &node);
    double calculateDistance(Body &body, Body &body_effect); 

    void calculateForce(Body &body, Node & node);
    void calculateForce(Body &body, Body &body_effect);
    
    void calculateAccelerations();

    void calculateAllForcesBody(Body & body, Node & node);
    void calculateForcesInBranch(Node & node);

    void printPositions(Node &node, double time, std::ofstream &file);
    void updateNodesAfterMove(Node &node, Body & body);

    void moveBodies(Node &node);

    void collectBodies(std::vector<double> &bodies, Node & node);

    std::vector<std::vector<Node*> > findLocalNodes(int rank, int nb_procs);
    void assignNode(std::vector<int> &bodies_per_node, std::vector<std::vector<Node*> > &node_assignment, Node &node);

};

std::ostream& operator<< (std::ostream & out, Quadtree const& qt);
std::ostream& operator<< (std::ostream & out, Node const& node);
std::ostream& operator<< (std::ostream & out, Body const& body);



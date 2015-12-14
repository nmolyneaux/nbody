/* ------------------------------------------------

Author: Nicholas Molyneaux
Date : 28 November 2015

Quad tree approach to the n-body problem. Implemented
or the course of program parallelization for PC clusters

 ------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <math.h>
#include "Quadtree.hpp"

// -------------------- BODY class ---------------------
Body::Body(){};

Body::Body(double x, double y, double m, double vx, double vy)
{
    pos_x = x;
    pos_y = y;
    mass = m;
    vel_x = vx;
    vel_y = vy;
    acc_x = 0;
    acc_y = 0;
};
Body::~Body(){};

void Body::printBody(std::ostream &os)
{
    //os << pos_x << "," << pos_y << "," << vel_x << "," << vel_y << "," << acc_x << "," << acc_y << "," <<  mass << std::endl;
    os << pos_x << "," << pos_y << std::endl;
};

// -------------------- NODE class ---------------------
Node::Node(){};

Node::Node(double x, double y, double w, double h)
{
    // Quadrant information
    center_x = x;
    center_y = y;
    width = w;
    height = h;
    is_leaf = true;
    contains_body = false;
    depth = 0;

    // bodies average
    nb_bodies = 0;
    pos_x = 0;
    pos_y = 0;
    mass = 0;
    vel_x = 0;
    vel_y = 0;
};

Node::Node(Node *parent_ptr, double x, double y, double w, double h)
{
    // Quadrant information
    center_x = x;
    center_y = y;
    width = w;
    height = h;
    is_leaf = true;
    contains_body = false;
    parent = parent_ptr;
    depth = parent->depth + 1;

    // bodies average
    nb_bodies = 0;
    pos_x = 0;
    pos_y = 0;
    mass = 0;
    vel_x = 0;
    vel_y = 0;
};

Node::~Node(){};

void Node::deleteNode()
{
    //std::cout << *this << std::endl;
    {
	if (!is_leaf)
	{
	    tr->deleteNode();
	    tl->deleteNode();
	    bl->deleteNode();
	    br->deleteNode();
	    is_leaf = true;
	}
	else if (is_leaf)
	{
	    if (depth != 0)
		delete this;
	}
    }
}

void Node::setBody(Body &body)
{
    body_local = body;
    pos_x = body.pos_x;
    pos_y = body.pos_y;
    mass = body.mass;
    //vel_x = body.vel_x;
    //vel_y = body.vel_y;
    contains_body = true;
};

void Node::updateNodeOnInsert(Body &body)
{
    nb_bodies += 1;
    pos_x = pos_x * (mass)/(mass+body.mass) + body.mass/(mass+body.mass) * (body.pos_x);
    pos_y = pos_y * (mass)/(mass+body.mass) + body.mass/(mass+body.mass) * (body.pos_y);
    mass = mass * (nb_bodies-1)/nb_bodies + (body.mass)/nb_bodies;
    //vel_x = vel_x * (nb_bodies-1)/nb_bodies + (body.vel_x)/nb_bodies;
    //vel_y = vel_y * (nb_bodies-1)/nb_bodies + (body.vel_y)/nb_bodies;
};

// -------------------- QUADTREE class ---------------------
Quadtree::Quadtree(){};

// Constructor of the QT, with the size of the universe
Quadtree::Quadtree(double x, double y, double w, double h, double time_step, int procs)
{
    nb_proc = procs;
    dt = time_step;
    root = Node(x, y, w, h);
};

Quadtree::Quadtree(double x, double y, double w, double h, double time_step)
{    
    dt = time_step;
    root = Node(x, y, w, h);
};

Quadtree::~Quadtree(){};

void Quadtree::updateRoot(double w, double h)
{
    root.width=w;
    root.height=h;
};

void Quadtree::empty()
{
    root.deleteNode();

    root.is_leaf = true;
    root.mass = 0.0;
    root.pos_x = 0.0;
    root.pos_y = 0.0;
    root.nb_bodies = 0;

}

// Finds in which quadrant of a node the body should go, then inserts it
void Quadtree::findQuadrant(Body &body, Node &node)
{
    if (body.pos_x < node.center_x)
    {
	if (body.pos_y < node.center_y)
	{
	    insertBody(body, *(node.bl));
	}
	else
	{
	    insertBody(body, *(node.tl));
	}
    }
    else
    {
	if (body.pos_y < node.center_y)
	{
	    insertBody(body, *(node.br));
	}
	else
	{
	    insertBody(body, *(node.tr));
	}
    }
};

// Actually inserts the body in the tree, if at empty leaf,
// otherwise moves down one level and calls findQuadrant again.
void Quadtree::insertBody(Body &body, Node &local_root)
{    
    // If the node is a leaf and contains a body, then leaves must be created and then 
    // the old body and new body must be replaced
    if (local_root.is_leaf && local_root.contains_body)
    {
	//std::cout << body << " in first" << std::endl;
	createLeaves(local_root);
	local_root.updateNodeOnInsert(local_root.body_local);
	findQuadrant(local_root.body_local, local_root);
	local_root.updateNodeOnInsert(body);
	findQuadrant(body, local_root);
    }
    // if an empty leaf, place body in node
    else if (local_root.is_leaf && !local_root.contains_body)
    {
	//std::cout << body << " in second" << std::endl;
	local_root.setBody(body);
    }
    // if not a leaf, move towards leaves by searching again for quadrant
    else if (!local_root.is_leaf)
    {
	//std::cout << body <<  " in third" << std::endl;
	local_root.updateNodeOnInsert(body);
	findQuadrant(body, local_root);
    }
};

// Populates the leaves for a node, and changes the status of the is_leaf and contains_body
void Quadtree::createLeaves(Node &node)
{
    node.is_leaf = false;
    node.contains_body = false;
    
    node.tr = new Node(&node, node.center_x + 0.25*node.width, node.center_y + 0.25*node.height, 0.5*node.height, 0.5*node.width);
    node.tl = new Node(&node, node.center_x - 0.25*node.width, node.center_y + 0.25*node.height, 0.5*node.height, 0.5*node.width);
    node.bl = new Node(&node, node.center_x - 0.25*node.width, node.center_y - 0.25*node.height, 0.5*node.height, 0.5*node.width);
    node.br = new Node(&node, node.center_x + 0.25*node.width, node.center_y - 0.25*node.height, 0.5*node.height, 0.5*node.width);
};

double Quadtree::calculateDistance(Body &body, Node &node)
{
    return sqrt(pow(body.pos_x - node.pos_x,2.0) + pow(body.pos_y - node.pos_y,2.0));
};

double Quadtree::calculateDistance(Body &body, Body &body_effect)
{
    return sqrt(pow(body.pos_x - body_effect.pos_x,2.0) + pow(body.pos_y - body_effect.pos_y,2.0));
};

void Quadtree::calculateForce(Body &body, Node &node)
{
    const double G = 6.674e-11;
    double r = calculateDistance(body, node);
    body.acc_x += (node.pos_x - body.pos_x) * G * node.mass / pow(r,3);
    body.acc_y += (node.pos_y - body.pos_y) * G * node.mass / pow(r,3);
};

void Quadtree::calculateForce(Body &body, Body &body_effect)
{
    const double G = 6.674e-11;
    double r = calculateDistance(body, body_effect);    
    body.acc_x += (body_effect.pos_x - body.pos_x) * G * body_effect.mass / pow(r,3);
    body.acc_y += (body_effect.pos_y - body.pos_y) * G * body_effect.mass / pow(r,3);
};

void Quadtree::updateNodesAfterMove(Node &node, Body & body)
{
    if (node.depth != 0)
    {
	node.pos_x = (dt * body.vel_x + dt*dt * body.acc_x)*body.mass/(node.mass+body.mass);
	node.pos_y = (dt * body.vel_y + dt*dt * body.acc_y)*body.mass/(node.mass+body.mass);
	updateNodesAfterMove(*node.parent, body);
    }   
}

void Quadtree::moveBodies(Node &node)
{
    if (node.is_leaf && node.contains_body)
    {
	node.body_local.pos_x += dt * node.body_local.vel_x + dt*dt * node.body_local.acc_x;
	node.body_local.pos_y += dt * node.body_local.vel_y + dt*dt * node.body_local.acc_y;
	node.body_local.vel_x += dt * node.body_local.acc_x;
	node.body_local.vel_y += dt * node.body_local.acc_y;
	node.body_local.acc_x = 0.0;
	node.body_local.acc_y = 0.0;
	//updateNodesAfterMove(*node.parent, node.body_local);
    }
    else if (!node.is_leaf)
    {
	moveBodies(*node.tr);
	moveBodies(*node.tl);
	moveBodies(*node.bl);
	moveBodies(*node.br);
    }
};

void Quadtree::calculateAllForcesBody(Body &body, Node &node)
{
    double theta = 0.5;
    double distance = calculateDistance(body, node);
    if (!node.is_leaf && (0.5*node.width + 0.5*node.height)/distance < theta)
	calculateForce(body, node);
    else if (!node.is_leaf && (0.5*node.width + 0.5*node.height)/distance > theta)    
    {
	calculateAllForcesBody(body, *node.tr);
	calculateAllForcesBody(body, *node.tl);
	calculateAllForcesBody(body, *node.bl);
	calculateAllForcesBody(body, *node.br);
    }
    else if (node.is_leaf && node.contains_body)
    {
	if (calculateDistance(body, node.body_local) != 0.0)
	    calculateForce(body, node.body_local);
    }
};

void Quadtree::calculateForcesInBranch(Node &node)
{
    if (node.contains_body && node.is_leaf)
    {
	calculateAllForcesBody(node.body_local, root);
    }
    else if (!node.is_leaf)
    {
	calculateForcesInBranch(*node.tr);
	calculateForcesInBranch(*node.tl);
	calculateForcesInBranch(*node.bl);
	calculateForcesInBranch(*node.br);
    }
}

void Quadtree::collectBodies(std::vector<double> &bodies, Node & node)
{
    if (node.contains_body && node.is_leaf)
    {
	bodies.push_back(node.body_local.pos_x);
	bodies.push_back(node.body_local.pos_y);
	bodies.push_back(node.body_local.mass);		
	bodies.push_back(node.body_local.vel_x);
	bodies.push_back(node.body_local.vel_y);	
    }
    else if (!node.is_leaf)
    {
	collectBodies(bodies, *node.tr);
        collectBodies(bodies, *node.tl);
	collectBodies(bodies, *node.bl);
	collectBodies(bodies, *node.br);
    }
}

void Quadtree::printPositions(Node &node, double time, std::ofstream &file)
{
    if (node.is_leaf && node.contains_body)
    {
	file << time << ",";
	node.body_local.printBody(file);
    }
    else if (!node.is_leaf)
    {
	printPositions(*node.tr, time, file);
	printPositions(*node.tl, time, file);
	printPositions(*node.bl, time, file);
	printPositions(*node.br, time, file);
    }
}

std::vector<std::vector < Node *> > Quadtree::findLocalNodes(int *bodies_per_node)
{
    max_block_size = 1000;//int(1.05 * (root.nb_bodies)/(nb_proc));
    min_block_size = 1000;//int(0.95 * (root.nb_bodies)/(nb_proc));
    min_bodies_per_node = int(0.95 * root.nb_bodies/(nb_proc));
    max_bodies_per_node = int(1.05 * root.nb_bodies/(nb_proc));

    for (int i = 0; i < nb_proc; i++)
	*(bodies_per_node+i) = 0;

    std::vector<std::vector< Node *> > node_assignment(nb_proc);
    assignNode(bodies_per_node, node_assignment, root);
    return node_assignment;
}

void Quadtree::assignNode(int *bodies_per_node, std::vector<std::vector<Node*> > &node_assignment, Node &node)
{
    
    if (node.nb_bodies <= max_block_size && !node.is_leaf)
    {
	int i = 0;
	while((( *(bodies_per_node+i) + node.nb_bodies) > max_bodies_per_node) && i < (nb_proc - 1) )
	{
	    i++;	    
	}
	node_assignment[i].push_back(&node);
	*(bodies_per_node+i) += node.nb_bodies;
    }
    else if (node.nb_bodies > max_block_size && !node.is_leaf)
    {
	assignNode(bodies_per_node, node_assignment, *node.tr);
	assignNode(bodies_per_node, node_assignment, *node.tl);
	assignNode(bodies_per_node, node_assignment, *node.bl);
	assignNode(bodies_per_node, node_assignment, *node.br);
    }
    else if(node.is_leaf && node.contains_body)
    {
	node_assignment[nb_proc-1].push_back(&node);
	*(bodies_per_node+nb_proc-1) += 1;;
    }

}

std::ostream& operator<< (std::ostream & out, Quadtree const& qt)
{
    out << "root: (x,y) = (" <<  qt.root.center_x << "," <<  qt.root.center_y << ") and (w,h) = (" <<  qt.root.width << "," <<  qt.root.height << "), mass = " << qt.root.mass << ", nb_bodies = " << qt.root.nb_bodies  << " depth=" << qt.root.depth << "\n";
    if (!qt.root.is_leaf)
    {
	out << *qt.root.tr << "\n";
	out << *qt.root.tl << "\n";
	out << *qt.root.bl << "\n";
	out << *qt.root.br << "\n";
    }
    return out;
};

std::ostream& operator<< (std::ostream & out, Node const& node)
{    
    out << "node: (x,y) = (" <<  node.center_x << "," <<  node.center_y << ") and (w,h) = (" << node.width << ","  << node.height << "), mass = " <<  node.mass <<  ", nb_bodies = " << node.nb_bodies << " "  << node.is_leaf << " " << node.contains_body << ", depth=" << node.depth << "\n";
    if (!node.is_leaf)
    {
	out << *node.tr << "\n";
	out << *node.tl << "\n";
	out << *node.bl << "\n";
	out << *node.br << "\n";
    }
    return out;
};

std::ostream& operator<< (std::ostream & out, Body const& body)
{    
    out << "body: (x,y) = (" <<  body.pos_x << "," <<  body.pos_y << "), mass = " << body.mass << " and acc: (x,y) = " << body.acc_x << "," << body.acc_y << "\n";
    return out;
};



#include <Eigen/Dense>
#include <iostream>
#include "conex/debug_macros.h"
#include "example_graphs.h"
using Eigen::VectorXd;
using Eigen::MatrixXd;

using Iterate = VectorXd;
using std::vector;
using Edge = std::pair<int, int>; 

namespace conex {
struct Node {
  std::vector<int> incoming_edges;
  std::vector<int> outgoing_edges;
};

struct ProblemData {
  int num_node_variables;
  int num_flow_variables;
  std::vector<Node> topology;
  std::vector<Edge> edges;
};

class FlowConstraints {
 public:
  FlowConstraints(const ProblemData* data) : data_(data) { AssembleConstraintMatrix(); }
  VectorXd Evaluate(const VectorXd& x) const { 
    VectorXd y(data_->num_node_variables);
    y.setZero();
    int i = 0;
    for (const auto& node : data_->topology) {
      for (const auto& edge : node.incoming_edges) {
        y(i) += x(edge);
      }
      for (const auto& edge : node.outgoing_edges) {
        y(i) -= x(edge);
      }
      ++i;
    }
    return y;
  }
    
  VectorXd EvaluateTranspose(const VectorXd& x) const { 
    int i = 0;
    VectorXd y(data_->num_flow_variables);
    y.setZero();
    for (const auto& edge : data_->edges) {
      y(i) =  x(edge.first) - x(edge.second);
      ++i;
    }
    return y;
  }

 private:
  void AssembleConstraintMatrix() {
    constraint_matrix_.resize(data_->num_node_variables, data_->num_flow_variables);
    constraint_matrix_.setZero();
    int i = 0;
    for (const auto& node : data_->topology) {
      for (const auto& edge : node.incoming_edges) {
        constraint_matrix_(i, edge) = 1;
      }
      for (const auto& edge : node.outgoing_edges) {
        constraint_matrix_(i, edge) = -1;
      }
      ++i;
    }
  }
  MatrixXd constraint_matrix_;
  const ProblemData* data_;
};

class NodeHessians {
 public:
  NodeHessians(const ProblemData* data) : data_(data) {}
  void AssembleAndFactorHessians(const Iterate& iterate)  { 
    constraint_matrix_.setIdentity(data_->num_flow_variables, data_->num_flow_variables ); 
    constraint_matrix_(1, 1) = 3;
    constraint_matrix_(0, 0) = 3;
  
  }
  VectorXd EvaluateInverse(const VectorXd& x) const { 
    return constraint_matrix_ * x; 
  }

 private:
  MatrixXd constraint_matrix_;
  const ProblemData* data_;
};


int SolveForDualVariables(const ProblemData& problem_data,   const Iterate& iterate, 
                          int iteration_limit) {

  FlowConstraints F(&problem_data);
  NodeHessians H(&problem_data); H.AssembleAndFactorHessians(iterate);


  auto f = [F, H](const VectorXd& s) -> VectorXd { 
    return F.Evaluate(H.EvaluateInverse(F.EvaluateTranspose(s)));
  };

  int num_rows = problem_data.num_node_variables;
  int max_iter = 3;
  VectorXd b(num_rows); b.setConstant(1);
  /* plant a feasible solution */
  VectorXd  x(num_rows); x.setConstant(0); x(0) = 1;
  b = f(x);

  VectorXd s(num_rows); s.setZero();
  VectorXd r(num_rows); r = b;
  VectorXd p(num_rows); p = r;


  for (int i = 0; i < iteration_limit; i++) {
    double norm_sqr_r  = r.dot(r);
    double alpha = norm_sqr_r/p.dot(f(p));
    s +=  alpha * p;
    r -=  alpha * f(p);
    double beta = r.dot(r)/norm_sqr_r;
    p = r + beta * p;
    std::cout << "norm: " << r.norm() << std::endl;
    if (r.norm() < 1e-8) {
      std::cout << "\nTerminating. ";
      break;
    }
  }
  DUMP(b);
  DUMP(f(s));
  DUMP(s);
}


/* 
 *  a --> b -->  c  --> d
 */
vector<Edge> LineGraph(int number_of_nodes) {
  vector<Edge> edges(number_of_nodes - 1);
  for (int i = 0; i < number_of_nodes - 1; i++) {
    edges.at(i) = Edge(i,  i+1);
  }
  edges.emplace_back(0, 2);
  edges.emplace_back(1, 2);
  return edges;
}


int GetMaxNodeIndex(const vector<Edge>& edge_list) {
  int max = 0;
  for (const auto& edge : edge_list) {
    if (edge.first > max) {
      max = edge.first;
    }
    if (edge.second > max) {
      max = edge.second;
    }
  }
  return max + 1;
}

vector<Node> MakeNodeList(const vector<Edge>& edge_list, int num_nodes) {
  vector<Node> y(num_nodes);
  int i = 0;
  for (const auto& edge : edge_list) {
    y.at(edge.first).incoming_edges.push_back(i);
    y.at(edge.second).outgoing_edges.push_back(i);
    ++i;
  }
  return y;
}

int DoMain() {
  ProblemData problem_data;
  /*problem_data.edges = LineGraph(3);*/
  problem_data.edges = PaperExample();
  problem_data.num_node_variables = GetMaxNodeIndex(problem_data.edges);
  problem_data.num_flow_variables = problem_data.edges.size();
  problem_data.topology = MakeNodeList(problem_data.edges, problem_data.num_node_variables);
  Iterate iterate;
  SolveForDualVariables(problem_data, iterate, 10);
}
}

int main() {
  conex::DoMain();
}

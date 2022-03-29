#include <Eigen/Dense>
#include <iostream>
#include "conex/debug_macros.h"
using Eigen::VectorXd;
using Eigen::MatrixXd;

using Iterate = VectorXd;
using std::vector;


struct Node {
  std::vector<int> incoming_edges;
  std::vector<int> outgoing_edges;
};

struct ProblemData {
  int num_node_variables;
  int num_flow_variables;
  std::vector<Node> topology;
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
    return constraint_matrix_.transpose() * x; }
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
  void AssembleAndFactorHessians(const Iterate& iterate)  { constraint_matrix_.setIdentity(data_->num_flow_variables, data_->num_flow_variables ); }
  VectorXd EvaluateInverse(const VectorXd& x) const { return constraint_matrix_ * x; }

 private:
  MatrixXd constraint_matrix_;
  const ProblemData* data_;
};


int SolveForDualVariables(const ProblemData& problem_data,   const Iterate& iterate) {

  FlowConstraints F(&problem_data);
  NodeHessians H(&problem_data); H.AssembleAndFactorHessians(iterate);


  auto f = [F, H](const VectorXd& s) -> VectorXd { 
    return F.Evaluate(H.EvaluateInverse(F.EvaluateTranspose(s)));
  };

  int num_rows = problem_data.num_node_variables;
  int max_iter = 3;
  VectorXd b(num_rows); b.setConstant(1);
  b(0) = 2;
  VectorXd s(num_rows); s.setZero();
  VectorXd r(num_rows); r = b;
  VectorXd p(num_rows); p = r;

  for (int i = 0; i < max_iter; i++) {
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
}

vector<Node> LineGraph(int number_of_nodes) {
  vector<Node> nodes(number_of_nodes);
  for (int i = 0; i < number_of_nodes; i++) {
    nodes[i].incoming_edges.push_back(i);
    nodes[i].outgoing_edges.push_back(i+1);
  }
  return nodes;
}

int DoMain() {
  ProblemData problem_data;
  problem_data.num_node_variables = 3;
  problem_data.num_flow_variables = 4;
  problem_data.topology = LineGraph(problem_data.num_node_variables);
  Iterate iterate;
  SolveForDualVariables(problem_data, iterate);
}

int main() {
  DoMain();
}

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "conex/debug_macros.h"
#include "example_graphs.h"
using Eigen::VectorXd;
using Eigen::MatrixXd;

using Iterate = VectorXd;
using std::vector;
using Edge = std::pair<int, int>; 

namespace conex {

struct Iterate {
  VectorXd flow_variables;
  double inv_sqrt_mu;
};

struct DualSolverConfig {
  int iteration_limit;
};

struct Node {
  std::vector<int> incoming_edges;
  std::vector<int> outgoing_edges;
};

struct ProblemData {
  int num_node_variables;
  int num_edges;
  std::vector<Node> topology;
  std::vector<Edge> edges;
  VectorXd linear_edge_cost;
};

class FlowConstraints {
 public:
  FlowConstraints(const ProblemData* data) : data_(data) { 
    if (!data) {
      throw std::runtime_error("Received null pointer.");
    }
    MatrixXd constraint_matrix = AssembleConstraintMatrix(*data);
    MatrixXd matrix = constraint_matrix * constraint_matrix.transpose();
    Eigen::SparseMatrix<double> sparse_matrix = matrix.sparseView();
    Eigen::AMDOrdering<int> ordering;
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> perm;
    ordering(sparse_matrix, perm); 
    DUMP(perm * matrix * perm.transpose());
    DUMP(perm.transpose() * matrix * perm);
    DUMP(constraint_matrix);

  }
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
    VectorXd y(data_->num_edges);
    y.setZero();
    for (const auto& edge : data_->edges) {
      y(i) =  x(edge.first) - x(edge.second);
      ++i;
    }
    return y;
  }

  static MatrixXd AssembleConstraintMatrix(const ProblemData& data) {
    MatrixXd constraint_matrix(data.num_node_variables, data.num_edges);
    constraint_matrix.setZero();
    int i = 0;
    for (const auto& node : data.topology) {
      for (const auto& edge : node.incoming_edges) {
        constraint_matrix(i, edge) = 1;
      }
      for (const auto& edge : node.outgoing_edges) {
        constraint_matrix(i, edge) = -1;
      }
      ++i;
    }
    return constraint_matrix;
  }
 
 private:
  const ProblemData* data_;
};

class NodeHessians {
 public:
  NodeHessians(const ProblemData* data) : data_(data) {}
  void AssembleAndFactorHessians(const Iterate& iterate)  { 
    constraint_matrix_ = iterate.flow_variables.array().square().matrix().asDiagonal();
  }

  VectorXd EvaluateInverse(const VectorXd& x) const { 
    return constraint_matrix_ * x; 
  }

 private:
  MatrixXd constraint_matrix_;

  const ProblemData* data_;
};

VectorXd SolveForNewtonDirection(const ProblemData& problem_data,   
                          double inv_sqrt_mu, 
                          const Iterate& iterate, 
                          const DualSolverConfig& config) {

  FlowConstraints F(&problem_data);
  NodeHessians H(&problem_data); H.AssembleAndFactorHessians(iterate);

  auto f = [F, H](const VectorXd& s) -> VectorXd { 
    return F.Evaluate(H.EvaluateInverse(F.EvaluateTranspose(s)));
  };

  int num_rows = problem_data.num_node_variables;
  int max_iter = 3;
  VectorXd b(num_rows); b.setConstant(0);
  /* plant a feasible solution */
  if (0) {
    VectorXd  x(num_rows); x.setConstant(0); x(0) = 1;
    b = f(x);
  } else {
    b(0) = 1;
    b.tail(1)(0) = -1;
  }

  VectorXd d = VectorXd::Constant(problem_data.num_edges, 1);
  VectorXd residual = inv_sqrt_mu *
            (b  + F.Evaluate(H.EvaluateInverse( problem_data.linear_edge_cost )) ) - 2 * F.Evaluate(iterate.flow_variables);
  {
    VectorXd s(num_rows); s.setZero();
    VectorXd r(num_rows); r = residual;
    VectorXd p(num_rows); p = r;

    for (int i = 0; i < config.iteration_limit; i++) {
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
    d -= iterate.flow_variables.cwiseProduct(inv_sqrt_mu * problem_data.linear_edge_cost -F.EvaluateTranspose(s));
  }
  return d;
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
  //problem_data.edges = LineGraph(3);
  problem_data.edges = PaperExample();
  problem_data.num_node_variables = GetMaxNodeIndex(problem_data.edges);
  problem_data.num_edges = problem_data.edges.size();
  problem_data.topology = MakeNodeList(problem_data.edges, problem_data.num_node_variables);
  Iterate iterate;
  iterate.flow_variables = VectorXd::Constant(problem_data.num_edges, 1);
  DualSolverConfig config; config.iteration_limit = problem_data.num_edges;
  problem_data.linear_edge_cost = VectorXd::Constant(problem_data.num_edges, 1);

  double inv_sqrt_mu = 100;
  for (int i = 0; i < 20; i++) {
    VectorXd d = SolveForNewtonDirection(problem_data, inv_sqrt_mu, iterate, config);

    double norminf = (d).array().abs().maxCoeff();
    if (norminf > 1) {
      d /= norminf;
    }
    iterate.flow_variables = iterate.flow_variables.cwiseProduct(d.array().exp().matrix());
    DUMP(norminf);
  }


}

void TestFlowConstraints() {

  ProblemData problem_data;
  problem_data.edges = PaperExample();
  problem_data.num_node_variables = GetMaxNodeIndex(problem_data.edges);
  problem_data.num_edges = problem_data.edges.size();
  problem_data.topology = MakeNodeList(problem_data.edges, problem_data.num_node_variables);

  FlowConstraints F(&problem_data);
  MatrixXd Fmat = FlowConstraints::AssembleConstraintMatrix(problem_data);
  double eps = 1e-15;
  for (int i = 0; i < 10; i++) {
    MatrixXd test(Fmat.cols(), 1); test.setRandom();
    double error = (Fmat * test - F.Evaluate(test) ).norm();
    if (error > eps) {
      throw std::runtime_error("Test failed");
    }
  }
  for (int i = 0; i < 10; i++) {
    MatrixXd test(Fmat.rows(), 1); test.setRandom();
    double error = (Fmat.transpose() * test - F.EvaluateTranspose(test) ).norm();
    if (error > eps) {
      throw std::runtime_error("Test failed");
    }
  }
}

}

/*
 *  min. c(f, x) 
 *
 *    Fx = b
 *    x >= 0
 *
 *
 *  Newton equations:
 *
 *  F e^v(1 + d) = b 
 *  e^{-v}(1 - d) = c - F'y
 *  
 * 
 * F Q(v) F' y = F e^v - (F' Q(v) c + b)
 *
 *
 * W(v)'
 */



int main() {
  conex::DoMain();
  conex::TestFlowConstraints();
}

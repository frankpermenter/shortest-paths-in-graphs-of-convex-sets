#include <vector>
#include "common.h"
#include <map>
using std::vector;
using std::string;
using std::map;
std::vector<Edge> PaperExample() {

/*edges = {
/*    's': ('p0', 'p1', 'p2'),
/*    'p0': ('e1',),
/*    'p1': ('p2', 'e0', 'e1'),
/*    'p2': ('p1', 'p3', 'e0'),
/*    'p3': ('t', 'p2', 'p4', 'e1'),
/*    'p4': ('t', 'e0'),
/*    'e0': ('p3', 'p4', 'e1'),
/*    'e1': ('t', 'p4', 'e0')
/*}
 */
//for u, vs in edges.items():
//    for v in vs:
//        G.add_edge(u, v, l)
//
//

map<string, int> variable_to_index;
map<string, vector<string>>  edges;

variable_to_index["s"] =  0;
variable_to_index["p0"] = 1;
variable_to_index["p1"] = 2;
variable_to_index["p2"] = 3;
variable_to_index["p3"] = 4;
variable_to_index["p4"] = 5;
variable_to_index["e0"] = 6;
variable_to_index["e1"] = 7;
variable_to_index["t"] = 8;
                               

edges["s"] =  {"p0", "p1", "p2"};
edges["p0"] = {"e1"}                ;
edges["p1"] = {"p2", "e0", "e1"}     ;
edges["p2"] = {"p1", "p3", "e0"}     ;
edges["p3"] = {"t", "p2", "p4", "e1"};
edges["p4"] = {"t", "e0"}            ;
edges["e0"] = {"p3", "p4", "e1"}     ;
edges["e1"] = {"t", "p4", "e0"};


std::vector<Edge> edge_list;
for (const auto& node_edge_pair : edges) {
  int outgoing_index = variable_to_index[node_edge_pair.first];
  for (const auto& incoming_edge : node_edge_pair.second) {
    int incoming_index = variable_to_index[incoming_edge];
    edge_list.push_back({outgoing_index, incoming_index});
  }
}


return edge_list;

}

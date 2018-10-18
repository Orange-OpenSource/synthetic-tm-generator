#!/usr/bin/env python
"""
Copyright 2018 Orange

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from __future__ import print_function
import yaml
import numpy as np
import sys
import random
import argparse
from subprocess import call
import tempfile
import re
import os

from pyomo.environ import *
from pyomo.opt import SolverFactory

import networkx as nx
from matplotlib.cm import autumn_r, hot_r
import matplotlib.pyplot as plt


"""Generates nb_nodes * (nb_nodes - 1) flow rates according to the Gravity model.
Args:
  nb_nodes: Number of nodes in the network.
  total_traffic: Total traffic going through the network, in Mbps.
Returns:
  The flow rates as a matrix indexed by node numbers. Cells in the diagonal are null.
"""
def generate_tm(nb_nodes, total_traffic):
    mean = 0.086

    # rate = 1 / mean
    tin = np.random.exponential(mean, nb_nodes)
    tout = np.random.exponential(mean, nb_nodes)
    tm = np.zeros((nb_nodes, nb_nodes))
    sum_tin = np.sum(tin)
    sum_tout = np.sum(tout)
    for i in range(nb_nodes):
        for j in range(nb_nodes):
            if i != j:
                tm[i][j] = total_traffic * tin[i] / sum_tin * tout[j] / sum_tout

    return tm


"""Defines the variables, the constraints, and the objective for use by Pyomo.
Args:
  nodes: List of node IDs.
  links: A dictionary of link capacities (in Mbps) indexed by link's edges.
  od_pairs: List of OD pair (as tuples).
  flow_rates: List of flow rates, in Mbps.
Returns:
  The Pyomo model.
"""
def define_model(nodes, links, od_pairs, flow_rates):
    L = range(len(flow_rates))
    model = ConcreteModel()

    model.mapping = Var(od_pairs, L, within=Binary)
    model.pair_rates = Var(od_pairs, within=NonNegativeReals)
    model.link_rates = Var(links.keys(), od_pairs, within=NonNegativeReals)

    model.z = Var(within=NonNegativeReals)

    def max_link_load(model):
        return model.z
    model.obj = Objective(rule=max_link_load, sense=minimize)
    def single_rate_per_pair(model, k1, k2):
        return sum(model.mapping[k1, k2, l] for l in L) == 1
    model.nine = Constraint(od_pairs, rule=single_rate_per_pair)
    def single_pair_per_rate(model, l):
        return sum(model.mapping[k1, k2, l] for (k1, k2) in od_pairs) == 1
    model.ten = Constraint(L, rule=single_pair_per_rate)
    def compute_pair_rates(model, k1, k2):
        return model.pair_rates[k1, k2] == sum(model.mapping[k1, k2, l] * flow_rates[l] for l in L)
    model.eleven = Constraint(od_pairs, rule=compute_pair_rates)
    def flow_conservation(model, u_, k1, k2):
        if u_ == k1:
            return sum(model.link_rates[v, u, k1, k2] for (v, u) in links.keys() if u_ == u) - sum(model.link_rates[u, v, k1, k2] for (u, v) in links.keys() if u_ == u) == -model.pair_rates[k1, k2]
        if u_ == k2:
            return sum(model.link_rates[v, u, k1, k2] for (v, u) in links.keys() if u_ == u) - sum(model.link_rates[u, v, k1, k2] for (u, v) in links.keys() if u_ == u) == model.pair_rates[k1, k2]
        return sum(model.link_rates[v, u, k1, k2] for (v, u) in links.keys() if u_ == u) - sum(model.link_rates[u, v, k1, k2] for (u, v) in links.keys() if u_ == u) == 0
    model.twelve = Constraint(nodes, od_pairs, rule=flow_conservation)
    def compute_max_link_load(model, u, v):
        return model.z >= (sum(model.link_rates[u, v, k1, k2] for (k1, k2) in od_pairs) / links[u, v])
    model.thirteen = Constraint(links.keys(), rule=compute_max_link_load)

    return model


"""Assigns flow rates to the given network using the ILP method.

Uses GLPK as the solver and Pyomo as the Python interface.

Args:
  data: The dictionary containing the topology and the node capacities. 
  tm: The matrix of generated flow rates.
  mipgap: The mipgap argument for the GLPK solver.
Returns:
  0
"""
def assign_flow_rates_ilp(data, tm, mipgap):
    nodes = [node['id'] for node in data['nodes']]
    links = {}
    for link in data['links']:
        links[(link['source'], link['destination'])] = link['capacity']
    od_pairs = []
    flow_rates = []
    for i in range(1, len(nodes) + 1):
        for j in range(1, len(nodes) + 1):
            if i != j:
                flow_rates.append(int(round(tm[i - 1, j - 1])))
                od_pairs.append((i, j))

    model = define_model(nodes, links, od_pairs, flow_rates)

    solver = SolverFactory('glpk')
    solver.set_options("mipgap=%f" % mipgap)
    results = solver.solve(model)
    print(results)
    if str(results['Solver'][0]['Termination condition']) == 'infeasible':
        return None
    model.solutions.load_from(results)

    for link in data['links']:
        for (s, d) in od_pairs:
            link['legit_load'] += model.link_rates[link['source'], link['destination'], s, d].value

    return 0


"""Assigns flow rates to the given network using the ILP method.

Uses GLPK as the solver but without relying on Pyomo as the Python interface.
Instead, it generates the model file with all the parameters and directly calls glpsol.

Args:
  data: The dictionary containing the topology and the node capacities.
  tm: The matrix of generated flow rates.
  mipgap: The mipgap argument for the GLPK solver.
Returns:
  The objective value, the maximum link load.
"""
def assign_flow_rates_ilp_glpk(data, tm, mipgap):
    # Generate model file from template GMPL program:
    nb_od_pairs = 1
    flow_rates_str = ""
    od_pairs_str = ""
    for i in range(1, len(nodes) + 1):
        for j in range(1, len(nodes) + 1):
            if i != j:
                flow_rates_str += "\n  %d %f" % (nb_od_pairs, tm[i - 1, j - 1])
                nb_od_pairs += 1
                od_pairs_str += "\n  %d %d" % (i, j)
    nb_od_pairs -= 1

    links_str = ""
    link_capacities_str = ""
    nb_links = 0
    for link in data['links']:
        links_str += "\n  %d %d" % (link['source'], link['destination'])
        link_capacities_str += "\n  %d %d %d" % (link['source'], link['destination'], link['capacity'])
        nb_links += 1

    with open("template.mod", 'r') as fh:
        mip_program = fh.read()
        mip_program = mip_program.replace("NB_NODES", str(len(data['nodes'])))
        mip_program = mip_program.replace("NB_OD_PAIRS", str(nb_od_pairs))
        mip_program = mip_program.replace("OD_PAIRS", od_pairs_str)
        mip_program = mip_program.replace("LINKS", links_str)
        mip_program = mip_program.replace("FLOW_RATES", flow_rates_str)
        mip_program = mip_program.replace("LINK_CAPABILITIES", link_capacities_str)
    _, tmp_model_file = tempfile.mkstemp()
    with open(tmp_model_file, 'w') as fh:
        fh.write(mip_program)

    # Run GLPK
    _, tmp_output_file = tempfile.mkstemp()
    call(["glpsol", "--mipgap", str(mipgap), "--model", tmp_model_file, "-o", tmp_output_file])

    # Retrieve variables' values:
    nb_constraints = 1 + 3 * nb_od_pairs + nb_od_pairs * len(data['nodes']) + nb_links
    nb_mapping_vars = nb_od_pairs * nb_od_pairs
    nb_pair_rate_vars = nb_od_pairs
    link_rates = {}
    objective = None
    with open(tmp_output_file, 'r') as fh:
        regexp = re.compile(r'^\s+\d+\s+link_rates\[(\d+),(\d+),\d+,\d+\]$')
        val_line = False
        for line in fh:
            if val_line:
                if (src, dst) in link_rates:
                    link_rates[src, dst] += float(line.split()[0])
                else:
                    link_rates[src, dst] = float(line.split()[0])
                val_line = False
            else:
                m = regexp.match(line)
                if m:
                    src = int(m.group(1))
                    dst = int(m.group(2))
                    val_line = True
                elif "Objective:  max_load" in line:
                    objective = float(line.split()[3])
    if nb_links != len(link_rates):
        print("Error: The number of link_rates variables retrieved from the glpsol output doesn't match the number of links. %d %d", (nb_links, len(link_rates)), file=sys.stderr)
        sys.exit(1)

    for link in data['links']:
        link['legit_load'] = link_rates[link['source'], link['destination']]

    os.remove(tmp_model_file)
    os.remove(tmp_output_file)
    return objective


"""Routes flows between OD pairs on the network according to all_paths.
Args:
  links: Dictionary of links, with the legitimate load value indexed under 'legit_load', indexed by link's edges.
  all_paths: Dictionary of tuples' lists, indexed by OD pair's edges. Tuples contain a path as a list of nodes to traverse and a weight for that path. 
  od_pairs: List of OD pairs (as tuples).
  flow_rates: List of flow rates. Flow rates i will be routed on OD pair i.
"""
def route_flows_multipaths(links, all_paths, od_pairs, flow_rates):
    for i in xrange(len(od_pairs)):
        (src, dst) = od_pairs[i]
        paths = all_paths[src, dst]
        for (path, weight) in paths:
            flow_rate = flow_rates[i] * weight
            for j in range(len(path) - 1):
                links[path[j], path[j + 1]]['legit_load'] += flow_rate


"""Computes the capacity of a given path.
Args:
  links: Dictionary of links, with the legitimate load value indexed under 'capacity', indexed by link's edges.
  path: Path, as a list of nodes.
Returns:
  The lowest link capacity on the path.
"""
def path_capacity(links, path):
    min_cap = float("inf")
    for i in range(len(path) - 1):
        min_cap = min(min_cap, links[path[i], path[i + 1]]['capacity'])
    return min_cap


"""Assigns flow rates to the given network using the heuristic method.



Args:
  data: The dictionary containing the topology and the node capacities.
  tm: The matrix of generated flow rates.
  mipgap: Not used. For compatibility with other methods' functions. 
Returns:
  The maximum link load.
"""
def assign_flow_rates_heuristic(data, tm, mipgap):
    flow_rates = []
    od_pairs = []
    for i in range(1, len(nodes) + 1):
        for j in range(1, len(nodes) + 1):
            if i != j:
                flow_rates.append(tm[i - 1, j - 1])
                od_pairs.append((i, j))

    # Compute multipath routing between all node pairs:
    links = {}
    for link in data['links']:
        links[link['source'], link['destination']] = {'legit_load': 0, 'capacity': link['capacity']}
    G = nx.DiGraph([(link['source'], link['destination']) for link in data['links']])
    all_paths = {}
    for (src, dst) in od_pairs:
        all_paths[src, dst] = []
        try:
            paths = [p for p in nx.all_shortest_paths(G, source=src, target=dst)]
        except nx.exception.NetworkXNoPath:
            print("Error: No path found between %d and %d." % (src, dst), file=sys.stderr)
            sys.exit(1)
        max_flow = 0
        path_capacities = []
        for path in paths:
            capacity = path_capacity(links, path)
            path_capacities.append(capacity)
            max_flow += capacity
        for i in range(len(paths)):
            weight = path_capacities[i] * 1.0 / max_flow
            all_paths[src, dst].append((paths[i], weight))

    # Compute sort information on OD pairs:
    node_infos = {}
    for node in data['nodes']:
        node_infos[node['id']] = {'fanout': 0, 'fanin': 0, 'connectivity': 0, 'nb_paths': 0}

    for link in data['links']:
        node_infos[link['destination']]['connectivity'] += 1
        node_infos[link['destination']]['fanin'] += link['capacity']
        node_infos[link['source']]['connectivity'] += 1
        node_infos[link['source']]['fanout'] += link['capacity']

    for (src, dst) in all_paths:
        for (path, _) in all_paths[src, dst]:
            for node in path[1:-1]:
                node_infos[node]['nb_paths'] += 1

    od_pair_infos = {}
    for (src, dst) in od_pairs:
        src_infos = node_infos[src]
        dst_infos = node_infos[dst]
        m1 = min(src_infos['fanout'], dst_infos['fanin'])
        m2 = min(src_infos['connectivity'], dst_infos['connectivity'])
        if src_infos['nb_paths'] == dst_infos['nb_paths'] == 0:
            m3 = float("Inf")
        else:
            m3 = 1.0 / max(src_infos['nb_paths'], dst_infos['nb_paths'])
        od_pair_infos[(src, dst)] = {'m1': m1, 'm2': m2, 'm3': m3}

    # Sort OD pairs:
    def make_comparator(od_pair_infos):
        def compare(od1, od2):
            pair1_infos = od_pair_infos[od1]
            pair2_infos = od_pair_infos[od2]
            if pair1_infos['m1'] == pair2_infos['m1']:
                if pair1_infos['m2'] == pair2_infos['m2']:
                    if pair1_infos['m3'] == pair2_infos['m3']:
                        return 0
                    elif pair1_infos['m3'] > pair2_infos['m3']:
                        return 1
                    else:
                        return -1
                elif pair1_infos['m2'] > pair2_infos['m2']:
                    return 1
                else:
                    return -1
            elif pair1_infos['m1'] > pair2_infos['m1']:
                return 1
            else:
                return -1
        return compare
    flow_rates.sort(reverse=True)
    od_pairs = sorted(od_pairs, cmp=make_comparator(od_pair_infos), reverse=True)

    # Route flows between OD pairs:
    route_flows_multipaths(links, all_paths, od_pairs, flow_rates)


    # Write link loads to YAML and compute objective value:
    objective = 0
    for link in data['links']:
        link['legit_load'] = links[link['source'], link['destination']]['legit_load']
        max_link_load = link['legit_load'] / link['capacity']
        if max_link_load > objective:
            objective = max_link_load

    return objective


"""Displays the graph of the network with legitimate loads on edges, in a new window.
Args:
  data: The dictionary containing the topology and the node capacities.
"""
def display_graph(data):
    links = [(link['source'], link['destination']) for link in data['links']]
    edge_labels = {}
    for link in data['links']:
        if link['legit_load'] > 0:
            edge_labels[(link['source'], link['destination'])] = link['legit_load']
    G = nx.DiGraph(links)
    pos = nx.spring_layout(G)
    nx.draw(G, pos=pos, with_labels=True, arrows=False)
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)
    plt.show()


"""Computes the mean load factor (load / capacity) of links in the network.
Args:
  data: The dictionary containing the topology and the node capacities.
Returns:
  The mean link load factor.
"""
def compute_mean_link_load(data):
    tot_link_load = 0
    nb_links = 0
    for link in data['links']:
        tot_link_load += link['legit_load'] / link['capacity']
        nb_links += 1
    return tot_link_load / nb_links


"""Scales down the legitimate loads of links in the network by factor.
Args:
  data: The dictionary containing the topology and the node capacities.
  factor: The value link loads are divided by.
"""
def scale_down_tm(data, factor):
    for link in data['links']:
        link['legit_load'] /= factor


"""Generates node and link items for the output YAML file.

For each link (u, v), a back link (v, u) of equal capacity is created.
Link capacities are all adapted according to the max_link_capacity.

Args:
  data: The dictionary containing the output information.
  template: The template information, directly from the YAML file (dictionary).
  max_node_class: The maximum node class to include in the output network. Used to restrict the network's size.
  max_link_capacity: The capacity of the largest links.
Returns:
  The list of selected nodes, according to the max_node_class argument.
"""
def generate_topology(data, template, max_node_class, max_link_capacity):
    classes = {}
    for class_ in template['classes']:
        classes[class_['name']] = class_

    # Collects neighbors for each node:
    all_neighbors = {}
    for node in template['nodes']:
        neighbors = []
        all_neighbors[node['id']] = neighbors
        for link in template['links']:
            if link['destination'] == node['id']:
                neighbors.append((link['source'], link['capacity']))
            elif link['source'] == node['id']:
                neighbors.append((link['destination'], link['capacity']))
        neighbors.sort(key=lambda tup: tup[0])

    # Collects links:
    links = {}
    for link in template['links']:
        if link['source'] < link['destination']:
            links[(link['source'], link['destination'])] = link['capacity']
        else:
            links[(link['destination'], link['source'])] = link['capacity']

    # Selects the nodes according to the max. node class wanted:
    nodes = []
    i = 0
    for node in template['nodes']:
        if node['class'] <= max_node_class:
            class_ = classes[node['class']]
            new_node = {'id': node['id'], 'cpus': class_['cpus'], 'memory': class_['memory']}
            data['nodes'].append(new_node)
            nodes.append(node['id'])
        else:
            # Removed node, need to bridge the neighbors:
            nb_neighbors = len(all_neighbors[node['id']])
            if nb_neighbors >= 2 and nb_neighbors <= 3:
                for i in range(nb_neighbors):
                    neighbor1 = all_neighbors[node['id']][i]
                    neighbor2 = all_neighbors[node['id']][(i + 1) % nb_neighbors]
                    if neighbor1[0] == neighbor2[0]:
                    # Only allow edges between different nodes.
                        continue
                    # Retrieves the max capacity between the already existing link, if any, and the new:
                    capacity = max(neighbor1[1], neighbor2[1])
                    if (neighbor1[0], neighbor2[0]) in links:
                        link2_capacity = links[(neighbor1[0], neighbor2[0])]
                        if link2_capacity > capacity:
                            capacity = link2_capacity
                        else:
                            continue
                    link = {'source': neighbor1[0], 'destination': neighbor2[0], 'capacity': capacity}
                    template['links'].insert(0, link)
                    # Removes the current node from neighbor lists:
                    all_neighbors[link['source']] = [(u, cap) for (u, cap) in all_neighbors[link['source']] if u != node['id']]
                    all_neighbors[link['destination']] = [(u, cap) for (u, cap) in all_neighbors[link['destination']] if u != node['id']]
                    # Adds the new neighbors:
                    all_neighbors[link['source']].append((link['destination'], link['capacity']))
                    all_neighbors[link['destination']].append((link['source'], link['capacity']))
                    if nb_neighbors == 2:
                        # If we continue we'll add a back-edge between the two neighbors. 
                        break
        i += 1

    # Selects the links according to the remaining nodes:
    cur_max_link_capacity = max(template['links'], key=lambda link: link['capacity'])['capacity']
    link_size_factor = max_link_capacity / cur_max_link_capacity
    for link in template['links']:
        if link['source'] in nodes and link['destination'] in nodes:
            already_added = sum([l['source'] == link['source'] and l['destination'] == link['destination'] for l in data['links']])
            if already_added == 0:
                link['capacity'] = link['capacity'] * link_size_factor
                data['links'].append(link)
                back_link = dict(link)
                back_link['source'] = link['destination']
                back_link['destination'] = link['source']
                data['links'].append(back_link)

    return nodes


"""Generates the legitimate link loads according to the Gravity model.

Displays a few info messages on stdout. 

Args:
  data: The dictionary containing the topology and the node capacities.
  nodes: The list of node IDs.
  mipgap: The mipgap argument for the GLPK solver.
  total_traffic: 
  target_mean_link_load: 
  flow_assign_method: The method to use to map generated flow rates into the network. One of 'heuristic', 'ilp', or 'ilp-glpk'.
"""
def generate_link_loads(data, nodes, mipgap, total_traffic, target_mean_link_load, flow_assign_method):
    methods = {'heuristic': assign_flow_rates_heuristic, 'ilp': assign_flow_rates_ilp, 'ilp-glpk': assign_flow_rates_ilp_glpk}
    tm = generate_tm(len(nodes), total_traffic)
    objective = methods[flow_assign_method](data, tm, mipgap)
    if objective > 1:
        print("Scaling down TM by %f to reach feasible routing." % objective, file=sys.stderr)
        scale_down_tm(data, objective)
    mean_link_load = compute_mean_link_load(data)
    factor =  mean_link_load / target_mean_link_load
    if factor < 1:
        print("Mean link load is at %f." % mean_link_load, file=sys.stderr)
    else:
        print("Scaling down TM by %f to reach %f mean link load." % (factor, args.mean_link_load), file=sys.stderr)
        scale_down_tm(data, factor)


"""Generates attacks between the given attackers and targets.

Attack loads follow an exponential distribution.
Attackers and targets are selected randomly. 

Args:
  data: The dictionary containing the topology and the node capacities.
  nodes: The list of node IDs.
  nb_attackers: The number attackers. 0 for all nodes.
  nb_targets: The number of targets.
  mean_attack_load: Mean attack load, in Mbps.
"""
def generate_attacks(data, nodes, nb_attackers, nb_targets, mean_attack_load):
    random.shuffle(nodes)
    if nb_attackers == 0:
        attackers = nodes
    else:
        attackers = nodes[0:nb_attackers]
    targets = nodes[0:nb_targets]
    for attacker in attackers:
        target = targets[random.randint(0, nb_targets - 1)]
        load = round(np.random.exponential(mean_attack_load))
        attack = {'source': attacker, 'destination': target, 'load': load}
        data['attacks'].append(attack)


"""Rounds values for the legitimate loads on links.
Args:
  data: The dictionary containing the topology and the node capacities.
"""
def round_link_loads(data):
    for link in data['links']:
        link['legit_load'] = round(link['legit_load'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a set of inputs for the VNF placement problem.')
    parser.add_argument('template_file',
                        help='Path to the template file to use.')
    parser.add_argument('output_file',
                        help='Path to the output file to use. - for stdout.')
    parser.add_argument('--method-flow-assign', choices=['ilp', 'ilp-glpk', 'heuristic'], default='heuristic',
                        help='Method for the assignment of flow rates to the network. ilp an dilp-glpk use the'
                        ' same solver, except ilp-glpk is tailored for GLPK and has a lower memory footprint. Defaults to %(default)s.')
    parser.add_argument('--max-node-class', type=int, default=3,
                        help='The maximum class of nodes to select (1, 2, or 3). Defaults to %(default)s to select all nodes from the template.')
    parser.add_argument('--max-link-capacity', type=int, default=100000,
                        help='The capacity assigned to the largest links, in Mbps. Defaults to 100 Gbps.')
    parser.add_argument('--total-traffic', type=int, default=1000000,
                        help='Total traffic going through the network, in Mbps. Used by the gravity model. Defaults to 1 Tbps.')
    parser.add_argument('nb_attackers', type=int,
                        help='Number of nodes from which the attacking flows arrive. 0 for all nodes.')
    parser.add_argument('nb_targets', type=int, help='Number of targets for the DDoS attacks.')
    parser.add_argument('--mean-attack-load', type=int,
                        default=10000, help='Mean attack load, in Mbps. Used to generate the attack loads. Defaults to 10 Gbps.')
    parser.add_argument('--mipgap', type=float, default=0.2,
                        help='mipgap parameter for GLPK. Defaults to %(default)s.')
    parser.add_argument('--mean-link-load', type=float, default=0.5, help='Mean link load. Defaults to %(default)s.')
    parser.add_argument('--display-network', action='store_true', help='Display the network graph with legitimate link loads in a new window.')
    args = parser.parse_args()

    with open(args.template_file, 'r') as fh:
        try:
            template = yaml.load(fh)
        except yaml.YAMLError as e:
            print(e, file=sys.stderr)
            sys.exit(1)

    data = {'nodes': [], 'links': [], 'attacks': [], 'vnfs': template['vnfs']}
    nodes = generate_topology(data, template, args.max_node_class, args.max_link_capacity)
    generate_link_loads(data, nodes, args.mipgap, args.total_traffic, args.mean_link_load, args.method_flow_assign)
    generate_attacks(data, nodes, args.nb_attackers, args.nb_targets, args.mean_attack_load)
    round_link_loads(data)

    if args.output_file == '-':
        print(yaml.dump(data, default_flow_style=False))
    else:
        with open(args.output_file, 'w') as fh:
            yaml.dump(data, fh, default_flow_style=False)

    if args.display_network:
        display_graph(data)

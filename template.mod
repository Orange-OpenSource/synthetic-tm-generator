# Copyright 2018 Orange
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#  http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Variables
############
set od_pairs, dimen 2;
set L := 1..NB_OD_PAIRS;
set nodes := 1..NB_NODES;
set links, dimen 2;

var z >= 0;
var mapping{od_pairs, L}, binary;
var pair_rates{od_pairs}, >= 0;
var link_rates{links, od_pairs}, >= 0;

param flow_rates {L};
param link_capacities {links};

# Objective
############
minimize max_load: z;

# Constraints
##############
s.t. single_rate_per_pair {(k1, k2) in od_pairs}:
  sum {l in L} mapping[k1, k2, l] = 1;

s.t. single_pair_per_rate {l in L}:
  sum {(k1, k2) in od_pairs} mapping[k1, k2, l] = 1;

s.t. compute_pair_rates {(k1, k2) in od_pairs}:
  pair_rates[k1, k2] == sum {l in L} mapping[k1, k2, l] * flow_rates[l];

s.t. flow_conservation1 {u_ in nodes, (k1, k2) in od_pairs : u_ != k1 and u_ != k2}:
  sum {(v, u) in links : u_ = u} link_rates[v, u, k1, k2] - sum {(u, v) in links : u_ = u} link_rates[u, v, k1, k2] = 0;

s.t. flow_conservation2 {u_ in nodes, (k1, k2) in od_pairs : u_ = k1}:
  sum {(v, u) in links : u_ = u} link_rates[v, u, k1, k2] - sum {(u, v) in links : u_ = u} link_rates[u, v, k1, k2] = - pair_rates[k1, k2];

s.t. flow_conservation3 {u_ in nodes, (k1, k2) in od_pairs : u_ = k2}:
  sum {(v, u) in links : u_ = u} link_rates[v, u, k1, k2] - sum {(u, v) in links : u_ = u} link_rates[u, v, k1, k2] = pair_rates[k1, k2];

s.t. compute_max_link_load {(u, v) in links}:
  z >= sum {(k1, k2) in od_pairs} link_rates[u, v, k1, k2] / link_capacities[u, v];


data;

set od_pairs := OD_PAIRS;
set links := LINKS;

param flow_rates := FLOW_RATES;
param link_capacities := LINK_CAPABILITIES;

end;

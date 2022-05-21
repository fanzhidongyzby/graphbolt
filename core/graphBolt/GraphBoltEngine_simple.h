// Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef GRAPHBOLT_ENGINE_SIMPLE_H
#define GRAPHBOLT_ENGINE_SIMPLE_H

#include "GraphBoltEngine.h"

// ======================================================================
// GRAPHBOLTENGINESIMPLE
// ======================================================================
template <class vertex, class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
class GraphBoltEngineSimple
    : public GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                             GlobalInfoType> {
public:
  GraphBoltEngineSimple(graph<vertex> &_my_graph, int _max_iter,
                        GlobalInfoType &_static_data, bool _use_lock,
                        commandLine _config)
      : GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>(_my_graph, _max_iter, _static_data,
                                        _use_lock, _config) {
    use_source_contribution = true;
  }

  // ======================================================================
  // TEMPORARY STRUCTURES USED BY THE SIMPLE ENGINE
  // ======================================================================
  void createTemporaryStructures() {
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::createTemporaryStructures();
  }
  void resizeTemporaryStructures() {
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::resizeTemporaryStructures();
    initTemporaryStructures(n_old, n);
  }
  void freeTemporaryStructures() {
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::freeTemporaryStructures();
  }
  void initTemporaryStructures() { initTemporaryStructures(0, n); }
  void initTemporaryStructures(long start_index, long end_index) {
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::initTemporaryStructures(start_index,
                                                             end_index);
  }
  // ======================================================================
  // TRADITIONAL INCREMENTAL COMPUTATION
  // ======================================================================
  // TODO : Currently, max_iterations = history_iterations.
  // Need to implement computation without history.
  int traditionalIncrementalComputation(int start_iteration) {
    timer iteration_timer, phase_timer;
    double misc_time, copy_time, phase_time, iteration_time;

    vertexSubset frontier_curr_vs(n, frontier_curr);
    bool use_delta = true;
    int iter = start_iteration;

    if (frontier_curr_vs.numNonzeros() == 0) {
      converged_iteration = start_iteration;

    } else {
      for (iter = start_iteration; iter < max_iterations; iter++) {
        // initialize timers
        {
          phase_timer.start();
          misc_time = 0;
          copy_time = 0;
        }

        // ========== COPY - Prepare curr iteration ==========
        if (iter > 0) {
          // Copy the aggregate and actual value from iter-1 to iter
          parallel_for(uintV v = 0; v < n; v++) {
            vertex_values[iter][v] = vertex_values[iter - 1][v];
            aggregation_values[iter][v] = aggregation_values[iter - 1][v];
            delta[v] = aggregationValueIdentity<AggregationValueType>();
          }
        }
        use_delta = shouldUseDelta(iter);

        // ========== MISC - count active edges for AE ==========
        phase_time = phase_timer.next();
        adaptive_executor.updateCopyTime(iter, phase_time);
        adaptive_executor.updateEdgesProcessed(iter, my_graph,
                                               frontier_curr_vs);
        misc_time = phase_timer.next();
        adaptive_executor.updateMiscTime(iter, phase_timer.next());

        // ========== EDGE COMPUTATION ==========
        if ((use_source_contribution) && (iter == 1)) {
          // Compute source contribution for first iteration
          parallel_for(uintV u = 0; u < n; u++) {
            if (frontier_curr[u]) {
              // compute source change in contribution
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  u, source_change_in_contribution[u],
                  vertexValueIdentity<VertexValueType>(),
                  vertex_values[iter - 1][u], global_info);
            }
          }
        }

        parallel_for(uintV u = 0; u < n; u++) {
          if (frontier_curr[u]) {
            // check for propagate and retract for the vertices.
            intE outDegree = my_graph.V[u].getOutDegree();
            granular_for(j, 0, outDegree, (outDegree > 1024), {
              uintV v = my_graph.V[u].getOutNeighbor(j);
              AggregationValueType contrib_change =
                  use_source_contribution
                      ? source_change_in_contribution[u]
                      : aggregationValueIdentity<AggregationValueType>();
#ifdef EDGEDATA
              EdgeData *edge_data = my_graph.V[u].getOutEdgeData(j);
#else
              EdgeData *edge_data = &emptyEdgeData;
#endif
              bool ret =
                  edgeFunction(u, v, *edge_data, vertex_values[iter - 1][u],
                               contrib_change, global_info);
              if (ret) {
                if (use_lock) {
                  vertex_locks[v].writeLock();
                  addToAggregation(contrib_change, delta[v], global_info);
                  vertex_locks[v].unlock();
                } else {
                  addToAggregationAtomic(contrib_change, delta[v], global_info);
                }
                if (!frontier_next[v])
                  frontier_next[v] = 1;
              }
            });
          }
        }

        phase_time = phase_timer.next();
        adaptive_executor.updateEdgeMapTime(iter, phase_time);

        // ========== VERTEX COMPUTATION ==========
        parallel_for(uintV v = 0; v < n; v++) {
          // Reset frontier for next iteration
          frontier_curr[v] = 0;
          // Process all vertices affected by EdgeMap
          if (frontier_next[v] ||
              forceComputeVertexForIteration(v, iter, global_info)) {

            frontier_next[v] = 0;
            // Update aggregation value and reset change received[v] (i.e.
            // delta[v])
            addToAggregation(delta[v], aggregation_values[iter][v],
                             global_info);
            delta[v] = aggregationValueIdentity<AggregationValueType>();

            // Calculate new_value based on the updated aggregation value
            VertexValueType new_value;
            computeFunction(v, aggregation_values[iter][v],
                            vertex_values[iter - 1][v], new_value, global_info);

            // Check if change is significant
            if (notDelZero(new_value, vertex_values[iter - 1][v], global_info)) {
              // change is significant. Update vertex_values
              vertex_values[iter][v] = new_value;
              // Set active for next iteration.
              frontier_curr[v] = 1;
            } else {
              // change is not significant. Copy vertex_values[iter-1]
              vertex_values[iter][v] = vertex_values[iter - 1][v];
            }
          }
          frontier_curr[v] =
              frontier_curr[v] ||
              forceActivateVertexForIteration(v, iter + 1, global_info);
          if (frontier_curr[v]) {
            if (use_source_contribution) {
              // update source_contrib for next iteration
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, source_change_in_contribution[v],
                  vertex_values[iter - 1][v], vertex_values[iter][v],
                  global_info);
            } else {
              source_change_in_contribution[v] =
                  aggregationValueIdentity<AggregationValueType>();
            }
          }
        }
        phase_time = phase_timer.stop();
        adaptive_executor.updateVertexMapTime(iter, phase_time);

        vertexSubset temp_vs(n, frontier_curr);
        frontier_curr_vs = temp_vs;
        misc_time += phase_timer.next();
        iteration_time = iteration_timer.stop();

        if (ae_enabled && iter == 1) {
          adaptive_executor.setApproximateTimeForCurrIter(iteration_time);
        }
        // Convergence check
        converged_iteration = iter;
        if (frontier_curr_vs.isEmpty()) {
          break;
        }
      }
    }
    if (ae_enabled) {
      adaptive_executor.updateEquation(converged_iteration);
    }
    return converged_iteration;
  }

  // ======================================================================
  // DELTACOMPUTE
  // ======================================================================
  void deltaCompute(edgeArray &edge_additions, edgeArray &edge_deletions) {
    timer iteration_timer, phase_timer, full_timer, pre_compute_timer;
    double misc_time, copy_time, phase_time, iteration_time, pre_compute_time;
    iteration_time = 0;
    // full_timer.start();

    // TODO : Realloc addition of new vertices
    n_old = n;
    if (edge_additions.maxVertex >= n) {
      processVertexAddition(edge_additions.maxVertex);
    }

    // Reset values before incremental computation
    parallel_for(uintV v = 0; v < n; v++) {
      frontier_curr[v] = 0;
      frontier_next[v] = 0;
      changed[v] = 0;

      vertex_value_old_prev[v] = vertexValueIdentity<VertexValueType>();
      vertex_value_old_curr[v] = vertexValueIdentity<VertexValueType>();
      // init vertex_value_old_next as iter=0 values
      initializeVertexValue<VertexValueType>(v, vertex_value_old_next[v],
                                             global_info);

      delta[v] = aggregationValueIdentity<AggregationValueType>();
      if (use_source_contribution) {
        source_change_in_contribution[v] =
            aggregationValueIdentity<AggregationValueType>();
      }
    }

    // ==================== UPDATE GLOBALINFO ===============================
    // deltaCompute/initCompute Save a copy of global_info before we lose any
    // relevant information of the old graph For example, In PageRank, we need
    // to save the outDegree for all vertices corresponding to the old graph
    global_info_old.copy(global_info);

    // Update global_info based on edge additions or deletions. This is
    // application specific. For example, for pagerank, the the outDegree of
    // vertices with edgeAddition will increase and those with edgeDeletions
    // will decrease
    global_info.processUpdates(edge_additions, edge_deletions);

    // ========== EDGE COMPUTATION - DIRECT CHANGES - for first iter ==========
    pre_compute_timer.start();

    // 从这里开始统计迭代时间，忽略图更新时间
    full_timer.start();

    // 如果是第一轮迭代，把新增边的source的上批次contribution算出来（相对于0），并传播到新增边的destination的delta，累加
    // Handle add edges changes, compute destination delta value
    parallel_for(long i = 0; i < edge_additions.size; i++) {
      uintV source = edge_additions.E[i].source;
      uintV destination = edge_additions.E[i].destination;

      // Update frontier and changed values
      hasSourceChangedByUpdate(source, edge_addition_enum,
                               frontier_curr[source], changed[source],
                               global_info, global_info_old);
      hasSourceChangedByUpdate(destination, edge_addition_enum,
                               frontier_curr[destination], changed[destination],
                               global_info, global_info_old);
      if (forceActivateVertexForIteration(source, 1, global_info_old)) {

        if (frontier_curr[source]) {
          // ensure to re-compute source and generate latest source_change_in_contribution[source]
          changed[source] = true;
        }
        if (frontier_curr[destination]) {
          // a bug, replace as following:
          changed[destination] = true;
          // changed[source] = true;
        }

        AggregationValueType contrib_change;
        if (use_source_contribution) {
          sourceChangeInContribution<AggregationValueType, VertexValueType,
                                     GlobalInfoType>(
              source, contrib_change, vertexValueIdentity<VertexValueType>(),
              vertex_values[0][source], global_info_old);
        }

// Do repropagate for edge source->destination.
#ifdef EDGEDATA
        EdgeData *edge_data = edge_additions.E[i].edgeData;
#else
        EdgeData *edge_data = &emptyEdgeData;
#endif
        bool ret =
            edgeFunction(source, destination, *edge_data,
                         vertex_values[0][source], contrib_change, global_info);
        if (ret) {
          if (use_lock) {
            vertex_locks[destination].writeLock();
            addToAggregation(contrib_change, delta[destination],
                             global_info_old);
            vertex_locks[destination].unlock();
          } else {
            addToAggregationAtomic(contrib_change, delta[destination],
                                   global_info_old);
          }
          if (!changed[destination])
            changed[destination] = true;
        }
      }
    }

    // 如果是第一轮迭代，把删除边的source的上批次contribution算出来（相对于0），并传播到删除边的destination的delta，累减
    // Handle delete edges changes, compute destination delta value
    parallel_for(long i = 0; i < edge_deletions.size; i++) {
      uintV source = edge_deletions.E[i].source;
      uintV destination = edge_deletions.E[i].destination;

      hasSourceChangedByUpdate(source, edge_deletion_enum,
                               frontier_curr[source], changed[source],
                               global_info, global_info_old);
      hasSourceChangedByUpdate(destination, edge_deletion_enum,
                               frontier_curr[destination], changed[destination],
                               global_info, global_info_old);
      if (forceActivateVertexForIteration(source, 1, global_info_old)) {
        // Update frontier and changed values
        if (frontier_curr[source]) {
          // ensure to re-compute source and generate latest source_change_in_contribution[source]
          changed[source] = true;
        }
        if (frontier_curr[destination]) {
          // a bug, replace as following:
          changed[destination] = true;
          // changed[source] = true;
        }

        AggregationValueType contrib_change;
        if (use_source_contribution) {
          sourceChangeInContribution<AggregationValueType, VertexValueType,
                                     GlobalInfoType>(
              source, contrib_change, vertexValueIdentity<VertexValueType>(),
              vertex_values[0][source], global_info_old);
        }

// Do retract for edge source->destination
#ifdef EDGEDATA
        EdgeData *edge_data = edge_deletions.E[i].edgeData;
#else
        EdgeData *edge_data = &emptyEdgeData;
#endif
        bool ret = edgeFunction(source, destination, *edge_data,
                                vertex_values[0][source], contrib_change,
                                global_info_old);
        if (ret) {
          if (use_lock) {
            vertex_locks[destination].writeLock();
            removeFromAggregation(contrib_change, delta[destination],
                                  global_info_old);
            vertex_locks[destination].unlock();
          } else {
            removeFromAggregationAtomic(contrib_change, delta[destination],
                                        global_info_old);
          }
          if (!changed[destination])
            changed[destination] = true;
        }
      }
    }
    pre_compute_time = pre_compute_timer.stop();

    // =============== INCREMENTAL COMPUTE - REFINEMENT START ================
    vertexSubset frontier_curr_vs(n, frontier_curr);
    bool should_switch_now = false;
    bool use_delta = true;

    if (ae_enabled && shouldSwitch(0, 0)) {
      should_switch_now = true;
    }

    // 开始迭代重放
    // The core of delta compute
    for (int iter = 1; iter < max_iterations; iter++) {
      // Perform switch if needed
      if (should_switch_now) {
        // 增量计算完毕，切换为全量计算
        converged_iteration = performSwitch(iter);
        break;
      }

      // initialize timers
      {
        iteration_timer.start();
        phase_timer.start();
        iteration_time = 0;
        misc_time = 0;
        copy_time = 0;
      }
      use_delta = shouldUseDelta(iter);

      // ================ COPY - PREPARE CURRENT ITERATION ================
      // Prepare old vertex values or switch to base graph compute
      {
        VertexValueType *temp1 = vertex_value_old_prev;
        vertex_value_old_prev = vertex_value_old_curr;
        vertex_value_old_curr = vertex_value_old_next;
        vertex_value_old_next = temp1;

        if (iter <= converged_iteration) {
          parallel_for(uintV v = 0; v < n; v++) {
            vertex_value_old_next[v] = vertex_values[iter][v];
          }
        } else {
          converged_iteration = performSwitch(iter);
          break;
        }
      }
      copy_time += phase_timer.next();
      // ========== EDGE COMPUTATION - TRANSITIVE CHANGES ==========
      // 第一轮迭代时，计算source的contribution变化（当前批的contribution - 上一批的contribution）
      // Reset source contribution for first iteration
      if ((use_source_contribution) && (iter == 1)) {
        parallel_for(uintV u = 0; u < n; u++) {
          if (frontier_curr[u]) {
            // compute source change in contribution
            // compute current contribution
            // vertex_values[iter - 1] point to vertex values in last iter on current graph
            AggregationValueType contrib_change =
                aggregationValueIdentity<AggregationValueType>();
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                u, contrib_change, vertexValueIdentity<VertexValueType>(),
                vertex_values[iter - 1][u], global_info);
            addToAggregation(contrib_change, source_change_in_contribution[u],
                             global_info);

            // retract previous contribution
            // vertex_value_old_curr point to vertex values in last iter on previous graph
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                u, contrib_change, vertexValueIdentity<VertexValueType>(),
                vertex_value_old_curr[u], global_info_old);
            removeFromAggregation(
                contrib_change, source_change_in_contribution[u], global_info);
          }
        }
      }

      // 把source的contribution变化累加传播到所有destination的delta（不用考虑destination是新的点还是已有的点，前边以后初始化好destination的delta了）
      // Compute destination delta value
      parallel_for(uintV u = 0; u < n; u++) {
        // frontier_curr[u] represents the add/delete edge's source
        if (frontier_curr[u]) {
          // check for propagate and retract for the vertices.
          intE outDegree = my_graph.V[u].getOutDegree();

          parallel_for(int i = 0; i < outDegree; i++) {
          // replace for debug: granular_for(i, 0, outDegree, (outDegree > 1024), {
            uintV v = my_graph.V[u].getOutNeighbor(i);
            bool ret = false;
            AggregationValueType contrib_change =
                use_source_contribution
                    ? source_change_in_contribution[u]
                    : aggregationValueIdentity<AggregationValueType>();

#ifdef EDGEDATA
            EdgeData *edge_data = my_graph.V[u].getOutEdgeData(i);
#else
            EdgeData *edge_data = &emptyEdgeData;
#endif
            ret = edgeFunction(u, v, *edge_data, vertex_values[iter - 1][u],
                               contrib_change, global_info);

            if (ret) {
              if (use_lock) {
                vertex_locks[v].writeLock();
                if (ret) {
                  addToAggregation(contrib_change, delta[v], global_info);
                }
                vertex_locks[v].unlock();

              } else {
                if (ret) {
                  addToAggregationAtomic(contrib_change, delta[v], global_info);
                }
              }
              if (!changed[v])
                changed[v] = 1;
            }
          //replace for debug: });
          }
        }
      }
      phase_time = phase_timer.next();

      // ========== VERTEX COMPUTATION  ==========
      bool use_delta_next_iteration = shouldUseDelta(iter + 1);
      parallel_for(uintV v = 0; v < n; v++) {
        // changed vertices need to be processed
        frontier_curr[v] = 0;
        if ((v >= n_old) && (changed[v] == false)) {
          changed[v] = forceComputeVertexForIteration(v, iter, global_info);
        }

        if (changed[v]) {
          frontier_curr[v] = 0;

          // 把delta值更新到aggregation_values
          // delta has the current cumulative change for the vertex.
          // Update the aggregation value in history
          addToAggregation(delta[v], aggregation_values[iter][v], global_info);

          // 根据aggregation_values执行compute计算
          VertexValueType new_value;
          computeFunction(v, aggregation_values[iter][v],
                          vertex_values[iter - 1][v], new_value, global_info);

          if (forceActivateVertexForIteration(v, iter + 1, global_info)) {
            frontier_curr[v] = 1;
          }

          // 重新计算目标点的contribution变化
          AggregationValueType contrib_change =
              aggregationValueIdentity<AggregationValueType>();
          source_change_in_contribution[v] =
              aggregationValueIdentity<AggregationValueType>();

          // Merge curr iter agg value
          if (notDelZero(new_value, vertex_values[iter - 1][v], global_info)) {
            // change is significant. Update vertex_values
            vertex_values[iter][v] = new_value;
            frontier_curr[v] = 1;
            if (use_delta_next_iteration) {
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, contrib_change, vertex_values[iter - 1][v],
                  vertex_values[iter][v], global_info);
            } else {
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, contrib_change, vertexValueIdentity<VertexValueType>(),
                  vertex_values[iter][v], global_info);
            }
            addToAggregation(contrib_change, source_change_in_contribution[v],
                             global_info);


          } else {
            // change is not significant. Copy vertex_values[iter-1]
            vertex_values[iter][v] = vertex_values[iter - 1][v];
          }

          // Retract next iter source_change_in_contribution to agg on previous graph
          if (notDelZero(vertex_value_old_next[v], vertex_value_old_curr[v],
                        global_info_old)) {
            // change is significant. Update v_change
            frontier_curr[v] = 1;
            if (use_delta_next_iteration) {
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, contrib_change, vertex_value_old_curr[v],
                  vertex_value_old_next[v], global_info_old);
            } else {

              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, contrib_change, vertexValueIdentity<VertexValueType>(),
                  vertex_value_old_next[v], global_info_old);
            }
            removeFromAggregation(contrib_change,
                                  source_change_in_contribution[v],
                                  global_info_old);
          }
        }
      }
      phase_time = phase_timer.next();

      // ========== EDGE COMPUTATION - DIRECT CHANGES - for next iter ==========
      bool has_direct_changes = false;

      // 准备下次迭代的destination的delta
      // Compute add edges' next iter destination delta
      parallel_for(long i = 0; i < edge_additions.size; i++) {
        uintV source = edge_additions.E[i].source;
        uintV destination = edge_additions.E[i].destination;
        AggregationValueType contrib_change;

        // The edge source is active(not converged last iter or force-activated next iter), compute destination delta
        if (notDelZero(vertex_value_old_curr[source],
                      vertex_value_old_next[source], global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old))) {
          if (use_delta_next_iteration) {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, vertex_value_old_curr[source],
                vertex_value_old_next[source], global_info_old);
          } else {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, vertexValueIdentity<VertexValueType>(),
                vertex_value_old_next[source], global_info_old);
          }
// Do repropagate for edge source->destination.
#ifdef EDGEDATA
          EdgeData *edge_data = edge_additions.E[i].edgeData;
#else
          EdgeData *edge_data = &emptyEdgeData;
#endif
          bool ret = edgeFunction(source, destination, *edge_data,
                                  vertex_values[0][source], contrib_change,
                                  global_info);

          if (ret) {
            if (use_lock) {
              vertex_locks[destination].writeLock();
              addToAggregation(contrib_change, delta[destination],
                               global_info_old);
              vertex_locks[destination].unlock();
            } else {
              addToAggregationAtomic(contrib_change, delta[destination],
                                     global_info_old);
            }
            if (!changed[destination])
              changed[destination] = 1;
            if (!has_direct_changes)
              has_direct_changes = true;
          }
        }
      }

      // Compute delete edges' next iter destination delta
      parallel_for(long i = 0; i < edge_deletions.size; i++) {
        uintV source = edge_deletions.E[i].source;
        uintV destination = edge_deletions.E[i].destination;
        AggregationValueType contrib_change;

        if (notDelZero(vertex_value_old_curr[source],
                      vertex_value_old_next[source], global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old))) {
          // Do repropagate for edge source->destination.
          if (use_delta_next_iteration) {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, vertex_value_old_curr[source],
                vertex_value_old_next[source], global_info_old);
          } else {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, vertexValueIdentity<VertexValueType>(),
                vertex_value_old_next[source], global_info_old);
          }
#ifdef EDGEDATA
          EdgeData *edge_data = edge_deletions.E[i].edgeData;
#else
          EdgeData *edge_data = &emptyEdgeData;
#endif
          bool ret = edgeFunction(source, destination, *edge_data,
                                  vertex_values[0][source], contrib_change,
                                  global_info);

          if (ret) {
            if (use_lock) {
              vertex_locks[destination].writeLock();
              removeFromAggregation(contrib_change, delta[destination],
                                    global_info_old);
              vertex_locks[destination].unlock();

            } else {
              removeFromAggregationAtomic(contrib_change, delta[destination],
                                          global_info_old);
            }
            if (!changed[destination])
              changed[destination] = 1;
            if (!has_direct_changes)
              has_direct_changes = true;
          }
        }
      }
      phase_time = phase_timer.next();

      vertexSubset temp_vs(n, frontier_curr);
      frontier_curr_vs = temp_vs;

      misc_time += phase_timer.next();
      iteration_time = iteration_timer.next();

      // Convergence check
      if (!has_direct_changes && frontier_curr_vs.isEmpty()) {
        // There are no more active vertices
        if (iter == converged_iteration) {
          break;
        } else if (iter > converged_iteration) {
          assert(("Missed switching to Traditional incremental computing when "
                  "iter == converged_iter",
                  false));
        } else {
          // Values stable for the changed vertices at this iteration.
          // But, the changed vertices might receive new changes. So,
          // continue loop until iter == converged_iteration vertices may
          // still not have converged. So, keep continuing until
          // converged_iteration is reached.
        }
      }
      if (iter == 1) {
        iteration_time += pre_compute_time;
      }

      if (ae_enabled && shouldSwitch(iter, iteration_time)) {
        should_switch_now = true;
      }
      misc_time += phase_timer.stop();
      iteration_time += iteration_timer.stop();
    }

    cout << "Finished batch : " << full_timer.stop() << "\n";
    cout << "Number of iterations : " << converged_iteration << "\n";
    // testPrint();
    printOutput();
  }

  // Refactor this in a better way
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::my_graph;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::config;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::max_iterations;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::history_iterations;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::converged_iteration;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::use_lock;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_locks;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::aggregation_values;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_values;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::n;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::global_info;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::delta;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::use_source_contribution;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::source_change_in_contribution;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::n_old;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::global_info_old;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_value_old_next;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_value_old_curr;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_value_old_prev;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::all;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::frontier_curr;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::frontier_next;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::changed;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::ingestor;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::current_batch;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::adaptive_executor;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::ae_enabled;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::testPrint;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::printOutput;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::shouldSwitch;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::performSwitch;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::processVertexAddition;
};
#endif

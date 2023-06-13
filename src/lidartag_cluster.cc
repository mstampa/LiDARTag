/* Copyright (C) 2013-2020, The Regents of The University of Michigan.
 * All rights reserved.
 * This software was developed in the Biped Lab (https://www.biped.solutions/)
 * under the direction of Jessy Grizzle, grizzle@umich.edu. This software may
 * be available under alternative licensing terms; contact the address above.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the Regents of The University of Michigan.
 *
 * AUTHOR: Bruce JK Huang (bjhuang@umich.edu)
 * WEBSITE: https://www.brucerobot.com/
 */

#include "lidartag.h"
#include "ultra_puck.h"

#include <ros/package.h>

#include <iostream>

using namespace std;

namespace BipedLab
{
/*
 * A function to cluster a single point into a new cluster or an existing cluster
 */
void LiDARTag::_clusterClassifier(const LiDARPoints_t& point, vector<ClusterFamily_t>& cluster_buff)
{
  // The first time to cluster the point cloud
  int ValidCluster = 1;  // Marker every cluster is valid and will be checked again later
  if (cluster_buff.size() == 0)
  {
    int bottom_ring = point.point.ring;
    int top_ring = point.point.ring;
    PointXYZIRT top_most_point = point.point;
    top_most_point.z = top_most_point.z + _linkage_threshold;
    PointXYZIRT bottom_most_point = point.point;
    bottom_most_point.z -= _linkage_threshold;

    PointXYZIRT front_most_point = point.point;
    front_most_point.x += _linkage_threshold;

    PointXYZIRT back_most_point = point.point;
    back_most_point.x -= _linkage_threshold;

    PointXYZIRT right_most_point = point.point;
    right_most_point.y -= _linkage_threshold;
    PointXYZIRT left_most_point = point.point;
    left_most_point.y += _linkage_threshold;

    ClusterFamily_t current_cluster = { 0,
                                        ValidCluster,
                                        top_ring,
                                        bottom_ring,
                                        top_most_point,
                                        bottom_most_point,
                                        front_most_point,
                                        back_most_point,
                                        right_most_point,
                                        left_most_point,
                                        point.point };

    MaxMin_t initial_value;  // = {(int)1e8, 0, -1};
    initial_value.min = (int)1e8;
    initial_value.average = (int)-1;
    initial_value.max = (int)-1;

    current_cluster.max_min_index_of_each_ring.resize(_beam_num, initial_value);
    current_cluster.max_min_index_of_each_ring[point.point.ring].max = point.index;
    current_cluster.max_min_index_of_each_ring[point.point.ring].min = point.index;

    current_cluster.max_intensity = point.point;
    current_cluster.min_intensity = point.point;

    current_cluster.edge_points.push_back(point);
    cluster_buff.push_back(current_cluster);
    return;
  }
  else
  {
    // Take a point to go through all the existing cluster to see if this
    // point belongs to any of them
    // Once it is found belonging to one of the clusters then return.
    // After looping though and couldn't find a belonging group then add it
    // to a new cluster
    TestCluster_t* new_cluster = new TestCluster_t{ 0 };

    for (int i = 0; i < cluster_buff.size(); ++i)
    {
      _updateCluster(point, cluster_buff[i], new_cluster);
      if (!(new_cluster->flag))
      {
        delete new_cluster;

        return;
      }
    }
    // Add a new cluster
    if (new_cluster->flag)
    {
      int Newcluster_id = cluster_buff.size();
      int top_ring = point.point.ring;
      int bottom_ring = point.point.ring;

      PointXYZIRT top_most_point = point.point;
      top_most_point.z += _linkage_threshold;
      PointXYZIRT bottom_most_point = point.point;
      bottom_most_point.z -= _linkage_threshold;

      PointXYZIRT front_most_point = point.point;
      front_most_point.x += _linkage_threshold;
      PointXYZIRT back_most_point = point.point;
      back_most_point.x -= _linkage_threshold;

      PointXYZIRT right_most_point = point.point;
      right_most_point.y -= _linkage_threshold;
      PointXYZIRT left_most_point = point.point;
      left_most_point.y += _linkage_threshold;

      new_cluster->new_cluster = { Newcluster_id,    ValidCluster,      top_ring,         bottom_ring,
                                   top_most_point,   bottom_most_point, front_most_point, back_most_point,
                                   right_most_point, left_most_point,   point.point };

      // To fill in points between initial points and end points later
      MaxMin_t initial_value;
      initial_value.min = (int)1e8;
      initial_value.average = (int)0;
      initial_value.max = (int)-1;

      new_cluster->new_cluster.max_min_index_of_each_ring.resize(_beam_num, initial_value);
      new_cluster->new_cluster.max_min_index_of_each_ring[point.point.ring].max = point.index;
      new_cluster->new_cluster.max_min_index_of_each_ring[point.point.ring].min = point.index;

      new_cluster->new_cluster.max_intensity = point.point;
      new_cluster->new_cluster.min_intensity = point.point;

      new_cluster->new_cluster.edge_points.push_back(point);
      cluster_buff.push_back(new_cluster->new_cluster);
    }
    delete new_cluster;
  }
}

/*
 * A function update some information about a cluster if this point belongs to
 * this cluster; if not belonging to this cluster then return and create a new
 * one
 */
void LiDARTag::_updateCluster(const LiDARPoints_t& point, ClusterFamily_t& old_cluster, TestCluster_t* new_cluster)
{
  // This point is outside of the current cluster
  if (!_isWithinCluster(point, old_cluster))
  {
    new_cluster->flag = 1;
    return;
  }
  else
  {
    new_cluster->flag = 0;
  }

  // This point is inside this cluster
  if (!new_cluster->flag)
  {
    // update the boundary of the old cluster
    if (point.point.ring < old_cluster.bottom_ring)
    {
      old_cluster.bottom_ring = point.point.ring;
    }
    if (point.point.ring > old_cluster.top_ring)
    {
      old_cluster.top_ring = point.point.ring;
    }
    if (point.point.x + _linkage_threshold > old_cluster.front_most_point.x)
    {
      old_cluster.front_most_point = point.point;
      old_cluster.front_most_point.x += _linkage_threshold;
    }
    if (point.point.x - _linkage_threshold < old_cluster.back_most_point.x)
    {
      old_cluster.back_most_point = point.point;
      old_cluster.back_most_point.x -= _linkage_threshold;
    }

    if (point.point.y + _linkage_threshold > old_cluster.left_most_point.y)
    {
      old_cluster.left_most_point = point.point;
      old_cluster.left_most_point.y += _linkage_threshold;
    }
    if (point.point.y - _linkage_threshold < old_cluster.right_most_point.y)
    {
      old_cluster.right_most_point = point.point;
      old_cluster.right_most_point.y -= _linkage_threshold;
    }

    if (point.point.z + _linkage_threshold > old_cluster.top_most_point.z)
    {
      old_cluster.top_most_point = point.point;
      old_cluster.top_most_point.z += _linkage_threshold;
    }
    if (point.point.z - _linkage_threshold < old_cluster.bottom_most_point.z)
    {
      old_cluster.bottom_most_point = point.point;
      old_cluster.bottom_most_point.z -= _linkage_threshold;
    }

    // update the max/min index of each ring in this cluster
    if (old_cluster.max_min_index_of_each_ring[point.point.ring].max < point.index)
      old_cluster.max_min_index_of_each_ring[point.point.ring].max = point.index;

    if (old_cluster.max_min_index_of_each_ring[point.point.ring].min > point.index)
      old_cluster.max_min_index_of_each_ring[point.point.ring].min = point.index;

    // update the max/min intensity of this cluster
    if (old_cluster.max_intensity.intensity < point.point.intensity)
      old_cluster.max_intensity = point.point;

    if (old_cluster.min_intensity.intensity > point.point.intensity)
      old_cluster.min_intensity = point.point;

    old_cluster.edge_points.push_back(point);
  }
}

bool LiDARTag::_isWithinCluster(const LiDARPoints_t& point, ClusterFamily_t& cluster)
{
  return (point.point.ring == cluster.bottom_ring || point.point.ring == (cluster.bottom_ring - 1)) &&
         _isWithinClusterHorizon(point, cluster, _linkage_threshold);
}

bool LiDARTag::_isWithinClusterHorizon(const LiDARPoints_t& point, ClusterFamily_t& cluster, double threshold)
{
  return (point.point.x < cluster.front_most_point.x + threshold) &&
         (cluster.back_most_point.x - threshold < point.point.x) &&
         (point.point.y < cluster.left_most_point.y + threshold) &&
         (cluster.right_most_point.y - _linkage_threshold < point.point.y);
}

}  // namespace BipedLab

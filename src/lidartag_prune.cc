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

#include "ultra_puck.h"
#include "lidartag.h"

#include <pcl/point_types.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <ros/package.h>
#include <velodyne_pcl/point_types.h>
#include <velodyne_pointcloud/pointcloudXYZIRT.h>

#include <cmath>
#include <fstream>

using namespace std;

namespace BipedLab
{
/*
 * A valid cluster, valid tag, the points from the original point cloud that belong to the cluster
 * could be estimated from the LiDAR system
 * Therefore, if the points on the tag is too less, which means it is not a valid
 * tag where it might just a shadow of a valid tag
 */
bool LiDARTag::_clusterPointsCheck(ClusterFamily_t& Cluster)
{
  auto distance = sqrt(pow(Cluster.average.x, 2) + pow(Cluster.average.y, 2));
  int maxPoints = LiDARTag::_areaPoints(distance, _payload_size * M_SQRT2, _payload_size * M_SQRT2);
  int minPoints = LiDARTag::_areaPoints(distance, _payload_size / M_SQRT2, _payload_size / M_SQRT2);

  return Cluster.data.size() >= minPoints;
}

/*
 * A function to get a number of points on a given-distance tag or object
 */
int LiDARTag::_areaPoints(const double& Distance, const double& ObjWidth, const double& ObjHeight)
{
  int NumOfHorizontalPoints = ceil(ObjWidth / (Distance * tan(0.1 * M_PI / 180)));
  double HalfVerticalAngle = atan((ObjHeight / 2) / abs(Distance)) * 180 / M_PI;

  int NumOfVerticalRing = 0;
  for (int i = 0; i < UltraPuckV2::beams; ++i)
  {
    if (HalfVerticalAngle > abs(UltraPuckV2::el[i]))
    {
      NumOfVerticalRing++;
    }
  }

  return NumOfVerticalRing * NumOfHorizontalPoints;
}

/*
 * A function to calculate the upper bound of points that can exist in a cluster
 * based on the payload size
 */
void LiDARTag::_maxPointsCheck(ClusterFamily_t& cluster)
{
  int ring = round(_beam_num / 2);
  double point_resolution = 2 * M_PI / _LiDAR_system.ring_average_table[ring].average;
  auto distance = sqrt(pow(cluster.average.x, 2) + pow(cluster.average.y, 2) + pow(cluster.average.z, 2));
  float payload_w = M_SQRT2 * _payload_size;
  int num_horizontal_points = ceil(payload_w / (distance * sin(point_resolution)));
  int num_vertical_ring = abs(cluster.top_ring - cluster.bottom_ring) + 1;
  int expected_points = num_vertical_ring * num_horizontal_points;

  if ((cluster.data.size() + cluster.edge_points.size()) > expected_points)
  {
    _result_statistics.cluster_removal.maximum_return++;
    _result_statistics.remaining_cluster_size--;
    if (_mark_cluster_validity)
      cluster.valid = false;
  }

  if (_debug_info)
  {
    ROS_DEBUG_STREAM("==== _maxPointsCheck ====");
    ROS_DEBUG_STREAM("Distance : " << distance << ", num_horizontal_points: " << num_horizontal_points);
    ROS_DEBUG_STREAM("Expected Points: " << expected_points);
    ROS_DEBUG_STREAM("Actual Points: " << cluster.data.size() + cluster.edge_points.size());
    ROS_DEBUG_STREAM("Status: " << (cluster.data.size() + cluster.edge_points.size() <= expected_points));
  }
}

/*
 * Fit a plane to a cluster. Returns false if unable to estimate a plane.
 * Otherwise, returns the number of inliers and the coefficients of the plane.
 */
bool LiDARTag::_rejectWithPlanarCheck(ClusterFamily_t& cluster, pcl::PointIndices::Ptr inliers,
                                      pcl::ModelCoefficients::Ptr coefficients, std::ostream& fplanefit)
{
  // Convert cluster family into pcl point cloud
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
  cloud->points.resize(cluster.data.size() + cluster.edge_points.size());
  for (std::size_t i = 0; i < cluster.edge_points.size(); ++i)
  {
    cloud->points[i].x = cluster.edge_points[i].point.x;
    cloud->points[i].y = cluster.edge_points[i].point.y;
    cloud->points[i].z = cluster.edge_points[i].point.z;
  }
  for (std::size_t i = cluster.edge_points.size(); i < cloud->points.size(); ++i)
  {
    cloud->points[i].x = cluster.data[i - cluster.edge_points.size()].point.x;
    cloud->points[i].y = cluster.data[i - cluster.edge_points.size()].point.y;
    cloud->points[i].z = cluster.data[i - cluster.edge_points.size()].point.z;
  }

  // Create segmentation object
  pcl::SACSegmentation<pcl::PointXYZ> seg;
  // Optional
  seg.setOptimizeCoefficients(true);
  // Mandatory
  seg.setModelType(pcl::SACMODEL_PLANE);
  seg.setMethodType(pcl::SAC_RANSAC);
  seg.setDistanceThreshold(_distance_to_plane_threshold);

  seg.setInputCloud(cloud);
  seg.segment(*inliers, *coefficients);
  cluster.inliers = inliers->indices.size();
  if (_debug_info)
  {
    ROS_DEBUG_STREAM("==== _rejectWithPlanarCheck ====");
    float distance = std::sqrt(pow(cluster.average.x, 2) + pow(cluster.average.y, 2) + pow(cluster.average.z, 2));
    ROS_DEBUG_STREAM("Distance: " << distance);
    ROS_DEBUG_STREAM("Actual Points: " << cluster.data.size() + cluster.edge_points.size());
    ROS_DEBUG_STREAM("Inliers: " << inliers->indices.size());
    ROS_DEBUG_STREAM("Outliers: " << cluster.data.size() - inliers->indices.size());
  }

  if (inliers->indices.size() == 0)
  {
    if (_debug_info)
    {
      ROS_DEBUG_STREAM("Status: " << false);
      ROS_WARN_STREAM("Failed to fit a plane model to the cluster.");
    }
    _result_statistics.cluster_removal.plane_fitting++;
    _result_statistics.remaining_cluster_size--;
    return false;
  }

  if (_debug_info)
    ROS_DEBUG_STREAM("Status: " << true);

  if (_log_data)
  {
    fplanefit << "Successfully fit plane!" << endl;
    fplanefit << "Cluster Size: " << cluster.data.size() << endl;
    fplanefit << "Inliers     : " << inliers->indices.size() << endl;
    fplanefit << "Outliers    : " << cluster.data.size() - inliers->indices.size() << endl;
  }

  return true;
}

}  // namespace BipedLab

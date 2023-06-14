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

#include <pcl/common/intersections.h>
#include <ros/package.h>

#include <iostream>
#include <string>

using namespace std;

namespace BipedLab
{
/*
 * A function to draw a point in rviz
 */
void LiDARTag::_assignMarker(visualization_msgs::Marker& Marker, const uint32_t Shape, const string NameSpace,
                             const double r, const double g, const double b, const PointXYZIRT& point, const int Count,
                             const double Size, const string Text)
{
  Marker.header.frame_id = _pub_frame;
  Marker.header.stamp = _current_scan_time;
  Marker.ns = NameSpace;
  Marker.id = Count;
  Marker.type = Shape;
  Marker.action = visualization_msgs::Marker::ADD;
  Marker.pose.position.x = point.x;
  Marker.pose.position.y = point.y;
  Marker.pose.position.z = point.z;
  Marker.pose.orientation.x = 0.0;
  Marker.pose.orientation.y = 0.0;
  Marker.pose.orientation.z = 0.0;
  Marker.pose.orientation.w = 1.0;
  Marker.text = Text;
  Marker.lifetime = ros::Duration(_sleep_time_for_vis);  // should disappear along with updateing rate
  Marker.scale.x = Size;
  Marker.scale.y = Size;
  Marker.scale.z = Size;
  Marker.color.r = r;
  Marker.color.g = g;
  Marker.color.b = b;
  Marker.color.a = 1.0;
}

void LiDARTag::_assignVectorMarker(visualization_msgs::Marker& Marker, const uint32_t Shape, const string NameSpace,
                                   const double r, const double g, const double b, const int Count, const double Size,
                                   Eigen::Vector3f edge_vector, const PointXYZIRT& centriod, const string Text)
{
  Marker.header.frame_id = _pub_frame;
  Marker.header.stamp = _current_scan_time;
  Marker.ns = NameSpace;
  Marker.id = Count;
  Marker.type = Shape;
  Marker.action = visualization_msgs::Marker::ADD;

  Marker.text = Text;
  Marker.lifetime = ros::Duration(_sleep_time_for_vis);  // should disappear along with updateing rate
  geometry_msgs::Point p1;
  p1.x = 0;  // centroid.x;
  p1.y = 0;  // centroid.y;
  p1.z = 0;  // centroid.z;
  Marker.points.push_back(p1);

  geometry_msgs::Point p2;
  p2.x = edge_vector[0];
  p2.y = edge_vector[1];
  p2.z = edge_vector[2];
  Marker.points.push_back(p2);
  Marker.color.r = r;
  Marker.color.g = g;
  Marker.color.b = b;
  Marker.scale.x = 0.02;
  Marker.scale.y = Size;
  Marker.scale.z = Size;
}
void LiDARTag::_plotIdealFrame()
{
  visualization_msgs::MarkerArray FrameMarkArray;
  for (int k = 0; k < _tag_size_list.size(); ++k)
  {
    visualization_msgs::Marker line_list;
    line_list.id = 0;
    line_list.header = _point_cloud_header;
    line_list.type = visualization_msgs::Marker::LINE_LIST;
    line_list.ns = "tag_size_" + to_string(_tag_size_list[k]);
    line_list.color.g = 1.0;
    line_list.color.a = 1.0;
    line_list.scale.x = 0.01;

    vector<vector<int>> vertex = { { 1, 1 }, { 1, -1 }, { -1, -1 }, { -1, 1 } };
    for (int i = 0; i < 4; ++i)
    {
      vector<int> v = vertex[i];
      geometry_msgs::Point p;
      p.x = -0.02;
      p.y = v[0] * _tag_size_list[k] / 2;  //_payload_size
      p.z = v[1] * _tag_size_list[k] / 2;
      line_list.points.push_back(p);
      p.x = 0.02;
      line_list.points.push_back(p);
    }

    for (int j = -1; j <= 1; j += 2)
    {
      for (int i = 0; i < 4; ++i)
      {
        vector<int> v = vertex[i];
        geometry_msgs::Point p;
        p.x = j * 0.02;
        p.y = v[0] * _tag_size_list[k] / 2;
        p.z = v[1] * _tag_size_list[k] / 2;
        line_list.points.push_back(p);
        v = vertex[(i + 1) % 4];
        p.y = v[0] * _tag_size_list[k] / 2;
        p.z = v[1] * _tag_size_list[k] / 2;
        line_list.points.push_back(p);
      }
    }
    FrameMarkArray.markers.push_back(line_list);
  }

  _ideal_frame_pub.publish(FrameMarkArray);
}

void LiDARTag::_plotTagFrame(const ClusterFamily_t& cluster)
{
  visualization_msgs::Marker line_list;
  line_list.id = cluster.cluster_id;
  line_list.header = _point_cloud_header;
  line_list.type = visualization_msgs::Marker::LINE_LIST;
  line_list.color.r = 0.0;
  line_list.color.g = 1.0;
  line_list.color.b = 1.0;
  line_list.color.a = 1.0;
  line_list.scale.x = 0.01;

  vector<vector<int>> vertex = { { 1, 1 }, { 1, -1 }, { -1, -1 }, { -1, 1 } };
  for (int i = 0; i < 4; ++i)
  {
    vector<int> v = vertex[i];
    Eigen::Vector4f corner_lidar1(-0.02, v[0] * cluster.tag_size / 2, v[1] * cluster.tag_size / 2, 1);
    Eigen::Vector4f corner_lidar2(0.02, v[0] * cluster.tag_size / 2, v[1] * cluster.tag_size / 2, 1);
    Eigen::Vector4f corner_tag1 = cluster.pose.homogeneous * corner_lidar1;
    Eigen::Vector4f corner_tag2 = cluster.pose.homogeneous * corner_lidar2;
    geometry_msgs::Point p;
    p.x = corner_tag1(0);
    p.y = corner_tag1(1);  //_payload_size
    p.z = corner_tag1(2);
    line_list.points.push_back(p);
    p.x = corner_tag2(0);
    p.y = corner_tag2(1);  //_payload_size
    p.z = corner_tag2(2);
    line_list.points.push_back(p);
  }

  for (int j = -1; j <= 1; j += 2)
  {
    for (int i = 0; i < 4; ++i)
    {
      vector<int> v = vertex[i];
      Eigen::Vector4f corner_lidar1(j * 0.02, v[0] * cluster.tag_size / 2, v[1] * cluster.tag_size / 2, 1);
      Eigen::Vector4f corner_tag1 = cluster.pose.homogeneous * corner_lidar1;
      geometry_msgs::Point p;
      p.x = corner_tag1(0);
      p.y = corner_tag1(1);
      p.z = corner_tag1(2);
      line_list.points.push_back(p);
      v = vertex[(i + 1) % 4];
      Eigen::Vector4f corner_lidar2(j * 0.02, v[0] * cluster.tag_size / 2, v[1] * cluster.tag_size / 2, 1);
      Eigen::Vector4f corner_tag2 = cluster.pose.homogeneous * corner_lidar2;
      p.x = corner_tag2(0);
      p.y = corner_tag2(1);
      p.z = corner_tag2(2);
      line_list.points.push_back(p);
    }
  }

  _tag_frame_pub.publish(line_list);
}

visualization_msgs::Marker LiDARTag::_visualizeVector(Eigen::Vector3f edge_vector, PointXYZIRT centriod, int ID)
{
  visualization_msgs::Marker edge;
  edge.id = ID;
  edge.header = _point_cloud_header;
  edge.type = visualization_msgs::Marker::LINE_STRIP;
  edge.color.g = 1.0;
  edge.color.a = 1.0;
  edge.scale.x = 0.02;

  geometry_msgs::Point p;
  p.x = centriod.x;
  p.y = centriod.y;
  p.z = centriod.z;
  edge.points.push_back(p);

  p.x += edge_vector[0];
  p.y += edge_vector[1];
  p.z += edge_vector[2];
  edge.points.push_back(p);

  return edge;
}

/*
 * A function to prepare for displaying results in rviz
 */
void LiDARTag::_clusterToPclVectorAndMarkerPublisher(
    const std::vector<ClusterFamily_t>& Cluster, pcl::PointCloud<PointXYZIRT>::Ptr OutCluster,
    pcl::PointCloud<PointXYZIRT>::Ptr OutEdgeCluster, pcl::PointCloud<PointXYZIRT>::Ptr OutPayload,
    pcl::PointCloud<PointXYZIRT>::Ptr OutPayload3D, pcl::PointCloud<PointXYZIRT>::Ptr OutTarget,
    pcl::PointCloud<PointXYZIRT>::Ptr OutInitialTarget, pcl::PointCloud<PointXYZIRT>::Ptr EdgeGroup1,
    pcl::PointCloud<PointXYZIRT>::Ptr EdgeGroup2, pcl::PointCloud<PointXYZIRT>::Ptr EdgeGroup3,
    pcl::PointCloud<PointXYZIRT>::Ptr EdgeGroup4, pcl::PointCloud<PointXYZIRT>::Ptr BoundaryPts,
    visualization_msgs::MarkerArray& ClusterArray)
{
  /* initialize random seed for coloring the marker*/
  srand(time(NULL));
  visualization_msgs::MarkerArray BoundMarkArray;
  visualization_msgs::MarkerArray PayloadMarkArray;
  visualization_msgs::MarkerArray IDMarkArray;

  // Used to identify between multiple clusters in a single point
  // cloud in the analysis file. The id being reset to 1 each time
  // the function is called is supposed to indicate in the output
  // file that the proceeding clusters belong to a new payload
  int cluster_pc_id = 1;

  int PointsInClusters = 0;
  int Clustercount = 0;
  for (int Key = 0; Key < Cluster.size(); ++Key)
  {
    if (Cluster[Key].valid != 1)
      continue;

    LiDARTag::_plotTagFrame(Cluster[Key]);

    // pick a random color for each cluster
    double r = (double)rand() / RAND_MAX;
    double g = (double)rand() / RAND_MAX;
    double b = (double)rand() / RAND_MAX;

    // Draw boundary marker of each cluster
    visualization_msgs::Marker BoundaryMarker;
    LiDARTag::_assignMarker(BoundaryMarker, visualization_msgs::Marker::CUBE,
                            "Boundary" + to_string(Cluster[Key].cluster_id), r, g, b, Cluster[Key].top_most_point, 0,
                            0.02);
    BoundMarkArray.markers.push_back(BoundaryMarker);
    LiDARTag::_assignMarker(BoundaryMarker, visualization_msgs::Marker::CUBE,
                            "Boundary" + to_string(Cluster[Key].cluster_id), r, g, b, Cluster[Key].bottom_most_point, 1,
                            0.02);
    BoundMarkArray.markers.push_back(BoundaryMarker);
    LiDARTag::_assignMarker(BoundaryMarker, visualization_msgs::Marker::CUBE,
                            "Boundary" + to_string(Cluster[Key].cluster_id), r, g, b, Cluster[Key].front_most_point, 2,
                            0.02);
    BoundMarkArray.markers.push_back(BoundaryMarker);
    LiDARTag::_assignMarker(BoundaryMarker, visualization_msgs::Marker::CUBE,
                            "Boundary" + to_string(Cluster[Key].cluster_id), r, g, b, Cluster[Key].back_most_point, 3,
                            0.02);
    BoundMarkArray.markers.push_back(BoundaryMarker);
    LiDARTag::_assignMarker(BoundaryMarker, visualization_msgs::Marker::CUBE,
                            "Boundary" + to_string(Cluster[Key].cluster_id), r, g, b, Cluster[Key].right_most_point, 4,
                            0.02);
    BoundMarkArray.markers.push_back(BoundaryMarker);
    LiDARTag::_assignMarker(BoundaryMarker, visualization_msgs::Marker::CUBE,
                            "Boundary" + to_string(Cluster[Key].cluster_id), r, g, b, Cluster[Key].left_most_point, 5,
                            0.02);
    BoundMarkArray.markers.push_back(BoundaryMarker);

    LiDARTag::_assignMarker(BoundaryMarker, visualization_msgs::Marker::SPHERE,
                            "AveragePoint" + to_string(Cluster[Key].cluster_id), 1, 0, 0, Cluster[Key].average, 1,
                            0.05);
    BoundMarkArray.markers.push_back(BoundaryMarker);

    LiDARTag::_assignMarker(
        BoundaryMarker, visualization_msgs::Marker::TEXT_VIEW_FACING, "Text" + to_string(Cluster[Key].cluster_id), 1, 1,
        1, Cluster[Key].average, 1, 0.05,
        string(to_string(Cluster[Key].cluster_id) + ", " + "\nAt: " + to_string(Cluster[Key].average.x) + ", " +
               to_string(Cluster[Key].average.y) + ", " + to_string(Cluster[Key].average.z) + "\nNormal vector: " +
               to_string(Cluster[Key].normal_vector(0)) + ", " + to_string(Cluster[Key].normal_vector(1)) + ", " +
               to_string(Cluster[Key].normal_vector(2)) + "\nActual points: " + to_string(Cluster[Key].data.size()) +
               ", " + "\nNumber of inliers: " + to_string(Cluster[Key].inliers) + ", " +
               "\nPercentages of inliers: " + to_string(Cluster[Key].percentages_inliers) + ", " +
               "\nBoundary points: " + to_string(Cluster[Key].boundary_pts) + ", " +
               "\nBoundary rings: " + to_string(Cluster[Key].boundary_rings) + ", " +
               "\nPayload points: " + to_string(Cluster[Key].payload_without_boundary) + ", " +
               "\nPose_xyz: " + to_string(Cluster[Key].pose_tag_to_lidar.translation[0]) + ", " +
               to_string(Cluster[Key].pose_tag_to_lidar.translation[1]) + ", " +
               to_string(Cluster[Key].pose_tag_to_lidar.translation[2]) +
               "\nPose_rpy: " + to_string(Cluster[Key].pose_tag_to_lidar.roll) + ", " +
               to_string(Cluster[Key].pose_tag_to_lidar.pitch) + ", " + to_string(Cluster[Key].pose_tag_to_lidar.yaw) +
               "\nIntensity: " + to_string(Cluster[Key].max_intensity.intensity) + ", " +
               to_string(Cluster[Key].min_intensity.intensity)));
    BoundMarkArray.markers.push_back(BoundaryMarker);

    LiDARTag::_assignVectorMarker(BoundaryMarker, visualization_msgs::Marker::ARROW,
                                  "NormalVector_z" + to_string(Cluster[Key].cluster_id), 0, 0, 1, 2, 0.01,
                                  Cluster[Key].principal_axes.col(2), Cluster[Key].average);
    BoundMarkArray.markers.push_back(BoundaryMarker);

    if (_id_decoding)
    {
      visualization_msgs::Marker IDMarker;
      LiDARTag::_assignMarker(IDMarker, visualization_msgs::Marker::TEXT_VIEW_FACING,
                              "Text" + to_string(Cluster[Key].cluster_id), 1, 1, 1, Cluster[Key].average, 1,
                              Cluster[Key].tag_size * 0.7, string(to_string(Cluster[Key].rkhs_decoding.id)));
      IDMarkArray.markers.push_back(IDMarker);
    }
    // Draw payload boundary marker
    visualization_msgs::Marker PayloadMarker;

    if (_adaptive_thresholding)
    {
      // Upper boundary
      for (int i = 0; i < Cluster[Key].tag_edges.upper_line.size(); ++i)
      {
        LiDARTag::_assignMarker(PayloadMarker, visualization_msgs::Marker::SPHERE,
                                "PayloadUpperBoundary_" + to_string(Cluster[Key].cluster_id), 0, 0, 1,
                                Cluster[Key].tag_edges.upper_line[i]->point, i, 0.015);
        PayloadMarkArray.markers.push_back(PayloadMarker);
      }

      // Lower boundary
      for (int i = 0; i < Cluster[Key].tag_edges.lower_line.size(); ++i)
      {
        LiDARTag::_assignMarker(PayloadMarker, visualization_msgs::Marker::SPHERE,
                                "PayloadLowerBoundary_" + to_string(Cluster[Key].cluster_id), 0, 0, 1,
                                Cluster[Key].tag_edges.lower_line[i]->point, i, 0.015);
        PayloadMarkArray.markers.push_back(PayloadMarker);
      }

      // Left boundary (green)
      for (int i = 0; i < Cluster[Key].tag_edges.left_line.size(); ++i)
      {
        LiDARTag::_assignMarker(PayloadMarker, visualization_msgs::Marker::SPHERE,
                                "PayloadLeftBoundary_" + to_string(Cluster[Key].cluster_id), 0, 1, 0,
                                Cluster[Key].tag_edges.left_line[i]->point, i, 0.015);
        PayloadMarkArray.markers.push_back(PayloadMarker);
      }

      // Right boundary (red)
      for (int i = 0; i < Cluster[Key].tag_edges.right_line.size(); ++i)
      {
        LiDARTag::_assignMarker(PayloadMarker, visualization_msgs::Marker::SPHERE,
                                "PayloadRightBoundary_" + to_string(Cluster[Key].cluster_id), 1, 0, 0,
                                Cluster[Key].tag_edges.right_line[i]->point, i, 0.015);
        PayloadMarkArray.markers.push_back(PayloadMarker);
      }
    }
    else
    {
      int count = 0;
      for (int i = 0; i < Cluster[Key].payload_boundary_ptr.size(); ++i)
      {
        LiDARTag::_assignMarker(PayloadMarker, visualization_msgs::Marker::SPHERE,
                                "PayloadBoundary_" + to_string(Cluster[Key].cluster_id), 1, 0, 0,
                                Cluster[Key].payload_boundary_ptr[i]->point, i, 0.015);
        PayloadMarkArray.markers.push_back(PayloadMarker);
        count++;
      }
    }

    // corner points and RANSAC line
    if (_adaptive_thresholding)
    {
      Eigen::Vector4f EigenPoint;
      PointXYZIRT point;  // just for conversion

      for (int i = 0; i < 4; ++i)  // 4 corners
      {
        pcl::lineWithLineIntersection(Cluster[Key].line_coeff[i], Cluster[Key].line_coeff[(i != 3) ? i + 1 : 0],
                                      EigenPoint, 1e-2);

        LiDARTag::_eigenVectorToPointXYZIRT(EigenPoint, point);
        LiDARTag::_assignMarker(PayloadMarker, visualization_msgs::Marker::SPHERE,
                                "Corner_" + to_string(Cluster[Key].cluster_id), 0, 1, 1, point, i, 0.02);
        PayloadMarkArray.markers.push_back(PayloadMarker);

        // RANSAC
        LiDARTag::_assignMarker(PayloadMarker, visualization_msgs::Marker::ARROW,
                                "RANSAC" + to_string(Cluster[Key].cluster_id), 1, 1, 0, point, i, 0.01);
        double Length = sqrt(pow(Cluster[Key].line_coeff[i][3], 2) + pow(Cluster[Key].line_coeff[i][4], 2) +
                             pow(Cluster[Key].line_coeff[i][5], 2));
        PayloadMarker.scale.x = 0.15;
        PayloadMarker.pose.orientation.x = Cluster[Key].line_coeff[i][3] / Length;
        PayloadMarker.pose.orientation.y = Cluster[Key].line_coeff[i][4] / Length;
        PayloadMarker.pose.orientation.z = Cluster[Key].line_coeff[i][5] / Length;
        PayloadMarkArray.markers.push_back(PayloadMarker);
      }
    }

    visualization_msgs::Marker Marker;
    for (int i = 0; i < Cluster[Key].data.size(); ++i)
    {
      if (Cluster[Key].data[i].valid != 1)
        continue;
      PointsInClusters++;
      OutCluster->push_back(Cluster[Key].data[i].point);

      double Intensity = Cluster[Key].data[i].point.intensity;
      LiDARTag::_assignMarker(Marker, visualization_msgs::Marker::SPHERE, to_string(Cluster[Key].cluster_id), Intensity,
                              Intensity, Intensity, Cluster[Key].data[i].point, i, 0.01);
      ClusterArray.markers.push_back(Marker);
    }
    if (_mark_cluster_validity)
    {
      for (int i = 0; i < Cluster[Key].edge_points.size(); ++i)
      {
        if (Cluster[Key].edge_points[i].valid != 1)
          continue;
        OutEdgeCluster->push_back(Cluster[Key].edge_points[i].point);
      }
      for (int i = 0; i < Cluster[Key].edge_group1.size(); ++i)
      {
        EdgeGroup1->push_back(Cluster[Key].edge_group1[i].point);
      }
      for (int i = 0; i < Cluster[Key].edge_group2.size(); ++i)
      {
        EdgeGroup2->push_back(Cluster[Key].edge_group2[i].point);
      }
      for (int i = 0; i < Cluster[Key].edge_group3.size(); ++i)
      {
        EdgeGroup3->push_back(Cluster[Key].edge_group3[i].point);
      }
      for (int i = 0; i < Cluster[Key].edge_group4.size(); ++i)
      {
        EdgeGroup4->push_back(Cluster[Key].edge_group4[i].point);
      }
      for (int ring = 0; ring < _beam_num; ++ring)
      {
        if (Cluster[Key].payload_right_boundary_ptr[ring] != 0)
        {
          BoundaryPts->push_back(Cluster[Key].payload_right_boundary_ptr[ring]->point);
        }
        if (Cluster[Key].payload_left_boundary_ptr[ring] != 0)
        {
          BoundaryPts->push_back(Cluster[Key].payload_left_boundary_ptr[ring]->point);
        }
      }
    }

    if (_id_decoding)
    {
      for (int i = 0; i < Cluster[Key].rkhs_decoding.associated_pattern_3d->cols(); ++i)
      {
        PointXYZIRT point;
        point.x = Cluster[Key].rkhs_decoding.associated_pattern_3d->col(i)(0);
        point.y = Cluster[Key].rkhs_decoding.associated_pattern_3d->col(i)(1);
        point.z = Cluster[Key].rkhs_decoding.associated_pattern_3d->col(i)(2);
        if (point.x >= 0)
        {
          point.intensity = 200;
        }
        else
        {
          point.intensity = 50;
        }
        OutPayload->push_back(point);
      }
      for (int i = 0; i < Cluster[Key].rkhs_decoding.template_points_3d.cols(); ++i)
      {
        PointXYZIRT point;
        point.x = Cluster[Key].rkhs_decoding.template_points_3d.col(i)(0);
        point.y = Cluster[Key].rkhs_decoding.template_points_3d.col(i)(1);
        point.z = Cluster[Key].rkhs_decoding.template_points_3d.col(i)(2);
        point.intensity = Cluster[Key].rkhs_decoding.template_points_3d.col(i)(3);
        OutPayload3D->push_back(point);
      }
    }
    if (_mark_cluster_validity)
    {
      for (int i = 0; i < Cluster[Key].rkhs_decoding.initial_template_points.cols(); ++i)
      {
        PointXYZIRT point;
        point.x = Cluster[Key].rkhs_decoding.initial_template_points.col(i)(0);
        point.y = Cluster[Key].rkhs_decoding.initial_template_points.col(i)(1);
        point.z = Cluster[Key].rkhs_decoding.initial_template_points.col(i)(2);
        point.intensity = Cluster[Key].rkhs_decoding.initial_template_points.col(i)(3);
        OutInitialTarget->push_back(point);
      }
      for (int i = 0; i < Cluster[Key].rkhs_decoding.template_points.cols(); ++i)
      {
        PointXYZIRT point;
        point.x = Cluster[Key].rkhs_decoding.template_points.col(i)(0);
        point.y = Cluster[Key].rkhs_decoding.template_points.col(i)(1);
        point.z = Cluster[Key].rkhs_decoding.template_points.col(i)(2);
        point.intensity = Cluster[Key].rkhs_decoding.template_points.col(i)(3);
        OutTarget->push_back(point);
      }
    }
    // Publish to a lidartag channel
    _detectionArrayPublisher(Cluster[Key]);
  }

  detectionsToPub.header = _point_cloud_header;
  detectionsToPub.frame_index = _point_cloud_header.seq;
  _boundary_marker_pub.publish(BoundMarkArray);
  _cluster_marker_pub.publish(ClusterArray);
  _payload_marker_pub.publish(PayloadMarkArray);
  _id_marker_pub.publish(IDMarkArray);
  _detectionArray_pub.publish(detectionsToPub);
}
}  // namespace BipedLab
